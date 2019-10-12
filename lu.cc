#include <cstdio>
#include <ctime>

#include "matrix.h"

/* Code taken from the GLIBC manual.
 *
 * Subtract the ‘struct timespec’ values X and Y,
 * storing the result in RESULT.
 * Return 1 if the difference is negative, otherwise 0.
 */
static int
timespec_subtract (struct timespec *result,
                   struct timespec *x,
                   struct timespec *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_nsec < y->tv_nsec) {
    int nsec = (y->tv_nsec - x->tv_nsec) / 1000000000 + 1;
    y->tv_nsec -= 1000000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_nsec - y->tv_nsec > 1000000000) {
    int nsec = (x->tv_nsec - y->tv_nsec) / 1000000000;
    y->tv_nsec += 1000000000 * nsec;
    y->tv_sec -= nsec;
  }

  /* Compute the time remaining to wait.
     tv_nsec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_nsec = x->tv_nsec - y->tv_nsec;

  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}


/* Global variables holding the matrix data. To complete this assignment
 * you are requested to only use arrays and access these arrays with
 * subscripts. Do not use pointers.
 */

const int max_n_elements = 131072;
const int max_n_rows = 16384;

static double values[max_n_elements];

static int col_ind[max_n_elements];
static int row_ptr_begin[max_n_rows];
static int row_ptr_end[max_n_rows];

double getValByCoords(int i, int j) {
  for (int rowEntryIndex = row_ptr_begin[i]; rowEntryIndex <= row_ptr_end[i]; rowEntryIndex++) {
    if (col_ind[rowEntryIndex] == j) {
      return values[rowEntryIndex];
    }
  }
  // if coords aren't listed, it's a 0 entry
  return 0;
}

void setValByCoords(int i, int j, double newVal) {
  for (int rowEntryIndex = row_ptr_begin[i]; rowEntryIndex <= row_ptr_end[i]; rowEntryIndex++) {
    if (col_ind[rowEntryIndex] == j) {
      values[rowEntryIndex] = newVal;
      break;
    }
  }
}

void swapRows(int r1, int r2) {
  int temp_begin = row_ptr_begin[r1];
  int temp_end = row_ptr_end[r1];
  row_ptr_begin[r1] = row_ptr_begin[r2];
  row_ptr_end[r1] = row_ptr_end[r2];
  row_ptr_begin[r2] = temp_begin;
  row_ptr_end[r2] = temp_end;
}

// interchange rows: max pivot <-> i
void orderForPivot(int i, int matrixSize) {
  double pivot = getValByCoords(i, i), pivotInd = i;
  for (int line = i + 1; line < matrixSize; line++) {
    double currVal = getValByCoords(line, i);
    if (currVal > pivot) {
      pivot = currVal;
      pivotInd = line;
    }
  }

  swapRows(i, pivotInd);
}

void luFactorization(int matrixSize) {
  for (int i = 0; i < matrixSize; i++) {
    double pivot = getValByCoords(i, i);
    if (pivot == 0) {
      // permute
      orderForPivot(i, matrixSize);
      pivot = getValByCoords(i, i);
      // printf("ZERO\n");
    }
    for (int j = i + 1; j < matrixSize; j++) {
      double mult = getValByCoords(j, i) / pivot;
      setValByCoords(j, i, mult);
      for (int k = i + 1; k < matrixSize; k++) {
        double x = getValByCoords(j, k) - mult * getValByCoords(i, k);
        setValByCoords(j, k, x);
      }
    }
  }
}   

int main(int argc, char **argv) {
  if (argc != 2)
    {
      fprintf(stderr, "usage: %s <filename>\n", argv[0]);
      return -1;
    }

  int nnz, n_rows, n_cols;
  bool ok(false);

  ok = load_matrix_market(argv[1], max_n_elements, max_n_rows,
                          nnz, n_rows, n_cols,
                          values, col_ind, row_ptr_begin, row_ptr_end);
  if (!ok)
    {
      fprintf(stderr, "failed to load matrix.\n");
      return -1;
    }

  // dump_nonzeros(n_rows, values, col_ind, row_ptr_begin, row_ptr_end);


  struct timespec start_time;
  clock_gettime(CLOCK_REALTIME, &start_time);

  /* Perform LU factorization here */
  luFactorization(n_rows);
  // solve();

  struct timespec end_time;
  clock_gettime(CLOCK_REALTIME, &end_time);

  dump_nonzeros(n_rows, values, col_ind, row_ptr_begin, row_ptr_end);

  struct timespec elapsed_time;
  timespec_subtract(&elapsed_time, &end_time, &start_time);

  double elapsed = (double)elapsed_time.tv_sec +
      (double)elapsed_time.tv_nsec / 1000000000.0;
  fprintf(stderr, "elapsed time: %f s\n", elapsed);

  return 0;
}
