#include <cstdio>
#include <ctime>
#include <utility> //std::pair
#include <vector>
#include <cmath>

#include "matrix.h"

using namespace std;

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

static vector<pair<int, int>> swappedRows;

void buildX1(int n, double sol[]) {
  for (int i = 0; i < n; n++) {
    sol[i] = 1;
  }
}

void buildX2(int n, double sol[]) {
  for (int i = 0; i < n; n++) {
    sol[i] = 0.1;
  }
}

// alternating -1, +1
void buildX3(int n, double sol[]) {
  for (int i = 0; i < n; n++) {
    sol[i] = (i % 2) * 2 - 1;
  }
}

// alternating -5, +5
void buildX4(int n, double sol[]) {
  for (int i = 0; i < n; n++) {
    sol[i] = ((i % 2) * 2 - 1) * 5;
  }
}

// alternating -100, +100
void buildX5(int n, double sol[]) {
  for (int i = 0; i < n; n++) {
    sol[i] = ((i % 2) * 2 - 1) * 100;
  }
}

void buildX(int n, double *solutions[5]) {
  buildX1(n, solutions[0]);
  buildX2(n, solutions[1]);
  buildX3(n, solutions[2]);
  buildX4(n, solutions[3]);
  buildX5(n, solutions[4]);
}

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
  swappedRows.push_back(make_pair(i, pivotInd));
}

void luFactorization(int matrixSize) {
  for (int i = 0; i < matrixSize; i++) {
    double pivot = getValByCoords(i, i);
    if (pivot == 0) {
      // permute
      orderForPivot(i, matrixSize);
      pivot = getValByCoords(i, i);
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

void calcRealB(int n, double *x[5], double *b[5]) {
  for (int solInd = 0; solInd < 5; solInd++) {
    for (int row = 0; row < n; row++) {
      b[solInd][row] = 0;
      for (int j = 0; j < n; j++) {
        b[solInd][row] += getValByCoords(row, j) * x[solInd][row];
      }
    }
  }
}

// Output: x
void solveByLU(int n, double *b[5], double *x[5]) {
  for (int solInd = 0; solInd < 5; solInd++) {
    // permutate B's
    for (const auto &swappedEntry : swappedRows) {
      swap(b[solInd][swappedEntry.first], b[solInd][swappedEntry.second]);
    }
    // calc y
    double y[5][n];
    for (int i = 0; i < n; i++) {
      y[solInd][i] = b[solInd][i];
      for (int j = i - 1; j >= 0; j--) {
        y[solInd][i] -= getValByCoords(i, j) * y[solInd][j];
      }
    }
    // calc x
    for (int i = n - 1; i >= 0; i--) {
      x[solInd][i] = b[solInd][i];
      for (int j = i + 1; j < n; j++) {
        x[solInd][i] -= getValByCoords(i, j) * x[solInd][j];
      }
      x[solInd][i] /= getValByCoords(i, i);
    }
  }
}

void calcErrors(int n, double **realX, double **solvedX, double errors[]) {
  double delta[5], realNorm[5];
  for (int sol = 0; sol < 5; sol++) {
    delta[sol] = 0;
    realNorm[sol] = 0;
    for (int i = 0; i < n; i++) {
      delta[sol] += (solvedX[sol][i] - realX[sol][i]) * (solvedX[sol][i] - realX[sol][i]);
      realNorm[sol] += realX[sol][i] * realX[sol][i];
    }
    printf("%f %f\n", delta[sol], realNorm[sol]);
    if (realNorm[sol] == 0) {
      errors[sol] = delta[sol];
    } else {
      errors[sol] = sqrt(delta[sol] / realNorm[sol]);
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

  double **realX = new double*[5];
  for (int i = 0; i < 5; i++) { realX[i] = new double[n_cols]; }
  buildX(n_cols, realX);
  printf("Built real X\n");

  double **b = new double*[5];
  for (int i = 0; i < 5; i++) { b[i] = new double[n_cols]; }
  calcRealB(n_cols, realX, b);
  printf("Calculated real B\n");

  struct timespec start_time;
  clock_gettime(CLOCK_REALTIME, &start_time);

  /* Perform LU factorization here */
  luFactorization(n_rows);

  struct timespec mid_time, elapsed_time;
  clock_gettime(CLOCK_REALTIME, &mid_time);
  if (timespec_subtract(&elapsed_time, &mid_time, &start_time) != 0) {
    printf("Invalid elapsed time\n");
  };
  double elapsed = (double)elapsed_time.tv_sec +
      (double)elapsed_time.tv_nsec / 1000000000.0;
  fprintf(stderr, "LU factorization elapsed time: %f s\n", elapsed);

  // Solve system
  double **solvedX = new double*[5];
  for (int i = 0; i < 5; i++) { solvedX[i] = new double[n_cols]; }
  solveByLU(n_cols, b, solvedX);

  struct timespec end_time;
  clock_gettime(CLOCK_REALTIME, &end_time);

  // dump_nonzeros(n_rows, values, col_ind, row_ptr_begin, row_ptr_end);

  if (timespec_subtract(&elapsed_time, &end_time, &mid_time) != 0) {
    printf("Invalid elapsed time\n");
  }

  elapsed = (double)elapsed_time.tv_sec +
      (double)elapsed_time.tv_nsec / 1000000000.0;
  fprintf(stderr, "System solved in: %f s\n", elapsed);

  double errors[5];
  calcErrors(n_cols, realX, solvedX, errors);
  for (int i = 0; i < 5; i++) {
    printf("Error of solution %d: %f\n", i, errors[i]);
  }

  delete[] realX;
  delete[] b;
  delete[] solvedX;

  return 0;
}
