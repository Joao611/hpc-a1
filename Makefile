CC = gcc-4.9
CXX = g++-4.9
CFLAGS = -Wall -Wextra -g -std=c++11

all:	lu

lu:	lu.o mmio.o matrix.o
	$(CXX) $(CFLAGS) -o $@ $^ -lrt


mmio.o:	mmio.c
	$(CC) $(CFLAGS) -c mmio.c

%.o:	%.cc
	$(CXX) $(CFLAGS) -c $<

clean:
	rm -f *o
	rm -f lu
