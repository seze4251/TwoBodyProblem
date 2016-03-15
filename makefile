# Makefile for Matrix Multiplication

.PHONY: clean
CC = gcc
CFLAGS = -Wall -O3  -c  
LFLAGS = -Wall -O3  
OBJS = twoBody.o matrixFix.o

twoBody: $(OBJS)
	$(CC) $(LFLAGS) -o $@ $^

$(OBJS): %.o: %.c
	$(CC) $(CFLAGS) $<

twoBody.o:

matrixFix.o: matrixFix.h

clean:
	rm -rf *.o matrixmult
