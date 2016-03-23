# Makefile for Matrix Multiplication

.PHONY: clean
CC = gcc
CFLAGS = -Wall -c  
LFLAGS = -Wall   
OBJS = twoBody.o functions.o

twoBody: $(OBJS)
	$(CC) $(LFLAGS) -o $@ $^

$(OBJS): %.o: %.c
	$(CC) $(CFLAGS) $<

twoBody.o: functions.h

functions.o:
clean:
	rm -rf *.o matrixmult
