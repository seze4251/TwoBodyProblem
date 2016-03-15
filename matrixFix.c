//
//  matrixFix.c
//  MatrixFix
//
//  Created by Seth on 2/2/16.
//  Copyright © 2016 Seth. All rights reserved.
//

#include "matrixFix.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


void matrixInit() {
    srand((unsigned int) time(NULL));
}


// Flexible
matrix * newMatrix(const int rows, const int cols) {
    
    if ( (rows <= 0) || (cols <= 0) ) {
        return NULL;
    }
    
    matrix * mtx = (matrix *) malloc(sizeof(matrix));
    mtx -> rows = rows;
    mtx -> cols = cols;
    mtx -> data = (double **) calloc(rows, sizeof(double *));
    
    // The first position on A points to the start of the vector length row*col
    mtx -> data[0] = (double*) calloc(rows * cols, sizeof(double));
    
    // Sets every A value to be pointing to the 0th position of the next row
    int n;
    
    for ( n=1; n<rows; ++n) {
        mtx -> data[n] = mtx -> data[n-1]+ cols;
    }
    
    return mtx;
}


void randomizeMatrix(const matrix * mtx) {
    int i, j;
    
    for ( i = 0; i < mtx -> rows; i++) {
        for ( j = 0; j < mtx -> cols; j++) {
            ELEM(mtx,i,j) = rand();
        }
    }
}


void constMatrix(const matrix * mtx, const double value) {
    int i, j;
	if( !mtx) {
	printf("ConstantMatrix: mtx points to NULL");
	}
    for ( i = 0; i < mtx -> rows; i++) {
        for ( j = 0; j < mtx -> cols; j++) {
            ELEM(mtx,i,j) = value;
        }
    }
}


void zeroMatrix(const matrix * mtx) {
    int i, j;
    
    for ( i = 0; i < mtx -> rows; i++) {
        for ( j = 0; j < mtx -> cols; j++) {
            ELEM(mtx,i,j) = 0;
        }
    }
}




// Deletes the Matrix
void deleteMatrix(matrix * mtx) {
    free(mtx -> data[0]);
    free(mtx -> data);
    free(mtx);
}


// Prints the matrix
void printMatrix(const matrix * mtx) {
    printf("\n");
	if (!mtx) {
	printf("printMatrix: mtx points to null");
	}
    int i, j;
    for ( i = 0; i < mtx -> rows; i++) {
        for ( j = 0; j < mtx -> cols ; j++) {
            printf("%11.2f", ELEM(mtx,i,j));
        }
        printf("\n");
    }
}

// Copy a Matrix
matrix * copyMatrix(const matrix * mtx) {
    matrix * copy = newMatrix(mtx -> rows, mtx -> cols);
    int i, j;
    for ( i = 0; i < mtx -> rows; i++) {
        for ( j = 0; j < mtx -> cols; j++) {
            ELEM(copy,i,j) = ELEM(mtx,i,j);
        }
    }
    return copy;
}

// Subtract a Matrix
int subtractMatrix(const matrix * mtx1, const matrix * mtx2) {
    int i, j;
    if (!mtx1 || !mtx2) {
		printf("subtractMatrix: mtx1 or mtx2 point to null\n");
		}

    if ( mtx1 -> cols != mtx2 -> cols || mtx1 -> rows != mtx2 -> rows ) {
        printf("Mtx1 and Mtx2 must be the same size\n");
		return -2;
    }
    
    for ( i = 0; i < mtx1 -> rows; i++) {
        for ( j = 0; j < mtx2 -> cols; j++) {
            if (ELEM(mtx1,i,j) - ELEM(mtx2,i,j) != 0) {
                return -3;
            }
        }
    }
    return 0;
}


int matrixProduct(const matrix * mtxA, const matrix * mtxB, matrix * mtxC) {
    // If mtxA, mtxB, mtxC is NULL, return -1
    if (!mtxA || !mtxB) {
        return -1;
    }
    
    // If mtxA and mtxB are not the correct size, return -2
    if (mtxA -> cols != mtxB -> rows) {
        printf("Matrix A cols must equal Matrix B rows");
        return -2;
    }
    
    // If mtxC is not the correct size, return -3
    if (mtxC -> rows != mtxA -> rows || mtxC -> cols != mtxB -> cols ) {
        printf("Matrix C must be the correct size");
        return -3;
    }
    
    int i, j, k;
    
    // Perform the matrix multiplication
    for (i = 0; i < mtxA -> rows; i++) {
        for ( j = 0; j < mtxB -> cols ; j++) {
            for ( k = 0; k < mtxA -> cols ; k++) {
                ELEM(mtxC,i,j) += ELEM(mtxA,i,k) * ELEM(mtxB,k,j);
            }
        }
    }
    
    // Return True
    return 0;
    
}

int matrixProductFix1(const matrix * mtxA, const matrix * mtxB, matrix * mtxC) {
    // If mtxA, mtxB, mtxC is NULL, return -1
    if (!mtxA || !mtxB) {
        return -1;
    }
    
    // If mtxA and mtxB are not the correct size, return -2
    if (mtxA -> cols != mtxB -> rows) {
        printf("Matrix A cols must equal Matrix B rows");
        return -2;
    }
    
    // If mtxC is not the correct size, return -3
    if (mtxC -> rows != mtxA -> rows || mtxC -> cols != mtxB -> cols ) {
        printf("Matrix C must be the correct size");
        return -3;
    }
    
    double space00, space01, space10, space11;
    int i, j, k;
    
    // Perform the matrix multiplication
    for ( i = 0; i < mtxA -> rows; i+=2) {
        for ( j = 0; j < mtxB -> cols ; j+=2) {
            
            space00 = space01 = space10 = space11 = 0;
            
            for ( k = 0; k < mtxA -> cols ; k++) {
                
                space00 += ELEM(mtxA,i,k) * ELEM(mtxB,k,j);
                space01 += ELEM(mtxA,i,k) * ELEM(mtxB,k,j+1);
                space10 += ELEM(mtxA,i+1,k) * ELEM(mtxB,k,j);
                space11 += ELEM(mtxA,i+1,k) * ELEM(mtxB,k,j+1);
                
            }
            
            ELEM(mtxC,i,j) = space00;
            ELEM(mtxC,i,j+1) = space01;
            ELEM(mtxC,i+1,j) = space10;
            ELEM(mtxC,i+1,j+1) = space11;
        }
    }
    
    // Return the result
    return 0;
    
}

int min(int a, int b) {
    if ( a > b ) {
        return a;
    } else {
        return b;
    }
}

int matrixProductFix2(const matrix * mtxA, const matrix * mtxB, matrix * mtxC) {
    // If mtxA, mtxB, mtxC is NULL, return -1
    if (!mtxA || !mtxB) {
        return -1;
    }
    
    // If mtxA and mtxB are not the correct size, return -2
    if (mtxA -> cols != mtxB -> rows) {
        printf("Matrix A cols must equal Matrix B rows");
        return -2;
    }
    
    // If mtxC is not the correct size, return -3
    if (mtxC -> rows != mtxA -> rows || mtxC -> cols != mtxB -> cols ) {
        printf("Matrix C must be the correct size");
        return -3;
    }
    
    double space00, space01, space10, space11;
    int ib = 20;
    int ii, j, i, k;
    
    // Perform the matrix multiplication
    for ( ii = 0; ii < mtxA -> rows; ii+=ib) {
        int counter = (ii + ib); //min( (ii + ib), mtxA -> rows );
        for ( j = 0; j < mtxB -> cols; j+=2) {
            for ( i = ii; i < counter; i+=2) {
                
                space00 = space01 = space10 = space11 = 0;
                
                for ( k = 0; k < mtxA -> cols; k++) {
                    
                    space00 += ELEM(mtxA,i,k) * ELEM(mtxB,k,j);
                    space01 += ELEM(mtxA,i,k) * ELEM(mtxB,k,j+1);
                    space10 += ELEM(mtxA,i+1,k) * ELEM(mtxB,k,j);
                    space11 += ELEM(mtxA,i+1,k) * ELEM(mtxB,k,j+1);
                    
                }
                
                ELEM(mtxC,i,j) = space00;
                ELEM(mtxC,i,j+1) = space01;
                ELEM(mtxC,i+1,j) = space10;
                ELEM(mtxC,i+1,j+1) = space11;
            }
        }
    }
    
    // Return correct
    return 0;
    
}


int matrixProductFix3(const matrix * mtxA, const matrix * mtxB, matrix * mtxC) {
    // If mtxA, mtxB, mtxC is NULL, return -1
    if (!mtxA || !mtxB) {
        return -1;
    }
    
    // If mtxA and mtxB are not the correct size, return -2
    if (mtxA -> cols != mtxB -> rows) {
        printf("Matrix A cols must equal Matrix B rows");
        return -2;
    }
    
    // If mtxC is not the correct size, return -3
    if (mtxC -> rows != mtxA -> rows || mtxC -> cols != mtxB -> cols ) {
        printf("Matrix C must be the correct size");
        return -3;
    }
    
    double space00, space01, space10, space11;
    int ib = 20;
    int kb = 50;
    int ii, kk, j, i, k;
    
    // Perform the matrix multiplication
    for (ii = 0; ii < mtxA -> rows; ii+=ib) {
        for (kk = 0; kk < mtxA -> cols; kk += kb) {
            for (j = 0; j < mtxB -> cols; j+=2) {
                for (i = ii; i < (ii + ib); i+=2) {
                    if (kk == 0) {
                        space00 = space01 = space10 = space11 = 0;
                    } else {
                        
                        space00 = ELEM(mtxC,i,j);
                        space01 = ELEM(mtxC,i,j+1);
                        space10 = ELEM(mtxC,i+1,j);
                        space11 = ELEM(mtxC,i+1,j+1);
                    }
                    
                    for (k = kk; k < (kk + kb); k++) {
                        
                        space00 += ELEM(mtxA,i,k) * ELEM(mtxB,k,j);
                        space01 += ELEM(mtxA,i,k) * ELEM(mtxB,k,j+1);
                        space10 += ELEM(mtxA,i+1,k) * ELEM(mtxB,k,j);
                        space11 += ELEM(mtxA,i+1,k) * ELEM(mtxB,k,j+1);
                        
                    }
                    
                    ELEM(mtxC,i,j) = space00;
                    ELEM(mtxC,i,j+1) = space01;
                    ELEM(mtxC,i+1,j) = space10;
                    ELEM(mtxC,i+1,j+1) = space11;
                }
            }
        }
    }
    
    // Return Correct
    return 0;
    
}

int max( int a, int b, int c) {
    if (a > b && a > c) {
        return a;
    } else if ( b > c) {
        return b;
    } else {
        return c;
    }
}


// Recursive Basic Cache Oblivious Algorithim

int matrixProductCacheObliv(const matrix * mtxA, const matrix * mtxB, matrix * mtxC, int startRA, int endRA, int startM, int endM, int startCB, int endCB) {
    
    const int maxDim = 20;
    
    // If mtxA, mtxB, mtxC is NULL, return -1
    if (!mtxA || !mtxB) {
        return -1;
    }
    
    // If mtxA and mtxB are not the correct size, return -2
    if (mtxA -> cols != mtxB -> rows) {
        printf("Matrix A cols must equal Matrix B rows\n");
        return -2;
    }
    
    // If mtxC is not the correct size, return -3
    if (mtxC -> rows != mtxA -> rows || mtxC -> cols != mtxB -> cols ) {
        printf("Matrix C must be the correct size\n");
        return -3;
    }
    
    if ( ((endRA - startRA) < maxDim) && ((endM - startM) < maxDim) && ((endCB -startCB) < maxDim) ) {
        int i, j, k;
        
        // Perform the matrix multiplication
        for (i = startRA; i < endRA; i++) {
            for ( k = startM; k < endM ; k++) {
                for ( j = startCB; j < endCB ; j++) {
                    ELEM(mtxC,i,j) +=  ELEM(mtxA,i,k) * ELEM(mtxB,k,j);
                }
            }
        }
        
        return 0;
        
    } else {
        
        if ( max( (endRA - startRA), (endM - startM), (endCB -startCB) ) ==  (endRA - startRA)) {
            // Split A horizonally
            
            int c = startRA + (endRA - startRA) /2; // integers automatically floor in C
            matrixProductCacheObliv(mtxA, mtxB, mtxC, startRA, c, startM, endM, startCB, endCB);
            matrixProductCacheObliv(mtxA, mtxB, mtxC, c, endRA, startM, endM, startCB, endCB);
            return 0;
            
        } else if (max( (endRA - startRA), (endM - startM), (endCB -startCB) ) ==  (endCB -startCB)) {
            // Split B Vertically
            
            int c = startCB + (endCB -startCB) /2;
            matrixProductCacheObliv(mtxA, mtxB, mtxC, startRA, endRA, startM, endM, startCB, c);
            matrixProductCacheObliv(mtxA, mtxB, mtxC, startRA, endRA, startM, endM, c, endCB);
            return 0;
            
        } else {
            
            // Split the Columns in A and the Rows in B
            int c = startM + (endM - startM) /2;
            matrixProductCacheObliv(mtxA, mtxB, mtxC, startRA, endRA, startM, c, startCB, endCB);
            matrixProductCacheObliv(mtxA, mtxB, mtxC, startRA, endRA, c, endM, startCB, endCB);
            return 0;
        }
    }
}


