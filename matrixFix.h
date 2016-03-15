//
//  matrixFix.h
//  MatrixFix
//
//  Created by Seth on 2/2/16.
//  Copyright © 2016 Seth. All rights reserved.
//

#ifndef matrixFix_h
#define matrixFix_h

typedef struct{
    int rows;
    int cols;
    double ** data;
} matrix;

#define ELEM(mtx, row, col) mtx->data[row][col]

void matrixInit();
matrix * newMatrix(const int rows, const int cols);
void randomizeMatrix(const matrix * mtx);
void constMatrix(const matrix * mtx, double value);
void deleteMatrix(matrix * mtx);
void printMatrix(const matrix * mtx);
matrix * copyMatrix(const matrix * mtx);
int subtractMatrix(const matrix * mtx1, const matrix * mtx2);
int matrixProduct(const matrix * mtxA, const matrix * mtxB, matrix * mtxC);
int matrixProductFix1(const matrix * mtxA, const matrix * mtxB, matrix * mtxC);
int matrixProductFix2(const matrix * mtxA, const matrix * mtxB, matrix * mtxC);
int matrixProductFix3(const matrix * mtxA, const matrix * mtxB, matrix * mtxC);
int max( int a, int b, int c);
int matrixProductCacheObliv(const matrix * mtxA, const matrix * mtxB, matrix * mtxC, int startRA, int endRA, int startM, int endM, int startCB, int endCB);
void zeroMatrix(const matrix * mtx);

#endif /* matrixFix_h */
