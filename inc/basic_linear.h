#ifndef BASIC_LINEAR_H
#define BASIC_LINEAR_H

#include "stdlib.h"
#include "stdio.h"
#include "math.h"

//#define MIN 0.000001
#define MIN 0
#define TYPE double

#ifdef __cplusplus
 extern "C" {
#endif

 /* return a nr*nc TYPE matrix
  */
TYPE **dmatrix(int nr, int nc);

 /* free a m*n TYPE matrix
  */
void freedmat(TYPE** A, int m, int n);

/* return a n*n eye matrix
 */
TYPE **deyes(int n);

/* return a m*n zero matrix
 */
TYPE **dzero(int m, int n);

/* make a matrix of m*m become an eyes 
 */
void makeeyes(TYPE** A, int m);

/* empty every elements in a matrix of m*n
 */
void makezero(TYPE** A, int m, int n);

/* print a nr*nc TYPE matrix
 */
void printdmat(TYPE **A, int nr, int nc);

/* A+B=C
 */
void matadd(TYPE **A, TYPE **B, TYPE **C, int p, int q);

/* A-B=C
 */
void matsub(TYPE **A, TYPE **B, TYPE **C, int p, int q);

/* A*B=C
 * C cannot be B or A
 */
void matmult(TYPE **A, TYPE **B, TYPE **C, int p, int q, int r);

/* s*A=B; s is a scaler
*/
void matscal(TYPE **A, TYPE **B, TYPE s, int m, int n);

/* B(nc, nr) = transpose of A
 */
void transpose(TYPE **A, TYPE **B, int nr, int nc);

/* "Cholesky-Banachiewicz" method 
 * A should be positive-definite n*n matrix
 * the function should return L (A = L*L^T)
 */
void Cholesky_getL(TYPE **A, TYPE **L, int n);

/* Solve lower triangular matrix (system)
 * Lx[p] = b[p]
 */
void Solve_LowTri_Sys(TYPE **L, TYPE **x, TYPE **b, int n, int p);

/* Solve upper triangular matrix (system)
 * Ux[p] = b[p]
 */
void Solve_UpTri_Sys(TYPE **U, TYPE **x, TYPE **b, int n, int p);

/* use Cholesky B. method to do inverse calculation
 * input: [1] positive-definite matrix A with n*n dimension
 * [2] TYPE matrix CX with n*n dimension
 * [3] TYPE matrix Ainv with n*n dimension
 * [4] TYPE matrix L with n*n dimension
 * [5] TYPE matrix LT with n*n dimension
 * [6] TYPE eye matrix I with n*n dimension
 * A = L*LT
 * I = L*CX
 * X = LT*Ainv
 */
void Cholesky_Inv(TYPE **A, TYPE **CX, TYPE **Ainv, TYPE **L, TYPE **LT, TYPE **I, int n);


#ifdef __cplusplus
 }
#endif

#endif
