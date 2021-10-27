#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 128;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    int i, t, j, k, iter;
    for (i=0; i<n; i++) { //Performing LU Factorization
        
        int maxIndex = i;
        double maxElement = fabs(A[i * n + i]);

        for (t=i+1; t<n; t++) { //Finding maximum element 
            if (fabs(A[t * n + i]) > maxElement) {
                maxIndex = t;
                maxElement = fabs(A[t * n + i]);
            }
        }

        if (maxElement == 0) {
            return -1;
        }
        else {
            if (maxIndex != i) {
                int temps = ipiv[i]; //Saving Pivoting Information
                ipiv[i] = ipiv[maxIndex];
                ipiv[maxIndex] = temps;

                double tempv[n]; //Swapping the rows
                for(iter=0; iter<n; iter++) {
                    tempv[iter] = A[i * n + iter];
                    A[i * n + iter] = A[maxIndex * n + iter];
                    A[maxIndex * n + iter] = tempv[iter];
                }
            }
        }

        for (j=i+1; j<n; j++) { //Performing Factorization
            A[j * n + i] = A[j * n + i] / A[i * n + i]; 
            for (k=i+1; k<n; k++) {
                A[j * n + k] = A[j * n + k] - A[j * n + i] * A[i * n + k];
            }
        }
    }

    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    int i, k;
    double xy[n];
    for (i=0; i<n; i++) {
        xy[i] = 0;
    }

    if (UPLO == 'L') {                         //Forward Substitution
        xy[0] = B[ipiv[0]];
        for (i=1; i<n; i++) {
            double sum = 0;
            for (k=0; k<i; k++) {
                sum += xy[k] * A[i * n + k];
            }
            xy[i] = B[ipiv[i]] - sum;
        }
    }

    else if (UPLO = 'U') {                      //Backward Substitution
        xy[n - 1] = B[n - 1] / A[(n * n) - 1];
        for (i=n-2; i>-1; i--) {
            double sum = 0;
            for (k=i+1; k<n; k++) {
                sum += xy[k] * A[i * n + k];
            }
            xy[i] = (B[i] - sum) / A[i * n + i];
        }
    }

    for (i=0; i<n; i++) 
        B[i] = xy[i];

    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    int i1, j1, k1;
    for (k=0; k<n; k+=b) { 
        for (i=0; i<n; i+=b) { 
            for (j=0; j<n; j+=b) { 
                /* B x B mini matrix multiplications */
                for (i1=i; i1<i+b; i1+=2) { 
                    for (j1=j; j1<j+b; j1+=2) { 
                        register double C00 = C[i1 * n + j1]; //2 X 2 Block Matrix
                        register double C01 = C[i1 * n + j1 + 1];
                        register double C10 = C[i1 * n + j1 + n];
                        register double C11 = C[i1 * n + j1 + n + 1];
                        for (k1=k; k1<k+b; k1+=2) {
                            register double A00 = A[i1 * n + k1];
                            register double A10 = A[i1 * n + k1 + n]; 
                            register double B00 = B[k1 * n + j1];
                            register double B01 = B[k1 * n + j1 + 1];

                            C00 += A00 * B00; C01 += A00 * B01;
                            C10 += A10 * B00; C11 += A10 * B01;

                            A00 = A[i1 * n + k1 + 1];     //A01
                            A10 = A[i1 * n + k1 + n + 1]; //A11
                            B00 = B[k1 * n + j1 + n];     //B10
                            B01 = B[k1 * n + j1 + 1 + n]; //B11

                            C00 += A00 * B00; C01 += A00 * B01;
                            C10 += A10 * B00; C11 += A10 * B01;
                        }
                        C[i1 * n + j1] = C00;
                        C[i1 * n + j1 + 1] = C01;
                        C[i1 * n + j1 + n] = C10;
                        C[i1 * n + j1 + n + 1] = C11;
                    }
                }
            }
        }
    }
    return;
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{

    return 0;
}

