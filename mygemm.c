#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

//Register Reuse part 1
void dgemm0(const double* A, const double* B, double* C, const int n)
{
    int i, j, k;
    for (i=0; i<n; i++) { //Iterator for C's Rows
        for (j=0; j<n; j++) { //Iterator for C's Columns
            for (k=0; k<n; k++) { //Iterator for A's Columns & B's Rows
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i=0; i<n; i++) { //Iterator for C's Rows
        for (j=0; j<n; j++) { //Iterator for C's Columns
            register double sum = C[i * n + j]; //Eliminates cache-miss for n reads
            for (k=0; k<n; k++) { //Iterator for A's Columns & B's Rows
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i=0; i<n; i+=2) { //Iterator for C's Rows

        for (j=0; j<n; j+=2) { //Iterator for C's Columns
        
            register double C00 = C[i * n + j]; //2 X 2 Block Matrix
            register double C01 = C[i * n + j + 1];
            register double C10 = C[i * n + j + n];
            register double C11 = C[i * n + j + n + 1];

            for (k=0; k<n; k+=2) { //Iterator for A's Columns & B's Rows
        
                register double A00 = A[i * n + k];
                register double A10 = A[i * n + k + n]; 
                register double B00 = B[k * n + j];
                register double B01 = B[k * n + j + 1];

                C00 += A00 * B00; C01 += A00 * B01;
                C10 += A10 * B00; C11 += A10 * B01;

                A00 = A[i * n + k + 1];     //A01
                A10 = A[i * n + k + n + 1]; //A11
                B00 = B[k * n + j + n];     //B10
                B01 = B[k * n + j + 1 + n]; //B11

                C00 += A00 * B00; C01 += A00 * B01;
                C10 += A10 * B00; C11 += A10 * B01;
            }

            C[i * n + j] = C00;
            C[i * n + j + 1] = C01;
            C[i * n + j + n] = C10;
            C[i * n + j + n + 1] = C11;
        }
    }
}
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i=0; i<n; i+=3) { //Iterator for C's Rows

        for (j=0; j<n; j+=3) { //Iterator for C's Columns
            
            register double C00 = C[i * n + j]; //3 X 3 Block Matrix
            register double C01 = C[i * n + j + 1];
            register double C02 = C[i * n + j + 2];
            register double C10 = C[i * n + j + n];
            register double C11 = C[i * n + j + n + 1];
            register double C12 = C[i * n + j + n + 2];
            register double C20 = C[i * n + j + 2 * n];
            register double C21 = C[i * n + j + 2 * n + 1];
            register double C22 = C[i * n + j + 2 * n + 2];

            for (k=0; k<n; k+=3) { //Iterator for A's Columns & B's Rows
                
                register double A00 = A[i * n + k];
                register double A10 = A[i * n + k + n]; 
                register double A20 = A[i * n + k + 2 * n];
                register double B00 = B[k * n + j];
                register double B01 = B[k * n + j + 1];
                register double B02 = B[k * n + j + 2];

                C00 += A00 * B00; C01 += A00 * B01; C02 += A00 * B02;
                C10 += A10 * B00; C11 += A10 * B01; C12 += A10 * B02;
                C20 += A20 * B00; C21 += A20 * B01; C22 += A20 * B02;

                A00 = A[i * n + k + 1];               //A01
                A10 = A[i * n + k + n + 1];           //A11
                A20 = A[i * n + k + 2 * n + 1];       //A21
                B00 = B[k * n + j + n];               //B10
                B01 = B[k * n + j + n + 1];           //B11
                B02 = B[k * n + j + n + 2];           //B12

                C00 += A00 * B00; C01 += A00 * B01; C02 += A00 * B02;
                C10 += A10 * B00; C11 += A10 * B01; C12 += A10 * B02;
                C20 += A20 * B00; C21 += A20 * B01; C22 += A20 * B02;

                A00 = A[i * n + k + 2];               //A02
                A10 = A[i * n + k + n + 2];           //A12
                A20 = A[i * n + k + 2 * n + 2];       //A22
                B00 = B[k * n + j + 2 * n];           //B20
                B01 = B[k * n + j + 2 * n + 1];       //B21
                B02 = B[k * n + j + 2 * n + 2];       //B22

                C00 += A00 * B00; C01 += A00 * B01; C02 += A00 * B02;
                C10 += A10 * B00; C11 += A10 * B01; C12 += A10 * B02;
                C20 += A20 * B00; C21 += A20 * B01; C22 += A20 * B02;
            }

            C[i * n + j] = C00;
            C[i * n + j + 1] = C01;
            C[i * n + j + 2] = C02;
            C[i * n + j + n] = C10;
            C[i * n + j + n + 1] = C11;
            C[i * n + j + n + 2] = C12;
            C[i * n + j + 2 * n] = C20;
            C[i * n + j + 2 * n + 1] = C21;
            C[i * n + j + 2 * n + 2] = C22;
        }
    }
}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n) 
{

}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void jik(const double *A, const double *B, double *C, const int n) 
{

}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void kij(const double *A, const double *B, double *C, const int n) 
{

}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{

}


void ikj(const double *A, const double *B, double *C, const int n) 
{

}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void jki(const double *A, const double *B, double *C, const int n) 
{

}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void kji(const double *A, const double *B, double *C, const int n) 
{

}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{

}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
    
}