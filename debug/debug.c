
//Perform debugging here
#include<stdio.h>

void printResult(char* statement, double* C, int n) {
    printf("\n%s:", statement);
    for (int i=0; i<n; i++) {
        printf("\n");
        for (int j=0; j<n; j++) {
            printf("%f ", C[i * n + j]);
        }
    }
}

void dgemm0(double* A, double* B, double* C, int n) {
    for (int i=0; i<n; i++) {
        printf("\n");
        for (int j=0; j<n; j++) {
            for (int k=0; k<n; k++) {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
    printResult("Naive Solution", C, n);
}

void dgemm1(double* A, double* B, double* C, int n) {
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            register double sum = C[i * n + j];
            for (int k=0; k<n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
    printResult("Register Re-use Solution", C, n);
}

void dgemm2(double* A, double* B, double* C, int n) {
    
    for (int i=0; i<n; i+=2) {

        for (int j=0; j<n; j+=2) {
        
            register double C00 = C[i * n + j];
            register double C01 = C[i * n + j + 1];
            register double C10 = C[i * n + j + n];
            register double C11 = C[i * n + j + n + 1];

            for (int k=0; k<n; k+=2) {
        
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
    printResult("Aggressive Register Reuse Solution", C, n);
}

void dgemm3(double* A, double* B, double* C, int n) {
    
    for (int i=0; i<n; i+=3) {

        for (int j=0; j<n; j+=3) {
            
            register double C00 = C[i * n + j];
            register double C01 = C[i * n + j + 1];
            register double C02 = C[i * n + j + 2];
            register double C10 = C[i * n + j + n];
            register double C11 = C[i * n + j + n + 1];
            register double C12 = C[i * n + j + n + 2];
            register double C20 = C[i * n + j + 2 * n];
            register double C21 = C[i * n + j + 2 * n + 1];
            register double C22 = C[i * n + j + 2 * n + 2];

            for (int k=0; k<n; k+=3) {

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
    printResult("More Aggressive Register Reuse Solution", C, n);
}

void bijk(double *A, double *B, double *C, int n, int b) 
{
    int i, j, k, i1, j1, k1;
    for (i=0; i<n; i+=b) { 
        for (j=0; j<n; j+=b) { 
            for (k=0; k<n; k+=b) { 
                /* B x B mini matrix multiplications */
                for (i1=i; i1<i+b; i1++) { 
                    for (j1=j; j1<j+b; j1++) { 
                        register double sum = C[i1 * n + j1];
                        for (k1=k; k1<k+b; k1++) {
                            sum -= A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = sum;
                    }
                }
            }
        }
    }
    printResult("Block IJK Version", C, n);
}

void bjik(double *A, double *B, double *C, int n, int b) 
{
    int i, j, k, i1, j1, k1;
    for (j=0; j<n; j+=b) { 
        for (i=0; i<n; i+=b) { 
            for (k=0; k<n; k+=b) { 
                /* B x B mini matrix multiplications */
                for (i1=i; i1<i+b; i1++) { 
                    for (j1=j; j1<j+b; j1++) { 
                        register double sum = C[i1 * n + j1];
                        for (k1=k; k1<k+b; k1++) {
                            sum -= A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = sum;
                    }
                }
            }
        }
    }
    printResult("Block JIK Version", C, n);
}

void bkij(double *A, double *B, double *C, int n, int b) 
{
    int i, j, k, i1, j1, k1;
    for (k=0; k<n; k+=b) { 
        for (i=0; i<n; i+=b) { 
            for (j=0; j<n; j+=b) { 
                /* B x B mini matrix multiplications */
                for (i1=i; i1<i+b; i1++) { 
                    for (j1=j; j1<j+b; j1++) { 
                        register double sum = C[i1 * n + j1];
                        for (k1=k; k1<k+b; k1++) {
                            sum -= A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = sum;
                    }
                }
            }
        }
    }
    printResult("Block KIJ Version", C, n);
}

void bikj(double *A, double *B, double *C, int n, int b) 
{
    int i, j, k, i1, j1, k1;
    for (i=0; i<n; i+=b) { 
        for (k=0; k<n; k+=b) { 
            for (j=0; j<n; j+=b) { 
                /* B x B mini matrix multiplications */
                for (i1=i; i1<i+b; i1++) { 
                    for (j1=j; j1<j+b; j1++) { 
                        register double sum = C[i1 * n + j1];
                        for (k1=k; k1<k+b; k1++) {
                            sum -= A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = sum;
                    }
                }
            }
        }
    }
    printResult("Block IKJ Version", C, n);
}

void bjki(double *A, double *B, double *C, int n, int b) 
{
    int i, j, k, i1, j1, k1;
    for (j=0; j<n; j+=b) { 
        for (k=0; k<n; k+=b) { 
            for (i=0; i<n; i+=b) { 
                /* B x B mini matrix multiplications */
                for (i1=i; i1<i+b; i1++) { 
                    for (j1=j; j1<j+b; j1++) { 
                        register double sum = C[i1 * n + j1];
                        for (k1=k; k1<k+b; k1++) {
                            sum -= A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = sum;
                    }
                }
            }
        }
    }
    printResult("Block JKI Version", C, n);
}

void bkji(double *A, double *B, double *C, int n, int b) 
{
    int i, j, k, i1, j1, k1;
    for (k=0; k<n; k+=b) { 
        for (j=0; j<n; j+=b) { 
            for (i=0; i<n; i+=b) { 
                /* B x B mini matrix multiplications */
                for (i1=i; i1<i+b; i1++) { 
                    for (j1=j; j1<j+b; j1++) { 
                        register double sum = C[i1 * n + j1];
                        for (k1=k; k1<k+b; k1++) {
                            sum -= A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = sum;
                    }
                }
            }
        }
    }
    printResult("Block JKI Version", C, n);
}

void main() {
    double A[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36};
    double B[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36};
    double C[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    int n = 6; //Assuming all are square matrix
    
    //dgemm0(A, B, C, n);
    //dgemm1(A, B, C, n);
    //dgemm2(A, B, C, n);
    //dgemm3(A, B, C, n);

    //bijk(A, B, C, n, 2);
    //bjik(A, B, C, n, 2);
    //bkij(A, B, C, n, 2);
    //bikj(A, B, C, n, 2);
    //bjki(A, B, C, n, 2);
    //bkji(A, B, C, n, 2);

    int i, j, k, i1, j1, k1, b = 2;
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
    printResult("Works", C, n);
}
