#include<stdio.h>
#include<math.h>
#include<stdlib.h>

void printMatrix(double *A, int n) {
    printf("Matix A: \n");
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            printf("%f ", A[i * n + j]);
        }
        printf("\n");
    }
} 

void printPivot(int* ipiv, int n) {
    printf("\nPivot Information: \n");
    for (int i=0; i<n; i++) {
        printf("%d ", ipiv[i]);
    }
}

void printVector(double* B, int n) {
    printf("\nVector Information: \n");
    for (int i=0; i<n; i++) {
        printf("%f ", B[i]);
    }
}

int mydgetrf(double *A, int *ipiv, int n) 
{
    int idg, ir, ic, ik;
 
    //BLAS-2 GEPP
    for (idg=0; idg<n; idg++) { //idg - Diagnol Iterator - start to end

        int maxIndex = idg;
        int idgr = idg * n; //idgr - Diagnol Row
        double maxElement = fabs(A[idgr + idg]);

        //Find Maximum Row
        for (ir=idg+1; ir<n; ir++) { //ir - Row Iterator - below diagnol to end
            int irr = ir * n; //irr - ir's Row
            if (fabs(A[irr + idg]) > maxElement) {
                maxIndex = ir;
                maxElement = fabs(A[irr + idg]);
            }
        }

        if (maxElement == 0)  //Singluar Matrix Check
            return -1;
        
        else {

            //Partial Pivoting
            if (maxIndex != idg) {

                //Store Pivot Information
                int temps = ipiv[idg]; //Saving Pivoting Information
                ipiv[idg] = ipiv[maxIndex];
                ipiv[maxIndex] = temps;

                //Swap Rows
                double tempv; 
                for(ic=0; ic<n; ic++) { //ic - Column Iterator - start to end
                    tempv = A[idgr + ic];
                    A[idgr + ic] = A[maxIndex * n + ic];
                    A[maxIndex * n + ic] = tempv;
                }
            }
        }

        //BLAS-2 Factorization
        for (ir=idg+1; ir<n; ir++) { //ir - Row Iterator - below diagnol to ends
            int irr = ir * n; //irr - ir's Row
            A[irr + idg] = A[irr + idg] / A[idgr + idg]; 
            for (ik=idg+1; ik<n; ik++) { //ik - Dual Iterator - start to end
                A[irr + ik] -= A[irr + idg] * A[idgr + ik];
            }
        }
    }
    printMatrix(A, n);
    printPivot(ipiv, n);
    return 0;
}

void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    int ir, ik;
    double xy[n];
    
    //Forward Substitution
    if (UPLO == 'L') {                         
        xy[0] = B[ipiv[0]];
        for (ir=1; ir<n; ir++) { //ir - Row Iterator - start to end
            int irr = ir * n; //irr - ir's Row
            double sum = 0;
            for (ik=0; ik<ir; ik++) //ik - Dual Iterator
                sum += xy[ik] * A[irr + ik];
            xy[ir] = B[ipiv[ir]] - sum;
        }
    }

    //Backward Substitution
    else if (UPLO == 'U') {                      
        xy[n - 1] = B[n - 1] / A[(n * n) - 1];
        for (ir=n-2; ir>-1; ir--) { //ir - Row Iterator - end to start
            int irr = ir * n; //irr - ir's Row
            double sum = 0;
            for (ik=ir+1; ik<n; ik++) //ik - Dual Iterator
                sum += xy[ik] * A[irr + ik];
            xy[ir] = (B[ir] - sum) / A[irr + ir];
        }
    }

    for (ir=0; ir<n; ir++) { B[ir] = xy[ir]; }
    printVector(B, n);
    return;
}

int debug_mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    for (int ib=0; ib<n; ib+=b) {
        int end = ib + b;
        for (int row_iter=ib; row_iter<end; row_iter++) { //Performing BLAS 2 LU Factorization
            
            int maxIndex = row_iter;
            double maxElement = fabs(A[row_iter * n + row_iter]);

            for (int t=row_iter + 1; t<n; t++) {
                if (fabs(A[t * n + row_iter]) > maxElement) { //Finding Max Element
                    maxIndex = t;
                    maxElement = fabs(A[t * n + row_iter]);
                }
            }

            if (maxElement == 0) {
                return -1;
            }
            else {
                if (maxIndex != row_iter) {
                    int temps = ipiv[row_iter]; //Saving Pivoting Information
                    ipiv[row_iter] = ipiv[maxIndex];
                    ipiv[maxIndex] = temps;

                    double tempv = 0; //Swapping the rows
                    for(int iter=0; iter<n; iter++) {
                        tempv = A[row_iter * n + iter];
                        A[row_iter * n + iter] = A[maxIndex * n + iter];
                        A[maxIndex * n + iter] = tempv;
                    }

                }
            }
            //printf("Each Block Factorization:\n");
            //printMatrix(A, n);

            for (int col_iter=row_iter+1; col_iter<n; col_iter++) { //Performing Factorization
                A[col_iter * n + row_iter] = A[col_iter * n + row_iter] / A[row_iter * n + row_iter]; 
                for (int k=row_iter+1; k<end; k++) {
                    A[col_iter * n + k] = A[col_iter * n + k] - A[col_iter * n + row_iter] * A[row_iter * n + k];
                }
             }
        } 

        

        double* LL;
        LL = (double *) malloc(sizeof(double) * b * b);

        //Extracting strict lower triangular matrix
        for (int i=0; i<b; i++) {
            for (int j=0; j<b; j++) {
                if (i == j) {LL[i * b + j] = 1; } // LL + I
                else if (j < i) { LL[i * b + j] = -A[(ib + i) * n + ib + j]; } //Reading A's element is different from reading block element
                else { LL[i * b + j] = 0; }
            }
        }

        printf("\nPrinting Lower Triangular Matrix: ");
        for (int i=0; i<b; i++) {
                printf("\n");
            for (int j=0; j<b; j++) {
                printf("%f ", LL[i * b + j]);
            }
        }

        //Performing Matrix Multiplication
        for (int i=0; i<b; i++) {
            for (int j=end; j<n; j++) {
                double sum = 0;
                for (int k=0; k<b; k++) {
                    sum += LL[i * b + k] * A[(ib + k) * n + j];
                }
                A[(ib + i) * n + j] = sum;
            }
        } 
        
        /*
       //Performing Backward Substitution
        for (int i=ib; i<end; i++) {
            for (int j=end; j<n; j++) {
                double sum = 0;
                for (int k=0; k<i; k++) {
                    sum += A[i * n + k] * A[k * n + j];
                }
                A[i * n + j] -= sum;
            }
        }
        
        //BLAS-3 Factorization
        for (int i=end; i<n; i++) {
            for (int j=end; j<n; j++) {
                double sum = 0;
                for (int k=ib; k<ib+b; k++) {
                    sum += A[i * n + k] * A[k * n + j];
                }
                A[i * n + j] -= sum;
            }
        }*/

        //BLAS-3 Factorization
        int i, j, k, i1, j1, k1;
        for (k=ib; k<end; k+=b) { 
            for (i=end; i<n; i+=b) { 
                for (j=end; j<n; j+=b) { 
                    /* B x B mini matrix multiplications */
                    for (k1=k; k1<k+b; k1++) { 
                        for (i1=i; i1<i+b; i1++) { 
                            register double r = A[i1 * n + k1];
                            for (j1=j; j1<j+b;j1++) {
                                A[i1 * n + j1] -= r * A[k1 * n + j1];
                            }
                        }
                    }
                }
            }
        }

    }
    printMatrix(A, n);
    printPivot(ipiv, n);
    return 0;
}

void mydgemm(double *A, int n, int ib, int end, int b)
{
    int ir, ic, ik, ibr, ibc, ibk;

    //BLAS-3 Factorization - Register Reuse with Cache Blocking
    for (ik=ib; ik<end; ik+=b) { //ik - Dual Iterator - Lower Column and Upper Row Iterator
        for (ir=end; ir<n; ir+=b) {  //ir - Row Iterator - Main and Lower Rectangular Matrix
            for (ic=end; ic<n; ic+=b) { //ic - Column Iterator - Main and Upper Rectangular Matrix 
                
                //3 X 3 Register Block Matrix
                for (ibr=ir; ibr<ir+b; ibr+=3) { //ibr - 
                    int ibrr = ibr * n; //ibrr - ibr's Row 
                    for (ibc=ic; ibc<ic+b; ibc+=3) { 
                        register double C00 = A[ibrr + ibc];            
                        register double C01 = A[ibrr + ibc + 1];
                        register double C02 = A[ibrr + ibc + 2];
                        register double C10 = A[ibrr + ibc + n];
                        register double C11 = A[ibrr + ibc + n + 1];
                        register double C12 = A[ibrr + ibc + n + 2];
                        register double C20 = A[ibrr + ibc + 2 * n];
                        register double C21 = A[ibrr + ibc + 2 * n + 1];
                        register double C22 = A[ibrr + ibc + 2 * n + 2];

                        for (ibk=ik; ibk<ik+b; ibk+=3) {
                            int ibkr = ibk * n; //ibkr - ibk's Row 
                            register double A00 = A[ibrr + ibk];
                            register double A10 = A[ibrr + ibk + n]; 
                            register double A20 = A[ibrr + ibk + 2 * n];
                            register double B00 = A[ibkr + ibc];
                            register double B01 = A[ibkr + ibc + 1];
                            register double B02 = A[ibkr + ibc + 2];

                            C00 -= A00 * B00; C01 -= A00 * B01; C02 -= A00 * B02;
                            C10 -= A10 * B00; C11 -= A10 * B01; C12 -= A10 * B02;
                            C20 -= A20 * B00; C21 -= A20 * B01; C22 -= A20 * B02;

                            A00 = A[ibrr + ibk + 1];             //A01
                            A10 = A[ibrr + ibk + n + 1];         //A11
                            A20 = A[ibrr + ibk + 2 * n + 1];     //A21
                            B00 = A[ibkr + ibc + n];             //B10
                            B01 = A[ibkr + ibc + n + 1];         //B11
                            B02 = A[ibkr + ibc + n + 2];         //B12


                            C00 -= A00 * B00; C01 -= A00 * B01; C02 -= A00 * B02;
                            C10 -= A10 * B00; C11 -= A10 * B01; C12 -= A10 * B02;
                            C20 -= A20 * B00; C21 -= A20 * B01; C22 -= A20 * B02;

                            A00 = A[ibrr + ibk + 2];               //A02
                            A10 = A[ibrr + ibk + n + 2];           //A12
                            A20 = A[ibrr + ibk + 2 * n + 2];       //A22
                            B00 = A[ibkr + ibc + 2 * n];           //B20
                            B01 = A[ibkr + ibc + 2 * n + 1];       //B21
                            B02 = A[ibkr + ibc + 2 * n + 2];       //B22

                            C00 -= A00 * B00; C01 -= A00 * B01; C02 -= A00 * B02;
                            C10 -= A10 * B00; C11 -= A10 * B01; C12 -= A10 * B02;
                            C20 -= A20 * B00; C21 -= A20 * B01; C22 -= A20 * B02;


                        }
                        A[ibrr + ibc] = C00;
                        A[ibrr + ibc + 1] = C01;
                        A[ibrr + ibc + 2] = C02;
                        A[ibrr + ibc + n] = C10;
                        A[ibrr + ibc + n + 1] = C11;
                        A[ibrr + ibc + n + 2] = C12;
                        A[ibrr + ibc + 2 * n] = C20;
                        A[ibrr + ibc + 2 * n + 1] = C21;
                        A[ibrr + ibc + 2 * n + 2] = C22;
                    }
                }
            }
        }
    }    
    return;
}



int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    int ib, idg, ir, ic, ik;
    for (ib=0; ib<n; ib+=b) { //ib - Block Iterator - 0 to n
        int end = ib + b;
        
        //BLAS-2 GEPP
        for (idg=ib; idg<end; idg++) { //idg - Diagnol Iterator - block start to block end

            int maxIndex = idg; 
            int idgr = idg * n; //idgr - idg's Row
            double maxElement = fabs(A[idgr + idg]);

            //Find Maximum Row
            for (ir=idg+1; ir<n; ir++) { //ir - Row Iterator - below diagnol to end
                int irr = ir * n; //irr - ir's Row
                if (fabs(A[irr + idg]) > maxElement) { 
                    maxIndex = ir;
                    maxElement = fabs(A[irr + idg]);
                }
            }

            if (maxElement == 0) //Singluar Matrix Check
                return -1;
            
            else {

                //Partial Pivoting
                if (maxIndex != idg) {
                    
                    //Store Pivot Information
                    int temps = ipiv[idg]; 
                    ipiv[idg] = ipiv[maxIndex];
                    ipiv[maxIndex] = temps;

                    //Swap Rows - Register Reuse
                    double tempv = 0; 
                    for(ic=0; ic<n; ic++) { //ic - Column Iterator - start to end
                        tempv = A[idgr + ic];
                        A[idgr + ic] = A[maxIndex * n + ic];
                        A[maxIndex * n + ic] = tempv;
                    }

                }
            }

            //Split Loop - BLAS-1 Factorization
            for (ir=idg+1; ir<n; ir++) { //ir - Row Iterator - below diagnol to ends
                int irr = ir * n; //irr - ir's Row
                A[irr + idg] = A[irr + idg] / A[idgr + idg]; 
            }
            //BLAS-2 Factorization
            for (ir=idg+1; ir<n; ir++) { //ir - Row Iterator - below diagnol to ends
                int irr = ir * n; //irr - ir's Row
                for (ik=idg+1; ik<end; ik++) { //ik - Col Iterator - diagnol row start to diagnol column end
                    A[irr + ik] -= A[irr + idg] * A[idgr + ik];
                }
            }
        } 

       //Backward Substitution - Register Reuse
        for (ir=ib; ir<end; ir++) { //ir - Row Iterator (L) - block start to block end
            int irr = ir * n; //irr - ir's Row
            for (ic=end; ic<n; ic++) { //ic - Col Iterator (U) - block end to end
                double sum = 0;
                for (ik=ib; ik<ir; ik++) { //ik - Dual Iterator (L, U) - block start to diagnol & block start to block end
                    sum += A[irr + ik] * A[ik * n + ic];
                }
                A[irr + ic] -= sum;
            }
        }

        /*
        //BLAS-3 Factorization - Register Reuse
        for (ir=end; ir<n; ir++) { //ir - Row Iterator - block end to end
            for (ic=end; ic<n; ic++) { //ic - Col Iterator - block end to end
                int irr = ir * n; //irr - ir's Row
                double sum = 0;
                for (ik=ib; ik<end; ik++) { //ik - Dual Iterator - block start to block end
                    sum += A[irr + ik] * A[ik * n + ic];
                }
                A[irr + ic] -= sum;
            }
        }*/
        mydgemm(A, n, ib, end, b);
    }
    printMatrix(A, n);
    return 0;
}

void main() {

    //int n = 3;
    //double A[] = {1, -1, 3, 4, -2, 1, -3, -1, 4};
    //double B[] = {13, 15, 8};
    //2.000000 -2.000000 3.000000 

    //int n = 4;
    //double A[] = {2, 1, 1, 0, 4, 3, 3, 1, 8, 7, 9, 5, 6, 7, 9, 8};
    //double B[] = {4, 6, 8, -2};
    //2.000000 -1.000000 1.000000 -2.000000 

    //int n = 4;
    //double A[] = {1, 2, 1, 1, 0, 3, 4, 1, 0, 1, 5, 1, 0, 0, 1, 2};
    //double B[] = {2, 3, 1, 2};
    //-0.700000 0.900000 -0.200000 1.100000 

    //int n = 4;
    //double A[] = {9, 1, 1, 0, 4, 9, 3, 1, 8, 2, 9, 5, 6, 7, 1, 8};
    //double B[] = {4, 6, 8, -2};
    //0.297993 0.284208 1.033854 -0.851409 

    int n = 6;
    double A[] = {12, 2, 4, 8, 3, 4, -1, 9, -2, 4, 2, 1, 8, 2, 3, 4, 0, 5, 4, 5, 5, 9, 2, 2, 4, 3, 5, 7, 4, 2, 1, 1, 0, 1, 1, 5};
    double B[] = {5, 4, 10, -2, -5, 5};
    //0.974653 0.972739 -0.563630 -0.815842 -1.343102 1.042310 

    //int n = 8;
    //double A[] = {12, 2, 4, 8, 3, 4, 9, 10, -1, 9, -2, 4, 2, 1, 4, 5, 8, 2, 3, 4, 0, 5, 3, 4, 4, 5, 5, 9, 2, 2, 1, 2, 4, 3, 5, 7, 4, 2, 5, 6, 1, 1, 0, 1, 1, 5, 9, 10, 4, 0, 3, 5, 10, 2, 6, 7, 5, 2, 0, 5, 1, 2, 4, 10};
    //double B[] = {5, 4, 10, -2, -5, 5, 2, -8};
    //0.866299 0.274856 -2.174060 0.038294 0.753902 2.883760 1.648389 -2.618765 

    int ipiv[n];
    for(int i=0; i<n; i++) { ipiv[i] = i; }
    
    /*if (mydgetrf(A, ipiv, n) != -1) {
        mydtrsv('L', A, B, n, ipiv);
        mydtrsv('U', A, B, n, ipiv);
    }*/

    if (mydgetrf_block(A, ipiv, n, 3) != -1) {
        mydtrsv('L', A, B, n, ipiv);
        mydtrsv('U', A, B, n, ipiv);
    }

    /*if (debug_mydgetrf_block(A, ipiv, n, 2) != -1) {
        mydtrsv('L', A, B, n, ipiv);
        mydtrsv('U', A, B, n, ipiv);
    }*/


    //printMatrix(A, n);
    //mydgetrf(A, ipiv, n);
    //mydgetrf_block(A, ipiv, n, 2);
    //debug_mydgetrf_block(A, ipiv, n, 2);

}