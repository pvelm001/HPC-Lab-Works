#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 64;
  
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
     int idg, ir, ic, ik;
 
    //BLAS-2 GEPP
    for (idg=0; idg<n; idg++) { //idg - Diagnol Iterator - start to end

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
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, int n, int ib, int end, int b)
{
    int ir, ic, ik, ibr, ibc, ibk;

    //BLAS-3 Factorization - Register Reuse with Cache Blocking
    for (ik=ib; ik<end; ik+=b) { //ik - Dual Iterator - Lower Column and Upper Row Iterator
        for (ir=end; ir<n; ir+=b) {  //ir - Row Iterator - Main and Lower Rectangular Matrix
            for (ic=end; ic<n; ic+=b) { //ic - Column Iterator - Main and Upper Rectangular Matrix 
                
                //2 X 2 Block Matrix
                for (ibr=ir; ibr<ir+b; ibr+=2) { //ibr - 
                    int ibrr = ibr * n; //ibrr - ibr's Row 
                    for (ibc=ic; ibc<ic+b; ibc+=2) { 
                        register double C00 = A[ibrr + ibc];            
                        register double C01 = A[ibrr + ibc + 1];
                        register double C10 = A[ibrr + ibc + n];
                        register double C11 = A[ibrr + ibc + n + 1];
                        for (ibk=ik; ibk<ik+b; ibk+=2) {
                            int ibkr = ibk * n; //ibkr - ibk's Row 
                            register double A00 = A[ibrr + ibk];
                            register double A10 = A[ibrr + ibk + n]; 
                            register double B00 = A[ibkr + ibc];
                            register double B01 = A[ibkr + ibc + 1];

                            C00 -= A00 * B00; C01 -= A00 * B01;
                            C10 -= A10 * B00; C11 -= A10 * B01;

                            A00 = A[ibrr + ibk + 1];     //A01
                            A10 = A[ibrr + ibk + n + 1]; //A11
                            B00 = A[ibkr + ibc + n];     //B10
                            B01 = A[ibkr + ibc + 1 + n]; //B11

                            C00 -= A00 * B00; C01 -= A00 * B01;
                            C10 -= A10 * B00; C11 -= A10 * B01;
                        }
                        A[ibrr + ibc] = C00;
                        A[ibrr + ibc + 1] = C01;
                        A[ibrr + ibc + n] = C10;
                        A[ibrr + ibc + n + 1] = C11;
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

            //BLAS-2 Factorization
            for (ir=idg+1; ir<n; ir++) { //ir - Row Iterator - below diagnol to ends
                int irr = ir * n; //irr - ir's Row
                A[irr + idg] = A[irr + idg] / A[idgr + idg]; 
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

        //BLAS-3 Factorization - Register Reuse
        mydgemm(A, n, ib, end, b);
    }

    return 0;
}