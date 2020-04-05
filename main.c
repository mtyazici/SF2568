//Ramiz Dundar
//Muhammet Tayyip Yazici
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#define MIN(a,b) ((a) < (b) ? (a) : (b))

//Necessary functions (go to code sectons for moe detailed explanation)
int comp (const void * elem1, const void * elem2);
void mergeArrays(double arr1[], double arr2[], int n1, 
                             int n2, double arr3[]);

//Simple function for get ceil of a double. ceill(2.4) = 3 for example
int ceil2(double num) {
    int inum = (int)num;
    if (num == (double)inum) {
        return inum;
    }
    return inum + 1;
}                            
int main(int argc, char **argv){
    //Initialization and necessary variables
    int P, myrank, N, I, i, evenprocess, step, n, rc;
    MPI_Status status;
    double start, end;
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &P); 
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
    /* Find problem size N from command line */ 
    if (argc < 2){
        printf("%d\n", argc);
        printf("No size N given");
        return 1;
    }
    N = atoi(argv[1]);
    /* Compute local indices for data distribution */
    int L = N/P;
    int R = N%P;
    I = (N+P-myrank-1)/P;      //(number of local elements)
    //int n = p*L+MIN(p,R)+i; //(global index for given (p,i)
    /* local size. Modify if P does not divide N */ 
    double x[I];
    double l[(N+P-myrank)/P];
    double r[(N+P-myrank-2)/P];
    double lz[I+(N+P-myrank)/P];
    double rz[I+(N+P-myrank-2)/P];
    start = MPI_Wtime();
    /* random number generator initialization */ 
    srandom(myrank+1);

    /* data generation */
    for (i = 0; i < I; i++)
        x[i] = ((double) random())/(INT_MIN);


    /* Initial local sort */
    qsort (x, sizeof(x)/sizeof(*x), sizeof(*x), comp);

    /* Global sorting phase*/
    for (step = 0; step < ceil2((double)P/2); step++) {

        /* Even phase */
        if(myrank%2 == 0 && myrank != P-1) {
            rc = MPI_Send(&x, I, MPI_DOUBLE, myrank+1, 0, MPI_COMM_WORLD);
            rc = MPI_Recv(&r, sizeof(r)/sizeof(r[0]) , MPI_DOUBLE, myrank+1, 0, MPI_COMM_WORLD, &status);
            /*Merge*/
            mergeArrays(x,r,I,sizeof(r)/sizeof(r[0]),rz);
            for( i = 0; i < I; i++){
                 x[i] = rz[i];
            }
               
        } else if(myrank%2 == 1 ){
            rc = MPI_Recv(&l, sizeof(l)/sizeof(l[0]) , MPI_DOUBLE, myrank-1, 0, MPI_COMM_WORLD, &status);
            rc = MPI_Send(&x, I, MPI_DOUBLE, myrank-1, 0, MPI_COMM_WORLD);
            /*Merge*/
            mergeArrays(x,l,I,sizeof(l)/sizeof(l[0]),lz);
            for( i = 0; i < I; i++){
                x[i] = lz[i+sizeof(l)/sizeof(l[0])];
            }
        }

        /* Odd phase */
         if(myrank%2 == 1 && myrank != P-1) {
            rc = MPI_Send(&x, I, MPI_DOUBLE, myrank+1, 0, MPI_COMM_WORLD);
            rc = MPI_Recv(&r, sizeof(r)/sizeof(r[0]) , MPI_DOUBLE, myrank+1, 0, MPI_COMM_WORLD, &status);
            /*Merge*/
            mergeArrays(x,r,I,sizeof(r)/sizeof(r[0]),rz);
            for( i = 0; i < I; i++)
                x[i] = rz[i];

        } else if(myrank%2 == 0 && myrank != 0){
            rc = MPI_Recv(&l, sizeof(l)/sizeof(l[0]) , MPI_DOUBLE, myrank-1, 0, MPI_COMM_WORLD, &status);
            rc = MPI_Send(&x, I, MPI_DOUBLE, myrank-1, 0, MPI_COMM_WORLD);
            /*Merge*/
            mergeArrays(x,l,I,sizeof(l)/sizeof(l[0]),lz);
            for( i = 0; i < I; i++)
                x[i] = lz[i+sizeof(l)/sizeof(l[0])];
        }

    }
    /* //Open this part for printing numbers
    for(i = 0; i < I; i++){
        printf("[%d] %e ",myrank*L+MIN(myrank,R)+i,x[i]);
    }
    printf("for p = %d\n",myrank);  
    */
        
    end = MPI_Wtime();

    //time calculation
    if(myrank == 0){
        printf("Execution time: %e\n",end-start);
    }
    MPI_Finalize();
    return 0;
}

//comparator for qsort used in initial sorting
int comp (const void * elem1, const void * elem2) 
{
    double f = *((double*)elem1);
    double s = *((double*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

//helper function for merging two sorted arrays into one sorted array
// original author: https://www.geeksforgeeks.org/merge-two-sorted-arrays/
void mergeArrays(double arr1[], double arr2[], int n1, 
                             int n2, double arr3[]) 
{ 
    int i = 0, j = 0, k = 0; 
  
    // Traverse both array 
    while (i<n1 && j <n2) 
    { 
        // Check if current element of first 
        // array is smaller than current element 
        // of second array. If yes, store first 
        // array element and increment first array 
        // index. Otherwise do same with second array 
        if (arr1[i] < arr2[j]) 
            arr3[k++] = arr1[i++]; 
        else
            arr3[k++] = arr2[j++]; 
    } 
  
    // Store remaining elements of first array 
    while (i < n1) 
        arr3[k++] = arr1[i++]; 
  
    // Store remaining elements of second array 
    while (j < n2) 
        arr3[k++] = arr2[j++]; 
} 