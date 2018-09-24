#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#define ARRAY_SIZE 128
#define PROC_NUM 8
#define LOCAL_ARRAY_SIZE 16 //ARRAY_SIZE/PROC_NUM
//Randomized quickSort

void swap(int *a, int *b)
{
    int t;
    t = *a; *a = *b; *b = t;
}

int partition( int a[], int l, int r) {
    int pivot, i, j, t;
    if(l > r)
    {
        t = rand() % (l - r);
        pivot = a[t];
        swap(&a[t], &a[l]);
    }
    else
    {
        pivot = a[l];
    }
    i = l; j = r+1;
    while( 1)
    {
        do ++i; while( a[i] <= pivot && i <= r );
        do --j; while( a[j] > pivot );
        if( i >= j ) break;
        swap(&a[i], &a[j]);
        //t = a[i]; a[i] = a[j]; a[j] = t;
    }
    swap(&a[l], &a[j]);
    //t = a[l]; a[l] = a[j]; a[j] = t;
    return j;
}

void quickSort( int a[], int l, int r)
{
    int j;
    if( l < r ) 
    {
    // divide and conquer
        j = partition( a, l, r);
        quickSort( a, l, j-1);
        quickSort( a, j+1, r);
    }
}

void compare_Split_Hi(int *A, int *B, int *C, int size)
{
    int i, j, k;
    i = 0; j = 0; k = 0;
    for(i = 0; i < size; i++)
        C[i] = A[i];
    //First, merge those arrays into the larger Array D
    i = size - 1;
    j = size - 1;
    for(k = size - 1; k >= 0; k--)
    {
        if(j < 0 || ((i >= 0) && C[i] > B[j])) 
            A[k] = C[i--];
        else 
            A[k] = B[j--];
    }
}

void compare_Split_Lo(int *A, int *B, int *C, int size)
{
    int i, j, k;
    i = 0; j = 0; k = 0;
    for(i = 0; i < size; i++)
        C[i] = A[i];
    //First, merge those arrays into the larger Array D
    i = 0;
    for(k = 0; k < size; k++)
    {
        if(j == size || (i < size && C[i] <= B[j])) 
            A[k] = C[i++];
        else 
            A[k] = B[j++];
    }
}

void main(int argc, char** argv)
{
    int node;
    FILE* outFile;
    int nodeCnt;
    int i, j;
    time_t start_time, end_time;
    //Input and Output Matrices
    int A[ARRAY_SIZE];
    int partialMatA[LOCAL_ARRAY_SIZE];
    outFile = fopen("output.txt", "w");
    //MPI Initialize
    MPI_Init(&argc, &argv);
    
    //Get number of nodes using MPI_Comm_rank Function
    MPI_Comm_rank(MPI_COMM_WORLD, &node);

    char hostname[HOST_NAME_MAX];
    if (! gethostname(hostname, sizeof hostname) == 0)
        perror("gethostname");
        
    //Main Process
    //Initialize A and B Arrays and sends data to different processes
    if(node == 0)
    {
        //Initialization of Input Matrices
        for(i = 0; i < ARRAY_SIZE; i++) {
            A[i] = rand() % 128;
        }
        start_time = time(NULL);
        nodeCnt = 1;
        while(nodeCnt < PROC_NUM)
        {
            //Send Matrices A 16 ints to each node from node 0
            MPI_Send(&A[LOCAL_ARRAY_SIZE*nodeCnt], LOCAL_ARRAY_SIZE, MPI_INT, nodeCnt,0,MPI_COMM_WORLD);
            nodeCnt++;
        }
        for(i = 0; i < LOCAL_ARRAY_SIZE;i++)
        {
            partialMatA[i] = A[i];
        }
        quickSort(partialMatA, 0 ,LOCAL_ARRAY_SIZE - 1);
        //printf("Hello From Node : %d\n", node);
    }
    else
    {
        //Create Partial Matrices A consisting of 16 ints
        //Receive Matrices A from Node 0
        MPI_Recv(&partialMatA, LOCAL_ARRAY_SIZE, MPI_INT, 0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        quickSort(partialMatA, 0, LOCAL_ARRAY_SIZE-1);
        //printf("Hello From Node : %d\n", node);
    }
     
    //Phase 1: Hypercube Compare Exchange
    //Define D for hypercube
    int Dhyp = ceil(log(PROC_NUM)/log(2));
    int k, tmp;
    int distance;
    int tempMatB[LOCAL_ARRAY_SIZE], tempMatC[LOCAL_ARRAY_SIZE];
    for(k = Dhyp - 1; k >= 0; k--)
    {
        //Distance between the processes that are exchanging data
        //and doing compare-exchange in hypercube algorithm
        distance = 1 << k;
        //fprintf(outFile, "K: %d, Node: %d, Node & Distance: %d\n", k, node, node & distance);
        
        if(0 != (node & (1 << k)))
        {
            //printf("Sending Node: %d\tReceiving Node: %d\n",node,node ^ distance);
            MPI_Sendrecv(&partialMatA, LOCAL_ARRAY_SIZE, MPI_INT, node ^ distance, 0, &tempMatB, LOCAL_ARRAY_SIZE, MPI_INT, node ^ distance, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            compare_Split_Hi(&partialMatA[0], &tempMatB[0], &tempMatC[0], LOCAL_ARRAY_SIZE);
        }
        else
        {
            //printf("Sending Node: %d\tReceiving Node: %d\n",node,node ^ distance);
            MPI_Sendrecv(&partialMatA, LOCAL_ARRAY_SIZE, MPI_INT, node | distance, 0, &tempMatB, LOCAL_ARRAY_SIZE, MPI_INT, node | distance, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            compare_Split_Lo(&partialMatA[0], &tempMatB[0], &tempMatC[0], LOCAL_ARRAY_SIZE);
        }    
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //Phase 2: Odd-Even Transposition
    int phaseNum = 1;
    int broadCastFlag = 0;
    while(phaseNum < ARRAY_SIZE)
    {
        broadCastFlag = 0;
        if((phaseNum % 2) != 0)
        {
            if(((node+1) % 2) != 0)
            {
                if(node + 1 <= PROC_NUM)
                {
                    MPI_Sendrecv(&partialMatA, LOCAL_ARRAY_SIZE, MPI_INT, node + 1, 0, &tempMatB, LOCAL_ARRAY_SIZE, MPI_INT, node + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    compare_Split_Lo(partialMatA, tempMatB, tempMatC, LOCAL_ARRAY_SIZE);
                    for(k = 0; k < LOCAL_ARRAY_SIZE; k++)
                    {
                        if(partialMatA[i] != tempMatC[i])
                        {
                            broadCastFlag = 1;
                            break;
                        }
                    }       
                }
            }
            else
            {
                if(node - 1 >= 0)
                {
                    MPI_Sendrecv(&partialMatA, LOCAL_ARRAY_SIZE, MPI_INT, node - 1, 0, &tempMatB, LOCAL_ARRAY_SIZE, MPI_INT, node - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    compare_Split_Hi(partialMatA, tempMatB, tempMatC, LOCAL_ARRAY_SIZE);
                    for(k = 0; k < LOCAL_ARRAY_SIZE; k++)
                    {
                        if(partialMatA[i] != tempMatC[i])
                        {
                            broadCastFlag = 1;
                            break;
                        }
                    }       
                }
            }
        }
        else
        {
            if(((node+1) % 2) == 0)
            {
                if(node + 1 < PROC_NUM)
                {
                    MPI_Sendrecv(&partialMatA, LOCAL_ARRAY_SIZE, MPI_INT, node + 1, 0, &tempMatB, LOCAL_ARRAY_SIZE, MPI_INT, node + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    compare_Split_Lo(partialMatA, tempMatB, tempMatC, LOCAL_ARRAY_SIZE);
                    for(k = 0; k < LOCAL_ARRAY_SIZE; k++)
                    {
                        if(partialMatA[i] != tempMatC[i])
                        {
                            broadCastFlag = 1;
                            break;
                        }
                    }       
                }
            }
            else
            {
                if(node - 1 >= 0)
                {
                    MPI_Sendrecv(&partialMatA, LOCAL_ARRAY_SIZE, MPI_INT, node - 1, 0, &tempMatB, LOCAL_ARRAY_SIZE, MPI_INT, node - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    compare_Split_Hi(partialMatA, tempMatB, tempMatC, LOCAL_ARRAY_SIZE);
                    for(k = 0; k < LOCAL_ARRAY_SIZE; k++)
                    {
                        if(partialMatA[i] != tempMatC[i])
                        {
                            broadCastFlag = 1;
                            break;
                        }
                    }       
                }
            }
        }
        ++phaseNum;
        nodeCnt = 1;
        int tempFromNeighbor = 0;
        if(node == 0)
        {
            while(nodeCnt < PROC_NUM)
            {
                MPI_Recv(&tempFromNeighbor,1, MPI_INT, nodeCnt,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                broadCastFlag = broadCastFlag + tempFromNeighbor;
                nodeCnt++;   
            }
        }
        else
        {
        MPI_Send(&broadCastFlag,1, MPI_INT, 0,0,MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&broadCastFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        if(broadCastFlag == 0)
        {
            break; 
        }
    }
   
    //Final Phase: Gather all subarrays into the Process 0 and print the final result
    if(node == 0)
    {
        //Receive output sorted partial array A from other nodes
        for(i = 0; i < LOCAL_ARRAY_SIZE;i++)
        {
            A[i] = partialMatA[i];
        }
        nodeCnt = 1;
        while(nodeCnt < PROC_NUM)
        {
            MPI_Recv(&A[LOCAL_ARRAY_SIZE*nodeCnt], LOCAL_ARRAY_SIZE, MPI_INT, nodeCnt,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            nodeCnt++;
        }
        
        //Write Result Matrix
        end_time = time(NULL);
        printf("Execution Time: %d\n", end_time - start_time);
        
        for(i = 0; i < ARRAY_SIZE; i++)
        {
            printf("%d\t", A[i]);
            if((i+1) % LOCAL_ARRAY_SIZE == 0)
                printf("\n");
        }
        //printf("Proces %d: done\n", node);
    }
    else
    {   
        //int i, j;
        //Send Calculated matrix C back to node 0
        MPI_Send(&partialMatA, LOCAL_ARRAY_SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD);    
        //printf("Proces %d: done\n", node);
    }
    //End of MPI Communications
    MPI_Finalize();
}

