#include<iostream>
#include <set>
#include <stdlib.h>
#include<time.h>
#include <chrono>

using namespace std;
using std::set;

#include <chrono>
using namespace std::chrono;

bool PRINT_MATRIX = false;


////////// Storage of Sparse Matrices////////////////////
/* To calculate the product of 2 matrices A, B such that C = AB,
A has been stored in CSR format and B in CSC format.
This makes the multiplication easier since each row elements of A are multiplied with
the corresponding elements of each column of B*/


///////////////FUNCTION TO CREATE THE SPARSE MATRICES//////////////
void populateMatrix(int m, int n, int p, float sp1, float sp2, int* &Adata, int* &Acol, int* &Arowptr, int* &Bdata, int* &Brow, int* &Bcolptr)
{
    srand(time(0));
    
    int Amax = m * n * sp1, Bmax = n * p * sp2;
    int i, sum, temp, step;
    float frac;
 
    Arowptr = new int[m+1]();
    Bcolptr = new int[p+1]();
  
    /*Random selection of number of elements in each row of A*/
    sum = 0;
    step = Amax/m;
    step = step>1? step : 2;
    while(sum < Amax)
    {
        temp = rand() % step;
        i = 1 + rand() % m;
        if(Arowptr[i] + temp < n)
        {
           Arowptr[i] += temp ;
           sum += temp;
        }
    }
    for(i=1 ; i <= m ; i++)
    { 
        Arowptr[i] += Arowptr[i-1];
    }
    Amax = sum;

 
    /*Random selection of number of elements in each column of B*/
    sum = 0;
    step = Bmax/p;
    step = step>1? step : 2;
    while(sum < Bmax)
    {
        temp = rand() % step;
        i = 1 + rand() % p;
        if(Bcolptr[i] + temp < n)
        {
           Bcolptr[i] += temp ;
           sum += temp;
        }
    }
    for(i=1 ; i <= p ; i++)
    {
        Bcolptr[i] += Bcolptr[i-1];
    }
    Bmax = sum; 
 

    /*Populating the sparse data in A and B matrices*/
    Adata = new int[Amax];
    Bdata = new int[Bmax];
 
    for(i=0 ; i< Amax ; i++)
      Adata[i] = rand()%3 +1;
    for(i=0 ; i< Bmax ; i++)
      Bdata[i] = rand()%3 +1;


 
    Acol = new int[Amax];
    Brow = new int[Bmax];
 
    set<int> indices;
 
    /*Random column assignments to elements of each row of A*/
    sum = 0;
    for(i=1 ; i <= m ; i++)
    {
      temp = Arowptr[i] - Arowptr[i-1] ;
      while (indices.size() != temp)
      {
          indices.insert(rand() % n);
      }
  
      for(auto &x : indices)
      {
          Acol[sum] = x;
          sum++;
      } 
     
      indices.clear();
    }
  
 
    /*Random row assignments to elements of each column of B*/
    sum = 0;
    for(i=1 ; i <= p ; i++)
    {
      temp = Bcolptr[i] - Bcolptr[i-1] ;
      while (indices.size() != temp)
      {
          indices.insert(rand() % n);
      }
  
      for(auto &x : indices)
      {
          Brow[sum] = x;
          sum++;
      } 
     
      indices.clear();
    } 

 
     
}



void printMatrix(int **C, int m, int n)
{
    for(int i = 0; i < m ; i++)
         {
            for(int j = 0; j < n; j++)
               cout<<C[i][j]<<"\t";
            cout<<endl;
         }
}



///////////////////////////FUNCTION TO MULTIPLY CSR X CSC MATRICES ///////////////////////
/* C = AB , A in CSR format and B in CSC format */
 void sparseMatmult (int ** C, int m, int n, int p, int* Adata, int* Acol, int* Arowptr, int* Bdata, int* Brow, int* Bcolptr )
 {
     int i, j, k;
     int tempsum, num1, num2;
     int Aindex = 0, Bindex = 0;
  
     cout<<"\n CSR X CSC matmult \n";
     auto start = high_resolution_clock::now();
     for(i = 0; i < m ; i++)
     {
         for(j = 0; j < p; j++)
         {  
             tempsum = 0;
             num1 = Arowptr[i+1];
             num2 = Bcolptr[j+1];
             Aindex = Arowptr[i];
             Bindex = Bcolptr[j];             
             
             while (Aindex < num1 && Bindex < num2)
             {
                 if(Acol[Aindex] == Brow[Bindex])
                 {
                     tempsum += Adata[Aindex] * Bdata[Bindex] ;
                     Aindex++;
                     Bindex++;
                 }
                 else if (Acol[Aindex] > Brow[Bindex])
                 {
                     Bindex++;
                 }
                 else
                 {
                     Aindex++;
                 }
             }
             C[i][j] = tempsum;             
          
         }

     }
     auto stop = high_resolution_clock::now();
     auto duration = duration_cast<microseconds>(stop - start);
     cout << "Time taken :" << duration.count() << endl;
  
     if (PRINT_MATRIX)
         printMatrix(C, m, p);
     
              
         /*////////////////////////////////////////////////////////////////////
         The outermost loop runs for m iterations
         The middle loop runs for p iterations per each outermost iteration
         The innermost loops runs for an average of n * max(sp1 , sp2) times per middle loop iteration 
         {In an arbitrary row of A, there will be (n * sp1) elements on average. Similarly for each column of B}
         Thus, complexity = O(m * n * p * (sp1 + sp2))
         ////////////////////////////////////////////////////////////////////*/
         
 }


///////////////////////////FUNCTION TO MULTIPLY NORMAL MATRICES ///////////////////////
/* C = AB */
 void normalMatmult (int m, int n, int p, int* Adata, int* Acol, int* Arowptr, int* Bdata, int* Brow, int* Bcolptr )
 {
     int i, j, k, sum;
     int **A = new int*[m];
     int **B = new int*[n];
     int **C = new int*[m];
     for( i = 0 ; i < m ; i++)
        A[i] = new int[n]();
     for( i = 0 ; i < n ; i++)
        B[i] = new int[p]();  
     for(i = 0 ; i < m ; i++)
        C[i] = new int[p]();
     
     
     /*Conversion of A from CSR to normal matrix*/
     sum = 0;
     for(i = 0 ; i < m ; i++)
     {
         for( j = Arowptr[i] ; j < Arowptr[i + 1] ; j++)
         {
             A[i][Acol[j]] = Adata[j];
             sum++;
             
         }
     }
  
     /*Conversion of B from CSC to normal matrix*/
     sum = 0;
     for(i = 0 ; i < p ; i++)
     {
         for( j = Bcolptr[i] ; j < Bcolptr[i + 1] ; j++)
         {
             B[Brow[j]][i] = Bdata[j];
             sum++;
         }
     }


     cout<<"\n Normal matmult \n";
     auto start = high_resolution_clock::now();
     for( i = 0 ; i < m ; i++)
     {
         for( j = 0 ; j < p ; j++)
         {
             sum = 0;
             for ( k = 0 ; k < n ; k++)
                sum += A[i][k] * B[k][j];
             C[i][j] = sum;
         }
     }
     auto stop = high_resolution_clock::now();
     auto duration = duration_cast<microseconds>(stop - start);
     cout << "Time taken :" << duration.count() << endl;
  
     if (PRINT_MATRIX)
         printMatrix(C, m, p);
     
     
              
         /*////////////////////////////////////////////////////////////////////
         The outermost loop runs for m iterations
         The middle loop runs for p iterations per each outermost iteration
         The innermost loops runs for an average of n times per middle loop iteration 
         Thus, complexity = O(m * n * p )
         ////////////////////////////////////////////////////////////////////*/
  
 }

int main()
{
    int m,n,p;
    float sp1, sp2;
 
    m=100;
    n=20;
    p=100;
    sp1=0.12;
    sp2=0.25;
 
    int *Adata, *Acol, *Arowptr, *Bdata, *Brow, *Bcolptr;
 
    int **Res = new int*[m];
    for(int i = 0 ; i < m ; i++)
    {
        Res[i] = new int[p]();
    }
    populateMatrix(m, n, p, sp1, sp2, Adata, Acol, Arowptr, Bdata, Brow, Bcolptr);

    sparseMatmult(Res, m, n, p, Adata, Acol, Arowptr, Bdata, Brow, Bcolptr);
 
    normalMatmult(m, n, p, Adata, Acol, Arowptr, Bdata, Brow, Bcolptr);
    
    /*////////////////////////////////////////////////////////////////////
    Complexity of normal MatMult        = O(m * n * p )
    Complexity of sparse matrix MatMult = O(m * n * p * (sp1 + sp2))
    So for large enough dimensions and small enough sparsities, this sparse matrix MatMult is expected to have lesser time of execution
    ////////////////////////////////////////////////////////////////////*/
 
}
