#include <stdio.h>
#include <acml.h>

#define MATRIX_IDX(n, i, j) j*n + i
#define MATRIX_ELEMENT(A, m, n, i, j) A[ MATRIX_IDX(m, i, j) ]

void init_matrix(double* A, int m, int n)
{
   double element = 0.0;
   for (int j = 0; j < n; j++)
   {
      for (int i = 0; i < m; i++)
      {
         MATRIX_ELEMENT(A, m, n, i, j) = 1.0 / (i + j + 1);
         element += 0.25;
      }
   }
}

void print_matrix(const double* A, int m, int n)
{
   for (int i = 0; i < m; i++)
   {
      for (int j = 0; j < n; j++)
      {
          printf("%8.4f", MATRIX_ELEMENT(A, m, n, i, j));
      }
      printf("\n");
   }
}

int main(int argc, char** argv)
{
   int m = 3;
   int n = 4;

   int p = (m < n ? m : n);

   double A[m * n];
   double U[m * m];
   double VT[n * n];
   double S[p * 1];

   int info;

   init_matrix(A, m, n);

   printf("Matrix A (%d x %d) is:\n", m, n);
   print_matrix(A, m, n);

   dgesvd('A', 'A', m, n, A, m, S, U, m, VT, n, &info);
   if (info != 0)
   {
      fprintf(stderr, "Warning: dgesvd returned with a non-zero status (info = %d)\n", info);
   }

   printf("\nMatrix U (%d x %d) is:\n", m, m);
   print_matrix(U, m, m);

   printf("\nVector S (%d x %d) is:\n", p, 1);
   print_matrix(S, p, 1);

   printf("\nMatrix VT (%d x %d) is:\n", n, n);
   print_matrix(VT, n, n);

   return 0;
}

