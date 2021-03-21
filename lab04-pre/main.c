#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#define N 5

void fprintMatrix(FILE * file, char * txt, float ** matrix, int size) 
{
	fprintf(file, "%s", txt);
	for(int i = 1; i <= size; i++) 
  {
		for(int j = 1; j <= size; j ++) 
    {
			fprintf(file, "%12g ", matrix[i][j]);
		}
		fprintf(file, "\n");
  }
}

void fprintVector(FILE * file,char * txt1,char * txt2,float * v)
{
    fprintf(file, "%s", txt1);
    for(int i=1;i<=N;i++)
    {
        fprintf(file,txt2,i, v[i]);
    }
    //fprintf(file,"\n\n");
}

float scalar(float *a, float*b)
{
    float x=0.;

    for(int i=1;i<=N;i++)
        x+= a[i] * b[i];
    return x;
}

void transponeMX(float **A)
{
    float ** T =matrix(1,N,1,N);
    for(int i=1;i<=N;i++)
        for(int j=1;j<=N;j++)
        {
            T[i][j]=A[j][i];
        }
    for(int i=1;i<=N;i++)
        for(int j=1;j<=N;j++)
        {
            A[i][j]=T[i][j];
        }
}

void multiplyMatrix(float**A, float **B, float **C)
{
    for(int i=1;i<=N;i++)
    {
        for(int j=1;j<=N;j++)
        {
            C[i][j]=0.;
            for(int k=1;k<=N;k++)
                C[i][j] += A[i][k]*B[k][j];
        }
    }
}

int main(void) {
  float **A=matrix(1,N,1,N);

  //// 1. uzupełnianie macierzy A
  for(int i=1;i<=N;i++)
  {
    for(int j=1;j<=N;j++)
    {
      A[i][j]=sqrt(i+j);
    }
  }
  FILE *f=fopen("wyniki.txt","w");
  if (!f)
  {
    printf("ERR while opening file! Leaving...\n");
    return 0;
  }
  /// Wypisuje macierz A
  fprintMatrix(f,"Macierz A:\n",A, N);

  float *d=vector(1,N);
  float *e=vector(1,N);

  /// 2. Redukuje macierz A do postaci trójdiagonajnej 
  tred2(A,N,d,e);

  /// 3. Wypisuje macierz przekształcenia P
  fprintMatrix(f,"\n\nMacierz przekształcenia P:\n",A, N);

  float **Y=matrix(1,N,1,N);
  for(int i=1;i<=N;i++)
  {
    Y[i][i]=1;
  }

  /// 4. Szukam wartości i wektory własne T
  tqli(d,e,N,Y);

  /// 5. Wypisuje wartości własne T
  fprintVector(f,"\n\nWartości własne macierzy T:\n","lambda%d: %f\n",d);

  ///Wypisuje wektory własne
  fprintf(f,"\nWektory własne macierzy T:\n");
  for(int i=1;i<=N;i++)
  {     
      fprintf(f,"Wektor wlasny y%d: ", i);
      for(int j=1;j<=N-1;j++)
      {
          fprintf(f,"%g, ", Y[j][i]);
      }
      fprintf(f,"%g", Y[N][i]);
      fprintf(f,"\n");
  }

  /// 6. Szukam wektorów własnych A
  float **X=matrix(1,N,1,N);
  float suma;
  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j <= N; j++)
    {
      suma = 0;
      for (int k = 1; k <= N; k++)
      {
        suma += A[i][k] * Y[k][j];
      }
      X[i][j] = suma;
    }
  }
  
  /// Zapisuje wektory własne A
  fprintf(f,"\n\nWektory własne macierzy A:\n");
  for(int i=1;i<=N;i++)
  {
      fprintf(f,"Wektor wlasny x%d: ", i);
      for(int j=1;j<=N-1;j++)
      {
          fprintf(f,"%g, ", X[j][i]);
      }
      fprintf(f,"%g", X[N][i]);
      fprintf(f,"\n");
  }

  /// 7. Szukam Beta
  float *beta = vector(1,N);
  for(int i=1;i<=N;i++)
        for(int j=1;j<=N;j++)
            A[i][j] = (float)(sqrt(i+j));
  
  float **A2 = matrix(1,N,1,N);
  multiplyMatrix(A, X, A2);
  transponeMX(X);
  transponeMX(A2);

  /// 8. Zapisuje wartości Beta
  fprintf(f,"\n\nBeta:\n");
  for(int i=1;i<=N;i++)
    {
        beta[i] = scalar(X[i], A2[i])/scalar(X[i], X[i]);
        fprintf(f,"beta %d = %g \n", i, beta[i]);
    }

  fclose(f);
  free_matrix(A,1,N,1,N);
  free_matrix(Y,1,N,1,N);
  free_matrix(X,1,N,1,N);
  free_matrix(A2,1,N,1,N);
  free_vector(d,1,N);
  free_vector(e,1,N);
  free_vector(beta,1,N);
  return 0;
}
