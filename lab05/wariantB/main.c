#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/tred2.c"
#include "numerical_recipes.c/nrutil.c"
#include "numerical_recipes.c/tqli.c"
#include "numerical_recipes.c/pythag.c"
#define N 7


float multiplyScalar(float ** x, float ** y) {
    float result = 0;
    for (int i = 1; i <=N; i++) {
        result += (*x)[i] * (*y)[i];
    }
    return result;
}

void Copy_matrix(float ***copy_to_matrix, float **copy_from_matrix)
{
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			(*copy_to_matrix)[i][j] = copy_from_matrix[i][j];
		}
	}
}

float norm(float **v) {
    return sqrt(multiplyScalar(v, v));
}

float ** multiplyMatrixByVector(float ** m, float * x, float ** result) {
    for (int i = 1; i <=N; i++) {
        result[i] = 0;
        for (int j = 1; j < N; j++) 
        (*result)[i] += m[i][j] * x[j];
    }
    return result;
}

void multiplyMatrixByMatrix(float ** m1, float ** m2, float ** res) {
    for (int i = 1; i <=N; i++) {
        for (int j = 1; j <=N; j++) {
            float resd = 0;
            for (int k =1; k <=N; k++) {
                resd += m1[i][k] * m2[k][j];
            }
            (res)[i][j] = resd;
        }
    }
}

void transposeMatrix(float *** m, float *** res) {
    for (int i = 1; i <=N; i++) {
        for (int j = 1; j <=N; j++) {
            (*res)[j][i] = (*m)[i][j];
        }
    }
}

void divideVectorByScalar(float ** x, float div, float ** res) {
    for (int i = 1; i <=N; i++) {
        (*res)[i] = (*x)[i] / div;
    }
}

void substractMatrixByScalar(float ***A, float scalar, float ***result) {
    for(int i = 1; i <=N; i++) {
        for (int j = 1; j <=N; j++) {
            (*result)[i][j] = (*A)[i][j] - scalar;
        }
    }
}

void fillMatrixColumn(float ***m, float **x, int k) {
    for (int i = 1; i <=N; i++) {
        (*m)[i][k] = (*x)[i];
    }
}

void refillWithTensorMultiply(float ***bef, float lambda, float **x, float ***res) {
    for (int i = 1; i <=N; i++) {
        for (int j = 1; j <=N; j++) {
            (*res)[i][j] = (*bef)[i][j] - (lambda * (*x)[i] * (*x)[j]);
        }
    }
}

float countElementForFill(int i, int j) {
    return (1.0) / (sqrt(2.0 + abs((i) - (j))));
}

void printVector(float ** x) {
    for (int i = 1; i <=N; i++) {
        printf("vector %d %f\n",i,(*x)[i]);
    }
}

void printArray(float **x) {
    for (int i = 1; i <=N; i++) {
        for (int j = 1; j <=N; j++) {
          printf("arr[%d][%d] = %f\n",i,j,x[i][j]);
        }
    }
}

void matrix_vector(float **matrix, float *_vector, float **result_vector)
{
	float suma;
	for (int i = 1; i <= N; i++)
	{
		suma = 0;
		for (int j = 1; j <= N; j++)
		{
			suma += matrix[i][j] * _vector[j];
		}
		(*result_vector)[i] = suma;
	}
}

int main() {
    const int IT_MAX = 12;

    float **A=matrix(1,N,1,N);

    for(int i = 1; i <=N; i++) {
        for(int j =1; j <=N; j++) {
            A[i][j] = countElementForFill(i, j);
        }
    }
   FILE * f = fopen("macierzA.dat", "w");
    float **w=matrix(1,N,1,N);  //skopiowac A
    Copy_matrix(&w,A);
    float **wn=matrix(1,N,1,N); 
    float **X=matrix(1,N,1,N);
    float **D=matrix(1,N,1,N);

 /*fprintf(f, "%s", "MACIERZ A");
    for(int i = 1; i <= N; i++) 
    {
      for(int j = 1; j <= N; j ++) 
      {
        fprintf(f, "%12g ", w[i][j]);
      }
      fprintf(f, "\n");
    }*/

   FILE * f1 = fopen("lambda.dat", "w");

    for(int k = 1; k <=N; k++) {
       float * xk0=vector(1,N); 
        for (int i=1;i<=N;i++){
          xk0[i]=1.;
        }
        float lambda = 0;
        for(int i = 0; i < IT_MAX; i++) {
            float * xn=vector(1,N); 
            matrix_vector(w, xk0, &xn);
            lambda = multiplyScalar(&xn, &xk0) / multiplyScalar(&xk0, &xk0);
            printf( "Lambda %d %8f\n",i,lambda);
            fprintf(f1, "%d %8g\n",i, lambda);
            divideVectorByScalar(&xn, norm(&xn), &xk0);
              
        }
        for (int i=1;i<=N;i++)
			        fprintf(f,"%d %8g\n",k,xk0[i]);
      fprintf(f,"\n");
        printf("\n");
        refillWithTensorMultiply(&w, lambda, &xk0, &wn);
        w = wn;
        fillMatrixColumn(&X, &xk0, k);
    }

    float **AX=matrix(1,N,1,N);
    float **XT=matrix(1,N,1,N);
    multiplyMatrixByMatrix(A, X, AX);
    transposeMatrix(&X, &XT);
    multiplyMatrixByMatrix(XT, AX, D);

    printf("Macierz D:");
    printArray(D);

    FILE * f3 = fopen("macierzD.dat", "w");

    fprintf(f3, "%s", "MACIERZ D");
    for(int i = 1; i <= N; i++) 
    {
      for(int j = 1; j <= N; j ++) 
      {
        fprintf(f3, "%12g ", D[i][j]);
      }
      fprintf(f3, "\n");
    }

    
  

    return 0;
}