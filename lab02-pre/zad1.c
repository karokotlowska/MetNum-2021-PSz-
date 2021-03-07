#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/nrutil.c"
#include "numerical_recipes.c/ludcmp.c"
#include "numerical_recipes.c/lubksb.c"
#define N 3

//Wypisywanie macierzy
void printMatrix(char * txt, float ** matrix, int size) {
	printf("%s", txt);
	for(int i = 1; i <= size; i++) {
		for(int j = 1; j <= size; j ++) {
			printf("%12g ", matrix[i][j]);
		}
		printf("\n");
	}
}
//Zapisywanie macierzy do pliku
void fprintMatrix(FILE * file, char * txt, float ** matrix, int size) {
	fprintf(file, "%s", txt);
	for(int i = 1; i <= size; i++) {
		for(int j = 1; j <= size; j ++) {
			fprintf(file, "%12g ", matrix[i][j]);
		}
		fprintf(file, "\n");
	}
}

void d_inicjalizuj_wektor(int **A)
{
    *A = ivector(1, N);
    for (int i = 1; i <= N; i++)
        (*A)[i] = 0;
}

float max_value(float **A)
{
    float max;
    if((A[1][1])>=0){
      max = (A[1][1]);
    }else{
      max = (-1)*(A[1][1]);
    }
    
    for (int i = 1; i <= N; i++){
        for (int j = 1; j <= N; j++){
            if ((A[i][j]) > max || (A[i][j])<(max*(-1))){
              if(A[i][j]>=0){
                max = A[i][j];
              }else max=(-1)*(A[i][j]);
            }
        }
    }
    //printf("%f\n",max);
    return max;
}

void pomnoz(float **A, float **B, float ***C)
{
    float suma;
    for (int i = 1; i <= N; i++)
    {
        for (int j = 1; j <= N; j++){
            suma = 0;
            for (int k = 1; k <= N; k++)
            {
                suma += A[i][k] * B[k][j];
            }
            (*C)[i][j] = suma;
        }
    }
}

///////////////////////////
int main(){

float **A=matrix(1,N,1,N);
float **B=matrix(1,N,1,N);

//uzupelniam macierze A i B
int k=1;
for (int i=1;i<=N;i++){
  for(int j = 1; j <= N; j ++) {
			A[i][j] = k;
			B[i][j] = k;
			k++;
		}
}
B[1][1]=1.1;
printMatrix("Macierz A:\n", A, N);
printMatrix("Macierz B:\n", B, N);

//tworze wektory premutacji typu int
int *indxA=ivector(1, N);
//d_inicjalizuj_wektor(&indxA);
int *indxB=ivector(1, N);
//d_inicjalizuj_wektor(&indxA);

//tworze zmienne typu float
float mA;
float mB;

//wykonuje procedurę ludcmp
ludcmp(A, N, indxA, &mA);
ludcmp(B, N, indxB, &mB);

//drukuje macierze po rozkładzie
printMatrix("Macierz A po rozkładzie:\n", A, N);
printMatrix("Macierz B op rozkładzie:\n", B, N);



//macierz L dla A trójkątna górna
float **LA=matrix(1, N, 1, N);

//macierz L dla B trójkątna górna
float ** LB = matrix(1, N, 1, N);

//macierz U dla A trójkątna dolna
float ** UA  = matrix(1, N, 1, N);

//macierz U dla B trójkątna dolna
float ** UB = matrix(1, N, 1, N);

//alokuje komorki maecierzy LA i LB
for (int i=1;i<=N;i++){
  for(int j = 1; j <= N; j ++) {
    if(i==j){
			LA[i][j] = 1;
			LB[i][j] = 1;
		}
    else if(i<j){
      LA[i][j] = 0;
			LB[i][j] = 0;
    }
    else{
      LA[i][j] = A[i][j];
			LB[i][j] = A[i][j];
    }
  }
}

//alokuje komorki maecierzy UA i UB
for (int i=1;i<=N;i++){
  for(int j = 1; j <= N; j ++) {
    if(i<=j){
      UA[i][j] = A[i][j];
			UB[i][j] = B[i][j];
    }
    else{
      UA[i][j] = 0;
			UB[i][j] = 0;
    }
  }
}

//drukuje macierze po rozkładzie
printMatrix("Macierz LA:\n", LA, N);
printMatrix("Macierz LB:\n", LB, N);
printMatrix("Macierz UA:\n", UA, N);
printMatrix("Macierz UB:\n", UB, N);

//tworze wektory wyrazów wolnych 
float *x1=vector(1, N);
x1[1]=1;x1[2]=0;x1[3]=0;
float *x2=vector(1, N);
x2[1]=0;x2[2]=1;x2[3]=0;
float *x3=vector(1, N);
x3[1]=0;x3[2]=0;x3[3]=1;

float *y1=vector(1, N);
y1[1]=1;y1[2]=0;y1[3]=0;
float *y2=vector(1, N);
y2[1]=0;y2[2]=1;y2[3]=0;
float *y3=vector(1, N);
y3[1]=0;y3[2]=0;y3[3]=1;
//
lubksb(A,N,indxA,x1);
lubksb(A,N,indxA,x2);
lubksb(A,N,indxA,x3);
lubksb(B,N,indxB,y1);
lubksb(B,N,indxB,y2);
lubksb(B,N,indxB,y3);



//towrze macierz A^-1 wynikowa
float **Awyn=matrix(1, N, 1, N);
//towrze macierz B^-1 wynikowa
float **Bwyn=matrix(1, N, 1, N);



//wpisuje wyniki
for(int i=1;i<=N;i++){
  Awyn[i][1]=x1[i];
  Awyn[i][2]=x2[i];
  Awyn[i][3]=x3[i];
  Bwyn[i][1]=y1[i];
  Bwyn[i][2]=y2[i];
  Bwyn[i][3]=y3[i];
}

//drukuje macierze po odwroceniu
printMatrix("Macierz A^-1':\n",Awyn, N);
printMatrix("Macierz B^-1':\n",Bwyn, N);

printf("maximum macierzy A %.2f\n",max_value(A));
printf("maximum macierzy A^-1 %.2f\n",max_value(Awyn));
printf("maximum macierzy B %.2f\n",max_value(B));
printf("maximum macierzy B^-1 %.2f\n",max_value(Bwyn));

//wskaznik uwarunkowania macierzy cond(A)=norm(A)*norm(A-1)
printf("\n\nwskaznik uwarunkowania macierzy %.2f\n",max_value(A)*max_value(Awyn));
printf("wskaznik uwarunkowania macierzy %.2f\n",max_value(B)*max_value(Bwyn));


k=1;
for (int i=1;i<=N;i++){
  for(int j = 1; j <= N; j ++) {
			A[i][j] = k;
			B[i][j] = k;
			k++;
		}
}
B[1][1]=1.1;

//ilocznyn A*A-1
float ** A_iloczyn = matrix(1, N, 1, N);
pomnoz(A,Awyn,&A_iloczyn);

//ilocznyn B*B-1
float ** B_iloczyn = matrix(1, N, 1, N);
pomnoz(B,Bwyn,&B_iloczyn);


//drukuje macierze po iloczynie
printMatrix("\n\nMacierz iloczynu AA^-1':\n",A_iloczyn, N);
printMatrix("Macierz iloczynu BB^-1':\n",B_iloczyn, N);


free_matrix(LA, 1, N, 1, N);
    free_matrix(UA, 1, N, 1, N);
    free_matrix(LB, 1, N, 1, N);
    free_matrix(UB, 1, N, 1, N);
    free_matrix(A, 1, N, 1, N);
    free_matrix(B, 1, N, 1, N);
    free_ivector(indxA, 1, N);
    free_ivector(indxB, 1, N);
    free_matrix(Awyn, 1, N, 1, N);
    free_matrix(Bwyn, 1, N, 1, N);
return 0;
}
