#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/nrutil.c"
#include "numerical_recipes.c/ludcmp.c"
#include "numerical_recipes.c/lubksb.c"
#define N 4


void printM(float **M){
  printf("\n");
  for(int i=1;i<=N;i++){
    for(int j=1;j<=N;j++){
      printf("%f\t",M[i][j]);
    }
    printf("\n");
  }
}

void fprint(float **M, FILE *out){
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            fprintf(out, "%f ", M[i][j]);
        }
        fprintf(out, "\n");
    }
}

void mnozenie(float **A,float **B, float **C){
  for(int i=1;i<=N;i++){
    for (int j=1;j<=N;j++){
      C[i][j]=0;
      for(int k=1;k<=N;k++){
        C[i][j]+=A[i][k]*B[k][j];
      }
    }
  }
}

float maximum(float **A){
  float max;
  if((A[1][1])>=0){
    max=A[1][1];
  }
  else{
    max=(-1)*A[1][1];
  }
  for(int i=1;i<=N;i++){
    for(int j=1;j<=N;j++){
      if(A[i][j]>max || A[i][j]<max*(-1)){
        if(A[i][j]>=0){
          max=A[i][j];
        }
        else{
          max=A[i][j]*(-1);
        }
      }
    }
  }
  return max;
}



int main(void) {

//otwieram plik
FILE *out=fopen("wynik.txt","w+");
if(out==NULL){
  return 1;
}

//towrzenie macierzy kwadratowej
float **A=matrix(1,N,1,N);

//uzupelnianie macierzy A
for(int i=1;i<=N;i++){
  for(int j=1;j<=N;j++){
    A[i][j]=1.0/(i+j);
  }
}
//drukuję macierz A
printf("\nMacierz A:");
printM(A);
fprintf(out,"\nMacierz A:\n");
fprint(A,out);

//----------------------

//tworzę wektor permutacji potrzebny do wywolania ludcmp
int *indxA=ivector(1,N);
//tworzę zmienną typu float potrzebna do wywolania ludcmp
float d;
//wywoluje ludcmp
ludcmp(A,N,indxA,&d);

//nadpisaną macierz A rozdzialam na macierze L i U
float **L=matrix(1,N,1,N);
float **U=matrix(1,N,1,N);
//petla dla macierzy L
for(int i=1;i<=N;i++){
  for (int j=1;j<=N;j++){
    if(i==j){
      L[i][j]=1;
    }
    else if(i<j){
      L[i][j]=0;
    }
    else{
      L[i][j]=A[i][j];
    }
  }
}
//petla dla macierzy U
for (int i=1;i<=N;i++){
  for (int j=1;j<=N;j++){
    if(i<=j){
      U[i][j]=A[i][j];
    }
    else{
      U[i][j]=0;
    }
  }
}
//drukuje obie macierze:
printf("\nMacierz L:");
printM(L);
fprintf(out,"\nMacierz L:\n");
fprint(L,out);
printf("\nMacierz U:");
printM(U);
fprintf(out,"\nMacierz U:\n");
fprint(U,out);

//zapisuje diagonalne elementy macierzy do pliku
fprintf(out,"\nElementy diagonalne:\n");
for(int i=1;i<=N;i++){
  fprintf(out,"\t%f\n",U[i][i]);
}

//zapisuje do pliku wyznacznik
float det=1;
for(int i=1;i<=N;i++){
    det=det*U[i][i];
}
printf("\n%e\n",det);
fprintf(out,"\nWyznacznik:\n");
fprintf(out,"%e",det);

//---------------------
//szukam macierzy odrwotnej
//tworzę wektory wyrazów wolnych
float *b1=vector(1,N);
b1[1]=1;b1[2]=0;b1[3]=0;b1[4]=0;
float *b2=vector(1,N);
b2[1]=0;b2[2]=1;b2[3]=0;b2[4]=0;
float *b3=vector(1,N);
b3[1]=0;b3[2]=0;b3[3]=1;b3[4]=0;
float *b4=vector(1,N);
b4[1]=0;b4[2]=0;b4[3]=0;b4[4]=1;

//wywoluje metode lubksb
lubksb(A,N,indxA,b1);
lubksb(A,N,indxA,b2);
lubksb(A,N,indxA,b3);
lubksb(A,N,indxA,b4);

//tworze macierz A-1 - wynikową
float **Awyn=matrix(1,N,1,N);

//wpisuje do maceirzy wyniki
for(int i=1;i<=N;i++){
  Awyn[i][1]=b1[i];
  Awyn[i][2]=b2[i];
  Awyn[i][3]=b3[i];
  Awyn[i][4]=b4[i];
}
//drukuje macierz
printf("\nMacierz A-1:");
printM(Awyn);
fprintf(out,"\n\nMacierz A-1:\n");
fprint(A,out);

//--------------------
//obliczam iloczyn AA-1
//tworze macierz w ktorej zapisze wynik
float **A_iloczyn=matrix(1,N,1,N);
//wracam do macierzy A, ktora zostala wczesniej nadpisana
for(int i=1;i<=N;i++) 
	{
		for(int j=1;j<=N;j++)
		{
			A[i][j]=1.0/(i+j);				
		}
	}
mnozenie(A,Awyn,A_iloczyn);
printf("\nMacierz AA-1:");
printM(A_iloczyn);
fprintf(out,"\nMacierz AA-1:\n");
fprint(A_iloczyn,out);

//liczę wskaznik uwarunkowania cond(A)=norm(A)*norm(A-1)
printf("\nMaximum macierzy A: \t");
printf("%f",maximum(A));
fprintf(out,"\nMaximum macierzy A: \t");
fprintf(out,"%f",maximum(A));

printf("\nMaximum macierzy A-1: \t");
printf("%f",maximum(Awyn));
fprintf(out,"\n\nMaximum macierzy A-1: \t");
fprintf(out,"%f",maximum(Awyn));

printf("\nWskaznik uwarunkowania macierzy:\n");
printf("%f",maximum(A)*maximum(Awyn));
fprintf(out,"\nWskaznik uwarunkowania macierzy:\n");
fprintf(out,"%f",maximum(A)*maximum(Awyn));

//zamykam plik
fclose(out);

//zwalniam pamięć
free_matrix(A,1,N,1,N);
free_matrix(L,1,N,1,N);
free_matrix(U,1,N,1,N);
free_matrix(Awyn,1,N,1,N);
free_matrix(A_iloczyn,1,N,1,N);
free_ivector(indxA, 1, N);
free_vector(b1,1,N);
free_vector(b2,1,N);
free_vector(b3,1,N);
free_vector(b4,1,N);
return 0;
}