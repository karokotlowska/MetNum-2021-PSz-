#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/tred2.c"
#include "numerical_recipes.c/nrutil.c"
#include "numerical_recipes.c/tqli.c"
#include "numerical_recipes.c/pythag.c"

#define n 7

void mnoz(float *tab1, float *tab2, float ***tab3, int x, int middle, int y)
{
	for (int i = 1; i <= n; i++){
		for (int j = 1; j <= n; j++){
			(*tab3)[i][j] = tab1[i]*tab2[j];
		}
	}
}

void lambdaRazyMacierz(float l,float ***m){
  for (int i = 1; i <= n; i++){
		for (int j = 1; j <= n; j++){
			(*m)[i][j]*=l;
		}
	}
}

void macierzRazyWektor(float **m,float *v,float **res){
  float sum;
  for(int i=1;i<=n;i++){
    sum=0;
    for(int j=1;j<=n;j++){
      sum+=m[i][j]*v[j];
    }
    (*res)[i]=sum;
  }
}

void hotteling(float **m1,float **m2,float ***m3){
  for (int i = 1; i <= n; i++){
		for (int j = 1; j <= n; j++){
			(*m3)[i][j] = m1[i][j]-m2[i][j];
		}
	}
}

float skalarny(float *v1,float *v2){
  float sum=0;
  for(int i=1;i<=n;i++){
    sum+=v1[i]*v2[i];
  }
  return sum;
}

void dzielWektor(float **v,float k){
  for (int i=1;i<=n;i++)
	{
		(*v)[i] /=k;
	}
}

float norma(float *v){
  float sum = 0;
	for (int i=1;i<=n;i++){
		sum += v[i] * v[i];
	}
	return sqrt(sum);
}

int main(void) {
//inicjalizuje macierze i wektory
float **A=matrix(1,n,1,n);
float **Z=matrix(1,n,1,n);
float **W=matrix(1,n,1,n);
float **temp=matrix(1,n,1,n);
float *d=vector(1,n);
float *e=vector(1,n);
float *x=vector(1,n);
float *stVec=vector(1,n);

FILE *f1;
//tu wpisuje wartosci wlasne z metody tred i tqli
f1=fopen("wartosciwlasne.dat","w");
FILE *f2;
//tu wpisuje macierz A
f2=fopen("macierzA.dat","w");
FILE *f3;
//tu wpisuje wyniki z metody iteracyjnej
f3=fopen("lambda.dat","w");

float lambda=1.;

//uzupelniam macierz A i wpisuje do pliku macierzA.dat
for(int i=1;i<=n;i++){
  for(int j=1;j<=n;j++){
    A[i][j]=sqrt(i+j);
    fprintf(f2,"%g\t",A[i][j]);
  }
}

//sprowadzam macierz do postaci trójdiagonalnej
tred2(A,n,d,e);
tqli(d,e,n,A);

//wypisuje wartosci wlasne do pliku na podstawie wetkroa d
for(int i=1;i<=n;i++){
  fprintf(f1,"%g\t",d[i]);
}


//uzupelniam macierz W macierzą A
for(int i=1;i<=n;i++){
  for(int j=1;j<=n;j++){
    W[i][j]=sqrt(i+j);
  }
}

//uzupelniam macierz Z
for(int i=1;i<=n;i++){
  for(int j=1;j<=n;j++){
    if(j==i){
      Z[i][j]=1;
    }
  }
}

//zaczynam pętle iteracyjną

for(int k=1;k<=n;k++){
  //uzueplniam wektor startowy jedynkami
  for(int i=1;i<=n;i++){
    stVec[i]=1;
  }
  fprintf(f3, "wartosc wlasna: %5g,\n", d[8 - k]);
  //wpisuje do pliku wektory wlasne
  fprintf(f3,"wektor wlasny: ");
	for (int i=1;i<=n;i++){
		  fprintf(f3,"%5g ",A[i][k]);
  }
  fprintf(f3,"\n\n");

  for(int i=0;i<8;i++){
    macierzRazyWektor(W,stVec,&x);
    float sk1=(skalarny(x,stVec));
    float sk2=(skalarny(stVec,stVec));
    lambda=sk1/sk2;
    float norm = norma(x);
    dzielWektor(&x,norm);
    for(int i=1;i<=n;i++){
      stVec[i]=x[i];
    }
    fprintf(f3, "l%d =  %8g \n", i + 1, lambda);
  }
  fprintf(f3,"\n\n");
  mnoz(stVec, stVec, &temp, n, 1, n);
  lambdaRazyMacierz(lambda,&temp);
  hotteling(W,temp,&Z);
  for (int i = 1; i <= n; i++){
		for (int j = 1; j <= n; j++){
			W[i][j] = Z[i][j];
		}
	}
}



//zwalniam pamięc
free_matrix(A, 1, n, 1, n);
free_matrix(Z, 1, n, 1, n);
free_matrix(W, 1, n, 1, n);
free_vector(d, 1, n);
free_vector(e, 1, n);
free_matrix(temp, 1, n,1,n);
free_vector(x,1,n);
free_vector(stVec,1,n);

return 0;
}