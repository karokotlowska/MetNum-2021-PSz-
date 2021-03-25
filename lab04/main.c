#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/tred2.c"
#include "numerical_recipes.c/nrutil.c"
#include "numerical_recipes.c/tqli.c"
#include "numerical_recipes.c/pythag.c"

//okreslam parametry z zadania
const int nx = 20;
const int ny = 20;
const int n = 400; //nx*ny=20*20
const int m = 10;
const float t = -0.021;

//procedura mnożąca macierze wpisująca wynik mnożenia do macierzy C
void pomnoz(float **A,float **B,float ***C){
  float sum=0;
  for(int i=1;i<=n;i++){
    for(int j=1;j<=n;j++){
      sum=0;
      for(int k=1;k<=n;k++){
        sum+=A[i][k]*B[k][j];
      }
      (*C)[i][j]=sum;
    }
  }
}

void sortujEnergie(float **d,int **indx){
  float e1,e2,l1,l2;
  for(int l=1;l<=n;l++) (*indx)[l]=l;
  for(int l=1;l<=n-1;l++){
    for(int k=n;k>=l+1;k--){
      e1=(*d)[k-1];
      e2=(*d)[k];
      l1=(*indx)[k-1];
      l2=(*indx)[k];
      if(e2<e1){ //wymieniamy energie i indeksy wektorów miejscami
        (*d)[k]=e1;
        (*d)[k-1]=e2;
        (*indx)[k]=l1;
        (*indx)[k-1]=l2; 
      }
    }
  }
}

int main(void) {
//inicjalizuje macierze i wektory
float **H=matrix(1,n,1,n);
float **Y=matrix(1,n,1,n);
float **X=matrix(1,n,1,n);
float *d=vector(1,n);
float *e=vector(1,n);
int *index=ivector(1,n);

//inicjalizuje macierz H, kodem z treści zadania
int l = 0;
for (int i = 1; i <= nx; i++){
  for (int j = 1; j <= ny; j++){
    l = j + (i - 1) * ny;
    for (int k = 1; k <= n; k++)
      H[l][k] = 0.;
      if (i > 1)
        H[l][l - ny] = t; //dla i=1 nie ma sasiada z lewej strony
      if (i < nx)
        H[l][l + ny] = t; //dla i=nx nie ma sasiada z prawej strony
        H[l][l] = -4 * t;
      if (j > 1)
        H[l][l - 1] = t; //dla j=1 nie ma sasiada ponizej siatki
      if (j < ny)
        H[l][l + 1] = t; //dla j=ny nie ma sasiada powyzej siatki
    }
  }

//wypelniam wektory e i d zerami
for(int i=0;i<=n;i++){
  e[i]=0;
  d[i]=0;
}
//uzupelniam macierz Y
for(int i = 1; i <=n; i++){
 	for(int j = 1; j <= n; j++){
 		if(i == j){
      Y[i][j] = 1;
    }
 		else Y[i][j] = 0;
 	}
}

//przekształcam macierz H do postaci trójdiagonalnej przy użyciu procedury tred2, zwraca macierz trójdiagonalną T zapisaną w postaci wektorów d(diagonalna) i e (pierwsza poddiagonalna)
tred2(H,n,d,e);
//diagonalizuje macierz T przy pomocy tqli
tqli(d,e,n,Y);
//przekształcamy wszzytskie wektory czyli dokonujemy mnożenia dwóch macierzy
pomnoz(H,Y,&X);
//sortuje energie
sortujEnergie(&d,&index);
//wpisujemy wektory własne do pliku
FILE *fp;
fp=fopen("dane.dat","w");
for(int i=1;i<=nx;i++){
  for(int j=1;j<=ny;j++){
    l=j+(i-1)*ny;
    fprintf(fp,"%6d %6d ",i,j);
    for(int k=1;k<=m;k++)
      fprintf(fp," %12.6f ",X[l][index[k]]);
      fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
}
fclose(fp);

//do drugiego pliku wpisuję m wartości własnych(wpisuje wektor d)
FILE *fp2;
fp2=fopen("dane2.dat","w");
fprintf(fp2,"m wartości własnych:\n");
for(int i=1;i<=m;i++){
  fprintf(fp2," %.9f\n",d[i]);
}
fclose(fp2);
//zwalniam pamięc
free_matrix(H, 1, n, 1, n);
free_matrix(Y, 1, n, 1, n);
free_matrix(X, 1, n, 1, n);
free_vector(d, 1, n);
free_vector(e, 1, n);
free_ivector(index,1,n);

return 0;
}