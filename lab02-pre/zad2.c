#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/nrutil.c"
#include "numerical_recipes.c/ludcmp.c"
#include "numerical_recipes.c/lubksb.c"

#define N 4

void printMatrix(float** Matrix) //wyswietla macierze
{
	int i,j;
	printf("\n");
	for(i = 1; i <=N; i++)
	{
		for(j = 1; j <= N; j++)
		{
			printf("%f ", Matrix[i][j]);
		}
		printf("\n");
	}
	
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

float ** cpyMatrix(float **A, float **B)		//kopiuje macierze
{
	int i,j;
	for(i=1;i<=N;i++) 
	{
		for(j=1;j<=N;j++)
		{
			B [i][j] = A[i][j];				
		}
	}
	return B;
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

int main(void)
{
	float ** A = matrix(1,N,1,N);	//macierz A
	
	
	FILE * data = fopen("Lab2_MN.txt", "w+");	

	if (data == NULL )
    {
        perror("Nie udalo sie otworzyc pliku");
        return 1;
    }
	
	int i,j;
	
/////UZUPELNIENIE MACIERZY A//////////////

	for(i=1;i<=N;i++) 
	{
		for(j=1;j<=N;j++)
		{
			A[i][j]=1.0/(i+j);				
		}
	}
	
	printMatrix(A);	
	
//tworze wektory premutacji typu int
int *indxA=ivector(1, N);
//d_inicjalizuj_wektor(&indxA);

//tworze zmienne typu float
float mA;

//wykonuje procedurę ludcmp
ludcmp(A, N, indxA, &mA);

//drukuje macierze po rozkładzie
printf("MACIERZ LU:");
printMatrix(A);


//tworze macierz LU do ktorej zapisze rozlonona macierz mA
float **LU = matrix(1,N,1,N);

///////ROZKALD MACIEZRY A NA LU////////////////////
  /*Kopiowanie macierzy A do macierzy LU*/
	cpyMatrix(A,LU);
	/*Rozkład macierzy A na L i U*/
	ludcmp(LU, N, indxA, &mA);
	
	

///wpisanie macierzy U i Lab2_MN//macierz L dla A trójkątna górna
float **LA=matrix(1, N, 1, N);


//macierz U dla A trójkątna dolna
float ** UA  = matrix(1, N, 1, N);


//alokuje komorki maecierzy LA 
for (int i=1;i<=N;i++){
  for(int j = 1; j <= N; j ++) {
    if(i==j){
			LA[i][j] = 1;
		}
    else if(i<j){
      LA[i][j] = 0;
    }
    else{
      LA[i][j] = A[i][j];
    }
  }
}

//alokuje komorki maecierzy UA
for (int i=1;i<=N;i++){
  for(int j = 1; j <= N; j ++) {
    if(i<=j){
      UA[i][j] = A[i][j];
    }
    else{
      UA[i][j] = 0;
    }
  }
}

//drukuje macierze po rozkładzie
printf("Macierz LA i UA:\n");
printMatrix( LA);
printMatrix( UA);



//////ELEMENTY DIAGOANLI MACIEZRY U//////////////
    printf("\n");
    for(i=1; i<=N; i++){
		printf("%f \n", UA[i][i]);	
	}
	
	fprintf(data,"Diagonala macierzy U:\n");	
	for(i=1; i<=N; i++)
	{
		fprintf(data,"\t%f\n", UA[i][i]);	
	}



//////LICZYMY WYZNACZNIK MACIEZRY A///////////////////////
float detA=1;
  for(i=1; i<=N; i++)
	{
		detA*=UA[i][i];	
	}
	//detA*=mA;
	printf("\n\n%e \n\n", detA);
	
	fprintf(data,"\nWyznacznik macierzy A:\t\t");
	fprintf(data,"%e\n", detA);
	
//////MACIERZ ODWOTNA DO MACIERZY A/////////////////////////

//tworze wektory wyrazów wolnych 
float *x1=vector(1, N);
x1[1]=1;x1[2]=0;x1[3]=0;x1[4]=0;
float *x2=vector(1, N);
x2[1]=0;x2[2]=1;x2[3]=0;x2[4]=0;
float *x3=vector(1, N);
x3[1]=0;x3[2]=0;x3[3]=1;x3[4]=0;
float *x4=vector(1, N);
x4[1]=0;x4[2]=0;x4[3]=0;x4[4]=1;

lubksb(A,N,indxA,x1);
lubksb(A,N,indxA,x2);
lubksb(A,N,indxA,x3);
lubksb(A,N,indxA,x4);

//towrze macierz A^-1 wynikowa
float **Awyn=matrix(1, N, 1, N);

//wpisuje wyniki
for(int i=1;i<=N;i++){
  Awyn[1][i]=x1[i];
  Awyn[2][i]=x2[i];
  Awyn[3][i]=x3[i];
  Awyn[4][i]=x4[i];
}
printf("Macierz A");
printMatrix(A);
printf("Macierz A-1");
printMatrix(Awyn);

	


////////////ILOCZYN////////////////
	for(i=1;i<=N;i++) 
	{
		for(j=1;j<=N;j++)
		{
			A[i][j]=1.0/(i+j);				
		}
	}
//ilocznyn A*A-1
float ** A_iloczyn = matrix(1, N, 1, N);
pomnoz(Awyn,A,&A_iloczyn);
    

//drukuje macierze po iloczynie
printf("\n\nMacierz iloczynu AA^-1':\n");
printMatrix(A_iloczyn);



/////WSKAZNIK UWARUNKOWANIA MACIERZY A//////////////////////////
printf("\nmaximum macierzy A %.2f\n",max_value(A));
printf("maximum macierzy A^-1 %.2f\n",max_value(Awyn));


//wskaznik uwarunkowania macierzy cond(A)=norm(A)*norm(A-1)
printf("\n\nwskaznik uwarunkowania macierzy %.2f\n",max_value(A)*max_value(Awyn));

/////ZAMYKANIE PLIKU I ZAWALNIANIE PAMIECI///////////////////	
	fclose(data); 

	free_matrix(A,1,N,1,N);
	free_matrix(LU,1,N,1,1);
	free_matrix(Awyn,1,N,1,1);
	
//	system("pause"); 
	
	return 0;
	
	

}
