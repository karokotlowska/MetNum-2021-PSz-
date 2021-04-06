#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define n 7
#define IT_Max 12

double iloczynSkalarny(double w1[n], double w2[n])
{
	double wynik = 0.0;
	for (int i = 0; i < n; i++)
		wynik += w1[i] * w2[i];

	return wynik;
}

void iloczynMacierz_Wektor(double a[n][n], double x[n], double tmp[n])
{
	for (int i = 0; i < n; i++)
	{
		tmp[i] = 0;
		for (int j =0; j < n; j++)
		{
			tmp[i] += (a[i][j] * x[j]);
		}
	}
}
void iloczynMacierz_Macierz(double A[n][n], double B[n][n], double C[n][n]){
    double s;
    for(int i = 0; i < n; i++ )
    for(int j = 0; j < n; j++ )
    {
      s = 0;
      for(int k = 0; k < n; k++ ) s += A[i][k] * B[k][j];
      C [i][j] = s;
    }
}

int main()
{
    
    int i;
    double A[n][n];
	double W[n][n];
    double X[n][n];
    double X_T[n][n];
    double D[n][n];
    double x_old[n];
    double x_new[n];
    double lambda[n];
    FILE *file1;
    file1 = fopen("wyniki1.dat", "w");
    FILE *file2;
    file2 = fopen("mtxD.dat", "w");

for(int i=0; i<n; i++){
  for(int j=0; j<n; j++){
    A[i][j]=(1.0+abs(i+j))/(1.0+abs(i-j));  
    W[i][j]=A[i][j];
    }
}
for(int k=0; k<n; k++){ 
        //fprintf(file1, "%d ", k+1);
  double c=0.;
  for(i=0;i<n;i++){
    x_old[i]=1.;
    c=c+x_old[i];
  }
        //normalizacja wektora
  for(i=0;i<n;i++)x_old[i]=x_old[i]/sqrt(c);
    for(int m=1; m<=IT_Max; m++)
        {
        iloczynMacierz_Wektor(W,x_old, x_new);
        lambda[k]=iloczynSkalarny(x_new,x_old);
            
        for(i=0;i<n;i++) 
        x_old[i]=x_new[i]/sqrt(iloczynSkalarny(x_new,x_new));
            
        fprintf(file1, "%d  % d  %7f \n", k,m, lambda[k]);
            //printf( "%d  % d  %7f \n", k,m, lambda[k]);
        }
        fprintf(file1, "\n");
        
        
        
        for(int i=0; i<n; i++){
         for(int j=0; j<n; j++){
            W[i][j]=W[i][j]-lambda[k]*x_old[i]*x_old[j];
         }
        }
        
        for(int i=0; i<n; i++)
        {
            X[i][k]=x_old[i];
        }   
    }

    double C[n][n]; 
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            X_T[j][i]=X[i][j];
            C[i][j]=0;
        }
    }
    
    iloczynMacierz_Macierz(X_T,A, C);
    iloczynMacierz_Macierz(C,X, D);

      for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            fprintf(file2, "%12.6g ", D[i][j]);
        }
        fprintf(file2, "\n");
    }

     fclose(file1);
     fclose(file2);
}