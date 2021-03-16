#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

#include "nrutil.c"
#include "gaussj.c"

/////tak naprawder to zadanie 3

#define n 1000
#define m 5
#define max(X,Y) ((X)>(Y)? (X):(Y))
#define min(X,Y) ((X)<(Y)? (X):(Y))
#define abs(X) ((X)>0? (X):-(X))

void set_A_b_x(float** _a, float* _b, float* _x);
void print_matrix(float** _a);
void multiply_A_x(float**A, float*x, float*y);
float il_skal(float* _r1, float* _r2);
void subtract(float *a, float *b, float *d);
void print_vector(float* v);
void addvector(float *a,float*b, float *d, float c)
{
    for(int i=0;i<n;i++)
    {
        d[i] = a[i] + c*b[i];
    }
}

float matrixvector(float** a, float* x, float *y){
	 int jmin, jmax;
    for(int i=0;i<n;i++){
        jmin=max(0,i-m);
        jmax=min(i+m,n-1);
        y[i]=0;
        for(int j=jmin;j<=jmax;j++)
            y[i]+=a[i][j]*x[j];
    }
}


float scalar(float *a,float *b,int flagt)
{
    float s=0.;
    if(flagt)
    {
        for(int i=0;i<n;i++)
            s+=a[i]*b[i];
    }
    else
    {
        for(int i=0;i<n;i++)
            s+=b[i]*a[n-1-i];
    }

    return s;

}


int main(){

    float **A = matrix(1, n,1,n);
    float* b=vector(1,n);
	  float* x =vector(1,n);
    float* r = vector(1,n);
    float *tmp2 =vector(1,n);
    float *tmp =vector(1,n);
    float alpha =0.;

    clock_t start, end;
    double time;

    FILE *f1 = fopen("zad_a.dat", "w");
    FILE *f2 = fopen("zad_b.dat", "w");


//######################################################################################################
// X0 = 0

    set_A_b_x(A, b, x);

    int k=1;
    double n1=0.;
    double n2=0.;
  for(int i = 1; i <=n; i++) {
		for(int j = 1; j <= n; j++) {
			if (abs(i - j) <=  m)
				A[i][j] =  1./(1 + (abs(i - j)));
			else
				A[i][j] =  0;
		}
		b[i] = i;
		x[i] = 0;
	}

    start=clock();

    do{

        matrixvector(A, x, tmp); //tp=A*x
        addvector(b, tmp, r, -1); //r = b - tmp = b - A*x
        matrixvector(A, r, tmp2); // tmp = A *r = A *( b-A*x)
        alpha = scalar(r,r,1) / scalar(r, tmp2, 1);
        addvector(x, r, x, alpha); // x = xx + alpha*r

        n1 = sqrt(scalar(r,r,1));
        n2 = sqrt(scalar(x,x,1));
        fprintf(f1, "%d %lf %lf %lf\n", k, n1, alpha, n2);

        k++;


    }while(n1>1e-3);

    end=clock();
    //time= static_cast<double>(end-start)/ static_cast<double>(CLOCKS_PER_SEC);
    printf("time 1: %lf\n", time);

set_A_b_x(A, b, x);
  //////////////////////////////////
  for(int i=0;i<n;i++)
    {
        b[i]=i+1;
        x[i]=1.;
        r[i]=0.;
        tmp[i]=0.;
    }

    k=1;
    float norma1=0.;
    float norma2=0.;

    start=clock();

    do{

        matrixvector(A, x, tmp); //tmp=A*x
        addvector(b, tmp, r, -1); //r = b - tmp = b - A*x
        matrixvector(A, r, tmp2); // tmp = A *r = A *( b-A*x)
        alpha = scalar(r,r,1) / scalar(r, tmp2, 1);
        addvector(x, r, x, alpha); // x = xx + alpha*r

        norma1 = sqrt(scalar(r,r,1));
        norma2 = sqrt(scalar(x,x,1));
        fprintf(f2, "%d %lf %lf %lf\n", k, norma1, alpha, norma2);
        k++;


    }while(norma1>1e-3);


}









void set_A_b_x(float** _a, float* _b, float* _x)
{
    for(int i=0; i<n; i++)
    {
        _a[i] = vector(1,n);
        _b[i] = i+1;
        _x[i] = 0;

        for(int j=0; j<n; j++)
        {
            if(abs(i-j) <= m)
                _a[i][j] = 1./(1 + abs(i-j));
            else
                _a[i][j] = 0;
        }
    }
}



void print_vector(float* v) 
{
    for (int i=0; i<n; i++) 
    {
        printf("%5.f ", v[i]);
    }
    printf("\n\n");
}

void multiply_A_x(float**A,float*x,float*y)
{   
    int jmin, jmax;
    for(int i=0; i<n; i++){
        jmin = max(0,i-m);
        jmax = min(i+m,n-1);
        y[i] = 0;
        for(int j=jmin; j<=jmax; j++)
            y[i] += A[i][j]*x[j];
    }
}

float il_skal(float* _r1, float* _r2)
{ 
    float s=0;

    for(int i=0; i<n; i++)
    {
        s += _r1[i]*_r2[i];
    }

    return (s);
}

void subtract(float *a, float *b, float *d)
{
    for(int i=0; i<n; i++)
    {
        d[i] = a[i] - b[i];
    }
}
