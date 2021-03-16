#include <stdio.h>
#include <math.h>
#include "nrutil.c"
#include "gaussj.c"


#define N 1000
#define m 10
#define max(X,Y) ((X)>(Y)? (X):(Y))
#define min(X,Y) ((X)<(Y)? (X):(Y))
#define abs(X) ((X)>0? (X):-(X))


void matrix_vector(float **A, float *x, float *y){
    int jmin, jmax, i,j;
    for(i=0;i<N;i++){
        jmin=max(0,i-m);
        jmax=min(i+m,N-1);
        y[i]=0;
        for(j=jmin;j<=jmax;j++)y[i]+=A[i][j]*x[j];
    }
}


float vector_vector(float *x, float *y){
    float s = 0;
    for(int i=0; i<N; i++){
        s+=x[i]*y[i];
    }
    return s;
}
void add_vector(float *a, float *b, float *c, float k){
     for(int i=0; i<N; i++){
            c[i] = a[i] + k*b[i];
        }
}

int main(){
    FILE *f1= fopen("a.dat", "w");
    FILE *f2 = fopen("b.dat", "w");
    float **A = malloc(N*sizeof(float*));
    for (int i=0; i<N; i++)
        A[i] = malloc(N*sizeof(float));
    float *b = malloc(N*sizeof(float));
    float *x = malloc(N*sizeof(float));
    float *r = malloc(N*sizeof(float));
    float *t = malloc(N*sizeof(float));


    for(int i=0; i<N; i++){
        b[i] = i;
        x[i] = 0;
        r[i] = 0;
        t[i] = 0;
        for(int j=0;j<N; j++){
            if(abs(i-j)<=m)
                A[i][j] = 1.0/(1.0+abs(i-j));
            else
                A[i][j] = 0;
        }
    }
    float norma_r=0;
    float norma_x=0;
    int k=0;
    float alpha = 0;
    do{
        k++;
        matrix_vector(A, x, t);//t = A*x
        for(int i=0; i<N; i++){
            r[i] = b[i] - t[i];
        }
        matrix_vector(A, r, t);
        alpha = vector_vector(r,r)/vector_vector(r,t);
        add_vector(x, r, x, alpha);
        norma_x = sqrt(vector_vector(x,x));
        norma_r = sqrt(vector_vector(r,r));
        fprintf(f1,"%d %f %f %f\n", k, norma_r, alpha, norma_x );

    }while(norma_r > 1e-6 && k<500);
    
    for(int i=0; i<N; i++){
        x[i] = 1;
        r[i] = 0;
        t[i]=0;
    }
    k=0;
    norma_x = 0;
    norma_r = 0;

    
    do{
        k++;
        matrix_vector(A, x, t);//t = A*x
        for(int i=0; i<N; i++){
            r[i] = b[i] - t[i];
        }
        matrix_vector(A, r, t);
        alpha = vector_vector(r,r)/vector_vector(r,t);
        add_vector(x, r, x, alpha);
        norma_x = sqrt(vector_vector(x,x));
        norma_r = sqrt(vector_vector(r,r));
        fprintf(f2,"%d %f %f %f\n", k, norma_r, alpha, norma_x );

    }while(norma_r > 1e-6 && k<500);
    
    fclose(f1);
    fclose(f2);
    free(b);
    for (int i=0; i<N; i++)
    free(A[i]);
    free(A);
    free(r);
    free(x);

  


    return 0;
}