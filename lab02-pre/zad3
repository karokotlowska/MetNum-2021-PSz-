#include <stdio.h>
#include "math.h"
#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/nrutil.c"
#define n 500
int main()
{
  
    double Xa=0.5, Xb=2;
    double h=2*Xb/(n-1);
	  float *a=vector(1,n+1); 
    float *c=vector(1,n+1); 
    float *d=vector(1,n+1);
    float *u=vector(1,n+1);
    float *l=vector(1,n+1);  
    float *y=vector(1,n+1); 
    float *x=vector(1,n+1); 
    float *b=vector(1,n+1);
    float *V_n=vector(1,n+1);
    float *V_t=vector(1,n+1);

	//V_n rozw numeryczne, V_t rozw teoretyczne   
    
     for(int i=1; i<=n; i++){
        x[i]=-Xb+h*(i-1);
         d[i]=-2/(pow(h,2));
         a[i]=1/(pow(h,2));
         c[i]=1/(pow(h,2));
        }
    for(int i=1; i<=n; i++){
      if(x[i]>=-Xb && x[i]<-Xa ||x[i]==0|| x[i]>Xa && x[i]<=Xb)
        b[i]=0;
      if(x[i]>=-Xa && x[i]<0) {
        b[i]=-1;
      }
     if(x[i]>0 && x[i]<=Xa){
       b[i]=1;
     } 
    }  

/*KOD KINGI
  for(int i=1; i<=n; i++){
        x[i]=-Xb+h*(i-1);
         d[i]=-2/(h*h);
         a[i]=1/(h*h);
         c[i]=1/(h*h);
      if(x[i]>=-Xb && x[i]<-Xa ||x[i]==0|| x[i]>Xa && x[i]<=Xb)
        b[i]=0;
      if(x[i]>=-Xa && x[i]<0) {
        b[i]=-1;
      }
     if(x[i]>0 && x[i]<=Xa){
       b[i]=1;
     } 
    }  
*/


  //warunki brzegowe
      d[1]=1;
      c[1]=0;
      b[1]=0;
      d[n]=1;
      a[n]=0;
      b[n]=0;
        
//liczenie macierzy Li U dla macierzy trójdiagonalnej
    u[1]=d[1];
    for(int i=2; i<=n; i++){
        l[i]=a[i]/u[i-1];
        u[i]=d[i]-l[i]*c[i-1];                         
    }

//rozwiązujemy układ Ly=b     
     y[1]=b[1];
     for(int i=2; i<=n; i++)
        y[i]=b[i]-l[i]*y[i-1];

//rozwiązujemy ukald Uv=y        
     V_n[n]=y[n]/u[n]; //wzór 11
    
     for(int i=n-1; i>=1; i--){
        V_n[i]=(y[i]-c[i]*V_n[i+1])/u[i];   //wzór 12
      }
 //rozwiązanie dokladne         
      for(int i=1; i<=n; i++){
        if(x[i]>=-Xb && x[i]<=Xa) {
          V_t[i]=x[i]/16 + 1./8;
        }
        if(x[i]>=-Xa && x[i]<=0){
           V_t[i]=-(x[i]*x[i])/2 - (7*x[i])/16;
        }
        if(x[i]>=0 && x[i]<=Xa) {
          V_t[i]=(x[i]*x[i])/2- (7*x[i])/16;
        }
        if(x[i]>=Xa   && x[i]<=Xb) {
          V_t[i]=x[i]/16 - 1./8;
        }
      }  
FILE *fp;
fp=fopen("wynik.txt","w");
     /* FILE *fp2;
      fp2=fopen("wynik2.txt","w");*/
      for(int i=1; i<=n; i++){
        fprintf(fp," %12.4g      %12.4g      %12.4g  \n",x[i], V_n[i],V_t[i]);
      }
      /*for(int i=1; i<=n; i++){
        fprintf(fp2," %12.4g      %12.4g  \n",x[i],  V_t[i]);
      }*/


      fclose(fp); 
      //fclose(fp2); 

      free_vector(a,1,n);
      free_vector(b,1,n);
      free_vector(d,1,n);
      free_vector(u,1,n);
      free_vector(l,1,n);
      free_vector(y,1,n);
      free_vector(x,1,n);
      free_vector(c,1,n);
      free_vector(V_n,1,n);
      free_vector(V_t,1,n);
return 0;
}
