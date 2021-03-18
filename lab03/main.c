#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//ustalam parametry początkowe
#define v0 0.1
#define x0 1.0
#define w 1.0
#define n 2000
#define h 0.02

static int version=0;

void Jakobi(float Beta,float F0,float W){

//ustalam stałe a
float a1=1.;
float a2=w*w*h*h-2.-Beta*h;
float a3=1.+Beta*h;

//tworzę wektory służace do rozwiązania układu
//najpierw 3 wektory d, bo macierz A jest trójprzekątna
float *d0=malloc((n+1)*sizeof(float));
float *d1=malloc((n+1)*sizeof(float));
float *d2=malloc((n+1)*sizeof(float));
//wektory nowego i starego przybliżenia, zwiąkszam ich rozmiar +2 ponieważ są indeksowane od -2
float *xn=malloc((n+2)*sizeof(float));
float *xs=malloc((n+2)*sizeof(float));
//tworze wektor wynikowy b
float *b=malloc((n+1)*sizeof(float));

//wypełniam odpowiednio wektory macierzy trójprzekatenj
for(int i=0;i<n;i++){
  d0[i]=a3;
  d1[i]=a2;
  d2[i]=a1;
  b[i]=F0*sin(W*h*i)*h*h;
}
d0[0]=1.;
d0[1]=1.;
d1[0]=0.;
d1[1]=-1.;
d2[0]=0.;
d2[1]=0.;

b[0]=1.;
b[1]=0.;

//wypełniam czymkolwiek wektory xn i xs
for(int i=0;i<=n;i++){
  xn[i]=0.;
  xs[i]=0.;
}
//xs[n+1]=0.;
xs[0]=1.; //losowe liczby
xs[1]=1.;

//sprawdzam poprawnosc parametrow a1 a2 a3
printf("Parametry a1 a2 a3: %f\t%f\t%f\n",a1,a2,a3);

//otwieram plik odpowiednio do przykladu
FILE *out;
if(version==0){
  out=fopen("wyniki1.txt","w+");
  version++;
}
else if(version==1){
  out=fopen("wyniki2.txt","w+");
  version++;
}
else if(version==2){
  out=fopen("wyniki3.txt","w+");
  version++;

}

//pętla iteracyjna
int k=0;
float sumaXn=0.;
float sumaXs=0.;
float krok=0.;
while(k<1e4){
  for(int i=0;i<=n;i++){
    xn[i]=(1./d0[i])*(b[i]-d1[i]*xs[i-1]-d2[i]*xs[i-2]); 
  }
  k++;
  for(int i=0;i<n;i++){
    xs[i]=xn[i];
  }
  
}
//wpisuje wyniki do pliku
 for(int i=0;i<n;i++){
    fprintf(out,"%f\t%f\n",krok,xn[i]);
    krok+=h;
  }

//zamykam plik

fclose(out);

free(d0);
free(d1);
free(d2);
free(xn);
free(xs);
free(b);
}


int main(void) {

///pierwszy przypadek
float Beta1=0.0;
float F01=0.0;
float W1=0.8; //Omega


///pierwszy przypadek
float Beta2=0.4;
float F02=0.0;
float W2=0.8; //Omega

///pierwszy przypadek
float Beta3=0.4;
float F03=0.1;
float W3=0.8; //Omega


//funkcje wywoluje zmieniajac parametry dla odpowiednich przypadków
Jakobi(Beta1,F01,W1);
Jakobi(Beta2,F02,W2);
Jakobi(Beta3,F03,W3);




return 0;
}