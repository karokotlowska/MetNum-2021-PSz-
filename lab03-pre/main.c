#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/nrutil.c"
#include "math.h"
#include <stdio.h>
#define V0 0
#define x0 1
#define w 1
#define N 2000
#define h 0.02
void Jakobi(float B, float F0, float omega);


//pamietajcie zeby kompilowac z -lm
int main()
{
  float Beta, F0, omega;
  //pierwsze dane
   Beta = 0.0;
   F0 = 0.0;
   omega = 0.8;
   Jakobi(Beta, F0, omega);

  //drugie dane
  // B = 0.4;
  // F0 = 0.0;
  // omega = 0.8;
  // Jakobi(B, F0, omega);

//3dane
  //Beta = 0.4;
  //F0 = 0.1;
  //omega = 0.8;
  //Jakobi(Beta, F0, omega);
  //return 0;
}

void Jakobi(float Beta, float F0, float omega)
{
  printf("B = %f,  F0 = %f, omega = %f \n",Beta,F0 , omega);

  float a1 = 1;
  float a2 =pow(w,2)*pow(h,2) - 2 - Beta * h;
  float a3 = 1 + Beta * h;

//utworzenie trzech wektorów, opnieważ macierz układu równań Ax=b jest trójprzekątna
  float *d0 =vector(1,N+1);
  float *d1 = vector(1,N+1);
  float *d2 =vector(1,N+1);

  float *b =vector(1,N+1);

  //inicjalizacja wektorow
  for (int i = 0; i < N; i++)
  {
    b[i] = F0 * sin(omega * h * i) * h * h;
    d0[i] = a3;
    d1[i] = a2;
    d2[i] = a1;
   
  }
  d0[0] = 1.0;
  d0[1] = 1.0;
  d1[0] = 0.0;
  d1[1] = -1.0;
  d2[0] = 0.0;
  d2[1] = 0.0;
  //uzupelniam dwie pierwsze wartosci wektora b
  b[0] = 1.0;
  b[1] = 0.0;

  //wypisz pierwsze 5 elementow wektorow
  //std::cout << "\nd0  d1  d2  b\n" ;
  for (int i = 0; i < 10; i++){
    printf ("%2f  %2f  %2f  %2f\n",d0[i],d1[i],d2[i],b[i]);
  }

  //aby wyznaczyc i-ty element nowego przybliżenia xn(i), tworze wektory xn i xs(starego przybliżeani) i wykonuje operacje (opisana równaniem 12)
  float *xn = vector(1,N+2);
  float *xs = vector(1,N+2);
  //zeruje wektory
  for (int i=0;i<N;i++)
  {
	  xn[i]=0;
	  xs[i]=0;
  }
  xs[N]=0;
  xs[N+1]=0;
  xn[N]=0;
  	
	xs[0]=888; // może to być dwowolna liczba
	xs[1]=888; // może to być dowolna liczba
	xs = xs+2;
	int iteration=1;
  
	while (iteration < 1e5)
	{
     for (int i = 0; i <= N; i++)
       xn[i]=(1/d0[i])*(b[i]-d1[i]*xs[i-1]-d2[i]*xs[i-2]);

	  iteration++;
	  double sum_xs=0;
	  double sum_xn=0;

	for (int i=0;i<N;i++)
	{
		sum_xs+=pow(xs[i],2);
		sum_xn+=pow(xn[i],2);
	}

  //sprawdzam warunek zebiżności: bezwzzględna różnica sum kwadratów przybliżenia nowego i tego z poprzedniej iteracji
	if (fabs(sum_xn-sum_xs) < 1e-6)
	break;

	for (int i = 0; i < N; i++)
      xs[i] = xn[i];
	  // std::cout<<"iteration: "<<iteration<<"\n";
}

	FILE *file_ptr;
	file_ptr = fopen ("Przypadek_3.dat","w");// tutaj zmieniamy realizowany przypadek
	if (!file_ptr)
	{
		printf("ERR while opening file! Leaving!\n"); 
		return;
	}
double start = 0.0;
  for (int i = 0; i < N; i++)
  {
    fprintf(file_ptr, "%f   %f\n", start, xn[i]);
    start += h;
  }

xs = xs-2;

 fclose(file_ptr);
}


