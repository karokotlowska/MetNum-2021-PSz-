#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/nrutil.c"
#include "numerical_recipes.c/gaussj.c"
/* Dyrektywy zakladajace, ze te trzy pliki sa skopiowane do aktualnego katalogu. */
//#include "nrutil.c" // To mozna usunac, jesli plik jest dodany w poleceniu kompilacji.
//#include "gaussj.c" // To tez mozna usunac, jesli plik jest dodany w poleceniu kompilacji.

/* Dyrektywy dla Taurusa (nie wymagaja kopiowania plikow, ale Taurus musi dzialac...) */
// #include "/opt/NR/numerical_recipes.c/nrutil.h"
// #include "/opt/NR/numerical_recipes.c/nrutil.c"
// #include "/opt/NR/numerical_recipes.c/gaussj.c"

#define N 200 // rozmiar macierzy M: NxN
void generate_plot(){
  FILE* gnuplot_pipe = popen ("gnuplot -persistent", "w");
  fprintf(gnuplot_pipe, "set term png \n");
  fprintf(gnuplot_pipe, "set out \"z1.png\" \n");
  fprintf(gnuplot_pipe, "set xl \"t\" \n");
  fprintf(gnuplot_pipe, "set yl \"x(t)\" \n");
  fprintf(gnuplot_pipe, "p cos(x), \"t.txt\" u 1:2 w p lt 3 pt 3 t \"h=0.1'\" \n");
  pclose(gnuplot_pipe);
}

int main(void)
{
  
	float **M, **b;
	//	Alokacja macierzy
	M = matrix(1, N, 1, N);
	b = matrix(1, N, 1, 1);

    float vo=0;
    float h=0.1;
    float omega=1.0;
	// 	Wypelnienie macierzy M i wektora b
  for (int i = 1; i <= N; ++i){
		b[i][1] = 0.0;
		for (int j = 1; j <= N; ++j)
			M[i][j] = 0.0;
	}
b[1][1] = 1;
b[2][1] = vo*h;
	for (int i = 1; i <= N; ++i)
	{
		for (int j = 1; j <= N; ++j){
      if(i==j){
        M[i][j] = 1.0;
        }
      if(i-1==j){
        M[i][j]=omega*omega*h*h-2;
        }
      if(i-2==j){
         M[i][j] = 1.0;
        }
		}

	}
	M[2][1] = -1;
 



	//	Rozwiazanie ukladu rownan Mx=b - wywolanie procedury:
	gaussj(M, N, b, 1);

	//	Wypisanie rozwiazania, ktore procedura gaussj(M, N, b, 1); zapisala w wektorze b.
	for (int i = 1; i <= N; ++i)
		printf("%g\n", b[i][1]);

	FILE *ptr;
  ptr = fopen("t.txt", "w");

    for (int i = 1; i <= N; ++i) {
        printf("%f\n", b[i][1]);
        fprintf(ptr,"%f %f\n",i*h, b[i][1]);
    }

    fclose(ptr);
generate_plot();
	//	Zwolnienie pamieci
	free_matrix(M, 1, N, 1, N);
	free_matrix(b, 1, N, 1, 1);

  

	return 0;
}
