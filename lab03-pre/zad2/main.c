#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "nrutil.h"


#define N 1000
#define max(X,Y) ((X)>(Y)? (X):(Y))
#define min(X,Y) ((X)<(Y)? (X):(Y))
#define abs(X) ((X)>0? (X):-(X))

//iloczyn macierzy
void multiplyAvec (float** A, float* b, float* result, int n, int m){
	for (int i = 1; i <= n; i++){
		result[i] = 0.0;
		for (int j = max(1, i - m); j <= min(n, i + m); j++)
			result[i] += A[i][j] * b[j];
	}
}

//iloczyn sklaarny
float multiplyScalar(float* a, float* b, int n){
	float sum = 0.0;
	for(int i =1; i <=n; i++)
		sum += a[i] * b[i];
	return sum;
}

void fprintMatrix(FILE * file, char * txt, float ** matrix, int size) {
	fprintf(file, "%s", txt);
	for(int i = 1; i <= size; i++) {
		for(int j = 1; j <= size; j ++) {
			fprintf(file, "%12g ", matrix[i][j]);
		}
		fprintf(file, "\n");
	}
}

int main(void) {
	float **A=matrix(1,N,1,N);
	float* b=vector(1,N);
	float* x =vector(1,N);
	int m = 5;
//ten czas nie jest jakos mega potrzebny, bardziej bajer
  time_t t1,t2;
  double t21;
  FILE * file = fopen("result.dat", "w+");

//wpisywyanie wartosci do macierzy A i wektorów b i i x
	for(int i = 1; i <=N; i++) {
		for(int j = 1; j <= N; j++) {
			if (abs(i - j) <=  m)
				A[i][j] =  1./(1 + (abs(i - j)));
			else
				A[i][j] =  0;
		}
		b[i] = i;
		x[i] = 0;
	}
  //fprintMatrix(file,"Macierz A:\n", A, N);

	float* r = vector(1,N);
	float* v = vector(1,N);
	float* y = vector(1,N);
	float* tmp = vector(1,N);

	float alfa, beta;

//zwraca zmodyfikowanego y
	multiplyAvec(A, x, y, N, m);

	for (int i = 1; i <=N; i++){
		v[i] = r[i] = b[i] - y[i];
	}
//określam czas - nie jest to jakos potrzene w zadaniu
	time(&t1);
	int k = 0;

  //pętla iteracyjna
	while(sqrt(multiplyScalar(r, r, N)) > pow(10, -6)) {

		double tmpScalar = multiplyScalar(r, r, N);
		double tmpScalar2 = multiplyScalar(x, x, N);
		multiplyAvec(A, v, tmp, N, m);
		alfa = tmpScalar / multiplyScalar(v, tmp, N);

//uzupelniam odpowiednie wektory
		for(int i = 1; i <= N; i++) {
			x[i] += alfa * v[i];
			r[i] -= alfa * tmp[i];
		}
    

		beta = multiplyScalar(r, r, N) / tmpScalar;

		for(int i = 1; i <= N; i++) {
			v[i] = v[i] * beta + r[i];
		}
if(k==0){
  //wypisuje pierwsza linijke, w ktorej nie ma być alfa i beta
  fprintf(file, "%d  %g %g %g %g\n", k, sqrt(tmpScalar),0., 0., sqrt((multiplyScalar(x, x, N))));
}
else{
  //wypisuje pozostałe linijki
		fprintf(file, "%d %g %g %g %g\n", k, sqrt(tmpScalar), alfa, beta, sqrt(tmpScalar2));
}
		k++;
	}

	fclose(file);

	time(&t2); 

	printf("Part two execution time: ");
  t21=difftime(t2,t1);
  printf("%f",t21);

	free_matrix (A,1,N,1,N);
	free_vector (b,1,N);
	free_vector (x,1,N);
	free_vector (v,1,N);
	free_vector (r,1,N);
	free_vector (y,1,N);
	free_vector (tmp,1,N);



	return 0;
}