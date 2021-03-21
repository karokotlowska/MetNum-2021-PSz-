#include "stdio.h"
#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/tred2.c"
#include "numerical_recipes.c/nrutil.c"
#include "numerical_recipes.c/tqli.c"
#include "numerical_recipes.c/pythag.c"
//przyjmuję parametry ktore są opisane w zadaniu
const int nx = 20;
const int ny = 20;
const int n = 400;
const int m = 10;
const float t = -0.021;




void pomnoz(float **A, float **B, float ***C)
{
    float suma;
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            suma = 0;
            for (int k = 1; k <= n; k++)
            {
                suma += A[i][k] * B[k][j];
            }
            (*C)[i][j] = suma;
        }
    }
}

void sort(float **d, int **index)
{
    float e1, e2, l1, l2;
    for (int l = 1; l <= n; l++)
        (*index)[l] = l; // inicjalizacja
    for (int l = 1; l <= n - 1; l++)
    {
        for (int k = n; k >= l + 1; k--)
        {
            e1 = (*d)[k - 1];
            e2 = (*d)[k];
            l1 = (*index)[k - 1];
            l2 = (*index)[k];
            if (e2 < e1)
            { //wymieniamy energie i indeksy wektorĂłw miejscami
                (*d)[k] = e1;
                (*d)[k - 1] = e2;
                (*index)[k] = l1;
                (*index)[k - 1] = l2;
            }
        }
    }
}


void file_print(float **X, int *index, FILE *fp)
{
    int l;
    for (int i = 1; i <= nx; i++)
    {
        for (int j = 1; j <= ny; j++)
        {
            l = j + (i - 1) * ny;
            fprintf(fp, "%6d %6d ", i, j);
            for (int k = 1; k <= m; k++)
                fprintf(fp, " %12.6f ", X[l][index[k]]);
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }
}

int main(void)
{
//inicjalizuje macierze i wektory
    float **H=matrix(1, n, 1, n);
    float **Y=matrix(1, n, 1, n);
    float **X=matrix(1, n, 1, n);
    float *d=vector(1, n);
    float *e=vector(1, n);
    int *index=ivector(1,n);
//otwieram plik
    FILE *file_ptr;
    file_ptr = fopen("dane.dat", "w");

    if (file_ptr)
        return 1;
    else
        return 0;

//inicjalizuje macierz H kodem ktory został podany w treści zadania, można sobie uprościć i wrzucić to do funkcji void (float ***matrix) i wtedy wywolac funkcję z maina np: inicjalizuj(&H)
    int l = 0;
    for (int i = 1; i <= nx; i++)
    {
        for (int j = 1; j <= ny; j++)
        {
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
//wypelniam wszystko inne zerami (chyba jest to ważne)
 for(int i = 1; i <= n; i++){
   //dwa wektory
    	e[i] = 0;
    	d[i] = 0;
    }
//oraz macierz Y
for(int i = 1; i <=n; i++){
    	for(int j = 1; j <= n; j++){
    		if(i == j) Y[i][j] = 1;
    		else Y[i][j] = 0;
    	}
    }

//przechodzę do punktu 4 z treści czyli "Przekształcamy macierz H do postaci trójdiagonalnej" - robimy to przy pomocy funkcji tred2, która zwraca macierz trójdiagolnalną T zapisaną w postaci wektorów d(diagonalna) i e(pierwsza poddiagonalną T), macierz H jest nadpisana przez macierz podobieństwa P (nie wiem za bardzo co to w tym moemencie, moze zaraz się wyjasni)
    tred2(H, n, d, e);
//teraz chodzi o to żeby zdiagonalizować macierz T, uzywamy procedury tqli, w Y zostaną wpisane wektory własne macierzy T, a wartości własne zostaną wpisane do d
    tqli(d, e, n, Y);
//wymnażamy macierze - punkt 6 w treści
    pomnoz(H, Y, &X);
//sortujemy energie - kod z tej fukncji jest podany w treści zadania
    sort(&d, &index);
//wpisujemy do pliku
    file_print(X, index, file_ptr);


    free_matrix(H, 1, n, 1, n);
    free_matrix(Y, 1, n, 1, n);
    free_matrix(X, 1, n, 1, n);
    free_vector(d, 1, n);
    free_vector(e, 1, n);
    free_ivector(index, 1, n);
    fclose(file_ptr);
}
