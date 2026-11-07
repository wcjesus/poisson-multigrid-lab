#include <stdio.h>
#include <math.h>
#include "matriz.h"
#include "residuo.h"

/*
Calcula o residuo:

   r = f - A*u


    A*u =
   (4u - vizinhos)/h^2

*/

void calcula_residuo(Matriz u, Matriz f,Matriz r,int n, double h)
{
int i;
double h2 =h*h; 
zera_matriz(r, n);

#ifdef DIM2
for (i = 1; i < n; i++)
{
    for (int j = 1; j < n; j++)
    {
        double Au;

        Au =
            (4.0 * u[i][j]
             - u[i-1][j]
             - u[i+1][j]
             - u[i][j-1]
             - u[i][j+1])
            / h2;

        r[i][j] = f[i][j] - Au;
    }
}
#elif defined(DIM3)
    for (i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            for (int k = 1; k < n; k++) 
	    {
                double Au = (6.0*u[i][j][k]
                           - u[i-1][j][k] - u[i+1][j][k]
                           - u[i][j-1][k] - u[i][j+1][k]
                           - u[i][j][k-1] - u[i][j][k+1]) / h2;
                r[i][j][k] = f[i][j][k] - Au;
            }
#else
for (i=1;i<n;i++)
{
	double Au;
	Au = 
		(2.0*u[i]
		 -u[i-1]
		 -u[i+1])
		/ h2;
	r[i] = f[i] -Au;
}
#endif


}


double norma_residuo(Matriz r,int n,double h)
{
int i;
double soma = 0.0;
#ifdef DIM2
for (i = 1; i < n; i++)
{
    for (int j = 1; j < n; j++)
    {
        soma += r[i][j] * r[i][j];
    }
}

return sqrt(h * h * soma);
#elif defined(DIM3)
    for (i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            for (int k = 1; k < n; k++)
                soma += r[i][j][k] * r[i][j][k];
    return sqrt(h*h*h*soma);
#else
    for (i=1;i<n;i++)
    	soma += r[i]*r[i];
    return sqrt(h*soma);
#endif
}
