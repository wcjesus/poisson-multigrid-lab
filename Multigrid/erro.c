#include <math.h>
#include "problema.h"
#include "matriz.h"
#include "erro.h"

double norma_erro(Matriz u,int n,double h)
{
int i;
double soma = 0.0;
double erro;
#ifdef DIM2
for (i = 1; i < n; i++)
{
    for (int j = 1; j < n; j++)
    {
        erro =u[i][j]- uex(i * h, j * h);
        soma += erro * erro;
    }
}

return sqrt(h * h * soma);
#elif DIM3
    for (i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            for (int k = 1; k < n; k++) 
	    {
                erro = u[i][j][k] - uex(i*h, j*h, k*h);
                soma += erro*erro;
            }
    return sqrt(h*h*h*soma);
#else 
    for (i=1;i<n;i++)
    {
	    erro = u[i]-uex(i*h);
	    soma += erro*erro;
    }
    return sqrt(h*soma);
#endif
}

double norma_infinito_erro(Matriz u, int n, double h)
{
    double maxe = 0.0;
#ifdef DIM2
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++) 
	{
            double e = fabs(u[i][j] - uex(i*h, j*h));
            if (e > maxe) maxe = e;
        }
#elif DIM3
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            for (int k = 1; k < n; k++) 
	    {
                double e = fabs(u[i][j][k] - uex(i*h, j*h, k*h));
                if (e > maxe) maxe = e;
            }
#else
    for (int i=1; i<n; i++)
    {
	double e = fabs(u[i] - uex(i*h));
	if (e>maxe) maxe = e;
    }	
#endif
    return maxe;
}
