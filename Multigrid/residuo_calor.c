#include <math.h>
#include "matriz.h"
#include "residuo_calor.h"

void calcula_residuo_calor(Matriz u, Matriz f, Matriz r, int n, double h, double dt)
{
    double alpha = dt/(h*h);
#ifdef DIM2
    double diag = 1.0 + 4.0*alpha;
#elif defined(DIM3)
    double diag = 1.0 + 6.0*alpha;
#else 
    double diag = 1.0 + 2.0*alpha;
#endif
    zera_matriz(r, n);
#ifdef DIM2
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++) {
            double Au = diag*u[i][j]
                      - alpha*(u[i-1][j] + u[i+1][j]
                             + u[i][j-1] + u[i][j+1]);
            r[i][j] = f[i][j] - Au;
        }
#elif defined(DIM3)
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            for (int k = 1; k < n; k++) {
                double Au = diag*u[i][j][k]
                          - alpha*(u[i-1][j][k] + u[i+1][j][k]
                                 + u[i][j-1][k] + u[i][j+1][k]
                                 + u[i][j][k-1] + u[i][j][k+1]);
                r[i][j][k] = f[i][j][k] - Au;
            }
#else
    for (int i=1; i<n; i++){
	    double Au = diag*u[i]
		      - alpha*(u[i-1]+u[i+1]);
	    r[i] = f[i]-Au;
    }

#endif
}

double norma_residuo_calor(Matriz r, int n, double h)
{
    double soma = 0.0;
#ifdef DIM2
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            soma += r[i][j]*r[i][j];
    return sqrt(h*h*soma);
#elif defined(DIM3)
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            for (int k = 1; k < n; k++)
                soma += r[i][j][k]*r[i][j][k];
    return sqrt(h*h*h*soma);
#else
    for (int i = 1; i<n; i++)
	    soma += r[i]*r[i];
    return sqrt(h*soma);
#endif
}
