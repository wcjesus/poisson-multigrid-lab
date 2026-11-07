#include "matriz.h"
#include "suavizacao_calor.h"

void suaviza_calor(Matriz u, Matriz u_novo, Matriz f, int n, double h, double dt)
{
    double alpha = dt/(h*h);
#ifdef DIM2
    double diag = 1.0 + 4.0*alpha;
#elif defined(DIM3)
    double diag = 1.0 + 6.0*alpha;
#else
    double diag = 1.0 + 2.0*alpha;
#endif
    copia_matriz(u_novo, u, n);
#ifdef DIM2
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++) {
            double viz = u[i-1][j] + u[i+1][j]
		       + u[i][j-1] + u[i][j+1];
            double uj = (f[i][j] + alpha*viz)/diag;
            u_novo[i][j] = (1.0-OMEGA_CALOR)*u[i][j] + OMEGA_CALOR*uj;
        }
#elif defined(DIM3)
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            for (int k = 1; k < n; k++) {
                double viz = u[i-1][j][k] + u[i+1][j][k]
                           + u[i][j-1][k] + u[i][j+1][k]
                           + u[i][j][k-1] + u[i][j][k+1];
                double uj = (f[i][j][k] + alpha*viz)/diag;
                u_novo[i][j][k] = (1.0-OMEGA_CALOR)*u[i][j][k] + OMEGA_CALOR*uj;
            }
#else
    for (int i = 1; i < n; i++)
    {
	    double viz = u[i-1] + u[i+1];
	    double uj = (f[i] + alpha*viz)/diag;
	    u_novo[i] = (1.0-OMEGA_CALOR)*u[i] + OMEGA_CALOR*uj;
    }
#endif
}

void suavizacoes_calor(Matriz u, Matriz u_novo, Matriz f, int n, double h, double dt, int numero)
{
    Matriz a = u;
    Matriz b = u_novo;
    for (int it = 0; it < numero; it++) {
        suaviza_calor(a, b, f, n, h, dt);
        Matriz tmp = a; a = b; b = tmp;
    }
    if (a != u)
        copia_matriz(u, a, n);
}
