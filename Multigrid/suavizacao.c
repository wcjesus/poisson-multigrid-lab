#include "matriz.h"
#include "suavizacao.h"

/*
Jacobi ponderado.

Primeiro calculamos o Jacobi normal:

   u_J =
   1/4 * (h^2*f
          + vizinhos)

Depois fazemos:

   u_novo =
   (1-omega)*u + omega*u_J

com omega = 2/3
*/

void suaviza(Matriz u, Matriz u_novo, Matriz f, int n, double h)
{
    double h2 = h*h;
    copia_matriz(u_novo, u, n);
#ifdef DIM2
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++) {
            double uj = 0.25*(h2*f[i][j] 
			    + u[i-1][j] + u[i+1][j]
                            + u[i][j-1] + u[i][j+1]);
            u_novo[i][j] = (1.0-OMEGA)*u[i][j] + OMEGA*uj;
        }
#elif defined(DIM3)
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            for (int k = 1; k < n; k++) {
                double uj = (h2*f[i][j][k]
                           + u[i-1][j][k] + u[i+1][j][k]
                           + u[i][j-1][k] + u[i][j+1][k]
                           + u[i][j][k-1] + u[i][j][k+1]) / 6.0;
                u_novo[i][j][k] = (1.0-OMEGA)*u[i][j][k] + OMEGA*uj;
            }
#else
    for (int i = 1; i < n; i++) {
	    double uj = (h2*f[i] 
			+ u[i-1] + u[i+1]) / 2.0;
	    u_novo[i] = (1.0-OMEGA)*u[i] + OMEGA*uj;
    }    
#endif
}

void suavizacoes(Matriz u, Matriz u_novo, Matriz f, int n, double h, int numero)
{
    Matriz a = u;
    Matriz b = u_novo;
    for (int it = 0; it < numero; it++) {
        suaviza(a, b, f, n, h); //b vai ser o a (que esta com u) apos aplicação de jacobi pond. 
        Matriz tmp = a; a = b; b = tmp;    /*a (#) b(@) tmp (!) -> (temp recebe a) -> a(#) b(@) tmp (#) -> (a recebe b) -> a(@) b(@) tmp(#)
                                             -> (b recebe tmp) -> a(@) b(#) tmp (#)*/
}
    if (a != u)
        copia_matriz(u, a, n); //a ficou com u_novo que estava com b, agora u vai virar u_novo (copiado de a)
}
