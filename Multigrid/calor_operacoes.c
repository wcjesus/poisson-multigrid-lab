#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "nucleo.h"
#include "calor.h"
#include "residuo_calor.h"
#include "suavizacao_calor.h"
#include "calor_interno.h"
void aplica_contorno_calor(Matriz u, int n, double h, double t)
{
#ifdef DIM2
    for (int i = 0; i <= n; i++) {
        u[i][0] = contorno_calor(i*h, 0.0, t);
        u[i][n] = contorno_calor(i*h, 1.0, t);
    }
    for (int j = 0; j <= n; j++) {
        u[0][j] = contorno_calor(0.0, j*h, t);
        u[n][j] = contorno_calor(1.0, j*h, t);
    }
#elif defined(DIM3)
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= n; j++) {
            u[i][j][0] = contorno_calor(i*h, j*h, 0.0, t);
            u[i][j][n] = contorno_calor(i*h, j*h, 1.0, t);
        }
    for (int i = 0; i <= n; i++)
        for (int k = 0; k <= n; k++) {
            u[i][0][k] = contorno_calor(i*h, 0.0, k*h, t);
            u[i][n][k] = contorno_calor(i*h, 1.0, k*h, t);
        }
    for (int j = 0; j <= n; j++)
        for (int k = 0; k <= n; k++) {
            u[0][j][k] = contorno_calor(0.0, j*h, k*h, t);
            u[n][j][k] = contorno_calor(1.0, j*h, k*h, t);
        }
#else
    (void)h;
    u[0] = contorno_calor(0.0, t);
    u[n] = contorno_calor(1.0, t);
#endif
}

void inicializa_calor(Matriz u, int n, double h)
{
#ifdef DIM2
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= n; j++)
            u[i][j] = inicial_calor(i*h, j*h);
#elif defined(DIM3)
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= n; j++)
            for (int k = 0; k <= n; k++)
                u[i][j][k] = inicial_calor(i*h, j*h, k*h);
#else 
    for (int i =0; i<=n; i++)
	    u[i] = inicial_calor(i*h);
#endif
}

void monta_rhs_calor(Matriz u, Matriz f, int n, double h, double dt, double t_novo)
{
    zera_matriz(f, n);
#ifdef DIM2
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            f[i][j] = u[i][j] + dt*func_calor(i*h, j*h, t_novo);
#elif defined(DIM3)
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            for (int k = 1; k < n; k++)
                f[i][j][k] = u[i][j][k]
                           + dt*func_calor(i*h, j*h, k*h, t_novo);
#else
    for (int i =1; i<n;i++)
	    f[i] = u[i] + dt*func_calor(i*h, t_novo);
#endif
}

void calcula_erro_calor(Matriz u, int n, double h, double t,
                               double *erro2, double *erroinf)
{
    double soma = 0.0, maxi = 0.0;
#ifdef DIM2
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++) {
            double e = u[i][j] - uex_calor(i*h, j*h, t);
            double ae = fabs(e);
            soma += e*e;
            if (ae > maxi) maxi = ae;
        }
    *erro2 = sqrt(h*h*soma);
#elif defined(DIM3)
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            for (int k = 1; k < n; k++) {
                double e = u[i][j][k] - uex_calor(i*h, j*h, k*h, t);
                double ae = fabs(e);
                soma += e*e;
                if (ae > maxi) maxi = ae;
            }
    *erro2 = sqrt(h*h*h*soma);
#else
    for (int i=1;i<n;i++){
	    double e = u[i] - uex_calor(i*h,t);
	    double ae = fabs(e);
	    soma += e*e;
	    if (ae>maxi) maxi = ae;
    }
    *erro2 = sqrt(h*soma);
#endif
    *erroinf = maxi;
}

