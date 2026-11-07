#define _USE_MATH_DEFINES
#include <math.h>
#include "calor.h"
#ifdef DIM2
/*
 solução exata:
 u(x,y,t)= e^-t sen(Kpix)sen(Kpiy)

Para:

u_t - Nabla u = f

temos:

u igual la em cima
e
-Nabla u = 2*(Kpi)^2 e^-t sen(Kpix)sen(Kpiy)

logo

f = (2K^2pi^2-1)e^-t sen(Kpix)sen(Kpiy)

condição inicial:

u(x,y,0) = sen(Kpix) sen(Kpiy)

contorno:
u=0
*/

double uex_calor(double x, double y, double t)
{
	return ( exp(-ALFA*t) * sin (K*M_PI*x) * sin (K*M_PI*y) );
}
double func_calor(double x, double y, double t)
{
	return ( ( 2.0*K*K*M_PI*M_PI - ALFA) * exp (-ALFA*t) * sin (K*M_PI*x) * sin (K*M_PI*y) );
}	
double contorno_calor(double x, double y, double t)
{
	
    	(void)x;
    	(void)y;
    	(void)t;
	return 0.0;
}
double inicial_calor(double x, double y)
{
	return uex_calor(x,y,0.0);
}
#elif defined(DIM3)
double uex_calor(double x, double y, double z, double t)
{
    return exp(-ALFA*t) * sin(K*M_PI*x) * sin(K*M_PI*y) * sin(K*M_PI*z);
}

double func_calor(double x, double y, double z, double t)
{
    return (3.0*K*K*M_PI*M_PI - ALFA)
         * exp(-ALFA*t) * sin(K*M_PI*x) * sin(K*M_PI*y) * sin(K*M_PI*z);
}

double inicial_calor(double x, double y, double z)
{
    return uex_calor(x, y, z, 0.0);
}

double contorno_calor(double x, double y, double z, double t)
{
    	(void)x;
    	(void)y;
	(void)z;
    	(void)t;
	return 0.0;
}
#else
double uex_calor(double x, double t)
{
	return exp(-ALFA*t)*sin(K*M_PI*x);
}

double func_calor(double x, double t)
{
	return (K*K*M_PI*M_PI-ALFA)*exp(-ALFA*t)*sin(K*M_PI*x);
}

double inicial_calor (double x)
{
	return uex_calor(x,0.0);
}

double contorno_calor(double x, double t)
{
	(void)x;
	(void)t;
	return 0.0;
}
#endif

