#include <math.h>
#include "problema.h"
#define _USE_MATH_DEFINES

double soma (double x)
{
	double s=0.0;
	int freq = K;

	for (int j=1; j<=EXP; j++)
	{
		s += sin(freq*M_PI*x);
		freq*=K;
	}
	return s;
}

double soma2 (double x)
{
	double s=0.0;
	int freq =K;

	for (int i=1; i<=EXP; i++)
	{
		s+= (double)freq*freq*sin(freq*M_PI*x);
		freq *=K;
	}
	return s;
}

#ifdef DIM2
/*
Solucao exata:

   u(x,y) = sen(K*pi*x) sen(K*pi*y)

Para:

   -Delta u = f

temos:

   f(x,y) = 2*(K*pi)^2*sen(K*pi*x)sen(K*pi*y)

*/

double uex(double x, double y)
{
return (soma(x)*soma(y));
}

double func(double x, double y)
{
return M_PI * M_PI*(soma2(x)*soma(y)+soma2(y)*soma(x));
}

double contorno(double x, double y)
{
	return 0.0*x*y;
}

#elif defined(DIM3)
double uex(double x, double y, double z)
{
    return soma(x)*soma(y)*soma(z);
}

double func(double x, double y, double z)
{
    return M_PI*M_PI*(soma2(x)*soma(y)*soma(z)+soma(x)*soma2(y)*soma(z)+soma(x)*soma(y)*soma2(z));
}

double contorno(double x, double y, double z)
{
    return 0.0*x*y*z;
}

#else
double uex(double x)
{
	return soma(x);
}

double func(double x)
{
	return (soma2(x)*M_PI*M_PI);
}

double contorno (double x)
{
	return 0.0*x;
}

#endif
