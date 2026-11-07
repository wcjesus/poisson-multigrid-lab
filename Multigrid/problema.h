#ifndef PROBLEMA_H
#define PROBLEMA_H
#define K 4
#define EXP 2
#define _USE_MATH_DEFINES

double soma (double x);
double soma2 (double x);

#ifdef DIM2


double uex(double x, double y);

double func(double x, double y);

double contorno(double x, double y);

#elif defined(DIM3)

double uex(double x, double y, double z);

double func(double x, double y, double z);

double contorno(double x, double y, double z);

#else


double uex(double x);

double func (double x);

double contorno (double x);

#endif

#endif
