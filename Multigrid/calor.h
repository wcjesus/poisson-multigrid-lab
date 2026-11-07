#ifndef CALOR_H
#define CALOR_H
#define K 10.0
#define ALFA 4.0
#define _USE_MATH_DEFINES
#define DT_PADRAO 1e-3
#define TEMPO_FINAL 1e-1
#ifdef DIM2
double uex_calor(double x, double y, double t);

double func_calor(double x, double y, double t);

double inicial_calor(double x, double y);

double contorno_calor(double x, double y, double t);
#elif defined(DIM3)
double uex_calor(double x, double y, double z, double t);

double func_calor(double x, double y, double z, double t);

double inicial_calor(double x, double y, double z);

double contorno_calor(double x, double y, double z, double t);
#else 
double uex_calor(double x, double t);

double func_calor(double x, double t);

double inicial_calor(double x);

double contorno_calor(double x, double t);
#endif
#endif
