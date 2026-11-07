#ifndef RESIDUO_CALOR_H
#define RESIDUO_CALOR_H
#include "matriz.h"
void calcula_residuo_calor(Matriz u, Matriz f, Matriz r, int n, double h, double dt);
double norma_residuo_calor(Matriz r, int n, double h);
#endif
