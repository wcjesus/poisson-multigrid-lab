#ifndef RESIDUO_H
#define RESIDUO_H
#include "matriz.h"
void calcula_residuo(Matriz u, Matriz f, Matriz r, int n, double h);
double norma_residuo(Matriz r, int n, double h);
#endif
