#ifndef SUAVIZACAO_H
#define SUAVIZACAO_H
#include "matriz.h"
#ifdef DIM2
#define OMEGA (4.0/5.0)
#elif defined(DIM3)
#define OMEGA (6.0/7.0)
#else
#define OMEGA (2.0/3.0)
#endif
void suaviza(Matriz u, Matriz u_novo, Matriz f, int n, double h);
void suavizacoes(Matriz u, Matriz u_novo, Matriz f, int n, double h, int numero);
#endif
