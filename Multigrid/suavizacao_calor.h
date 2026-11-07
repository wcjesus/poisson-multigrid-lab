#ifndef SUAVIZACAO_CALOR_H
#define SUAVIZACAO_CALOR_H
#include "matriz.h"

#ifdef DIM2
#define OMEGA_CALOR (4.0/5.0)
#elif defined(DIM3)
#define OMEGA_CALOR (6.0/7.0)
#else
#define OMEGA_CALOR (2.0/3.0)
#endif

void suaviza_calor(Matriz u, Matriz u_novo, Matriz f, int n, double h, double dt);
void suavizacoes_calor(Matriz u, Matriz u_novo, Matriz f, int n, double h, double dt, int numero);
#endif
