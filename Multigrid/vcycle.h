#ifndef VCYCLE_H
#define VCYCLE_H
#include "matriz.h"
#define COARSE_ITER 500
#define NU1 2
#define NU2 2
void v_cycle(Matriz u, Matriz f, int n, double h, double dt,
             Matriz *f_niveis, Matriz *r_niveis,
             Matriz *e_niveis, Matriz *u_novo_niveis,
             int nivel_atual, int um_nivel);
#endif
