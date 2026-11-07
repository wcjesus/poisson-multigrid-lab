#ifndef CALOR_INTERNO_H
#define CALOR_INTERNO_H
#include "nucleo.h"
void aplica_contorno_calor(Matriz u, int n, double h, double t);
void inicializa_calor(Matriz u, int n, double h);
void monta_rhs_calor(Matriz u, Matriz f, int n, double h, double dt, double t_novo);
void calcula_erro_calor(Matriz u, int n, double h, double t, double *erro2, double *erroinf);
int fmg_calor(Matriz u, Hierarquia *H, double dt, double t_novo, int *iter_grosso);
extern Matriz captura_solucao_temporal;
#ifdef PLOT
void salva_plot_calor(Matriz u, int n, double t, const char *instante);
#endif
#endif
