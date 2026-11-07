#ifndef EQCALOR_SOLVER_H
#define EQCALOR_SOLVER_H
int Multigrid(int n, int modo, int um_nivel, double dt_solicitado, int passos_solicitados, double custom_tol, double *out_h, double *out_dt, double *out_erro_inf);
int valida_ordem_temporal_solucoes(int n, int modo, int um_nivel, int passos1);
int valida_ordem_espacial_solucoes(int n1, int n2, int n3, int modo, int um_nivel, double dt);
#endif
