#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "nucleo.h"
#include "calor.h"
#include "residuo_calor.h"
#include "suavizacao_calor.h"
#include "calor_interno.h"
#include "eqcalor.h"
extern Matriz captura_solucao_espacial;
/* CALOR */
int Multigrid(int n, int modo, int um_nivel,
              double dt_solicitado, int passos_solicitados,
              double custom_tol, double *out_h, double *out_dt,
              double *out_erro_inf)
{
    if (!valida_malha(n, um_nivel)) return 1;
    if (passos_solicitados <= 0 && (!isfinite(dt_solicitado) || dt_solicitado <= 0.0
        || TEMPO_FINAL/dt_solicitado > INT_MAX)) {
        printf("Erro: dt deve ser positivo.\n");
        return 1;
    }
    double h = 1.0/n;
    double tol = custom_tol > 0.0 ? custom_tol : T;
    int passos;
    double dt;
    if (passos_solicitados > 0) {
        passos = passos_solicitados;
        dt = TEMPO_FINAL/passos;
    } else {
        passos = (int)ceil(TEMPO_FINAL/dt_solicitado);
        if (passos < 1) passos = 1;
        dt = TEMPO_FINAL/passos;
    }

    Matriz u = aloca_matriz(n), f = aloca_matriz(n), r = aloca_matriz(n);
    if (!u || !f || !r) {
        printf("Erro de alocacao de memoria.\n");
        libera_matriz(u); libera_matriz(f); libera_matriz(r);
        return 1;
    }
    inicializa_calor(u, n, h);

    Hierarquia H;
    if (!aloca_hierarquia(&H, n, f, um_nivel)) {
        libera_matriz(u); libera_matriz(f); libera_matriz(r);
        return 1;
    }
    
    #ifdef DEBUG
FILE *faioul = fopen("resultados.txt", "a");

if (faioul == NULL) {
    perror("resultados.txt");
    libera_hierarquia(&H);
    libera_matriz(u);
    libera_matriz(f);
    libera_matriz(r);
    return 1;
}
#endif

    printf("\nEQUACAO DO CALOR\n\n%s\n", modo ? "V-CYCLE MULTIGRID" : "FULL MULTIGRID");
#ifdef DIM2
    printf("Dimensao = 2D\nnx = %d\nny = %d\n", n, n);
#elif defined(DIM3)
    printf("Dimensao = 3D\nnx = %d\nny = %d\nnz = %d\n", n, n, n);
#else
    printf("Dimensao = 1D\nnx = %d\n", n);
#endif
    printf("h = %.12e\ndt = %.12e\nTempo final = %.12e\n", h, dt, (double)TEMPO_FINAL);
    printf("Passos de tempo = %d\nOmega = %.12f\n", passos, OMEGA_CALOR);
    printf("Pre-smoothing = %d\nPos-smoothing = %d\nTolerancia = %.12e\n", NU1, NU2, tol);
    if (um_nivel) printf("V-cycle limitado a um engrossamento.\n");
    printf("\n");

struct timespec inicio, fim;

double t = 0.0;
double residuo = 0.0;

long long ciclos_totais = 0;
long long iter_grosso_total = 0;

int ciclos_ultimo = 0;

#ifdef PLOT

/*
 Se passos for par, isso corresponde exatamente a T/2.
 Se for ímpar, escolhe o passo temporal mais próximo.
 */
int passo_meio = (passos + 1)/2;
int intervalo_quadros = passos/60;
if (intervalo_quadros < 1) intervalo_quadros = 1;
int numero_quadro = 0;
char sufixo[16];

/*
  u ainda é exatamente a condição inicial.
 */
salva_plot_calor(u, n, t, "inicio");
snprintf(sufixo, sizeof(sufixo), "%05d", numero_quadro++);
salva_plot_calor(u, n, t, sufixo);

#endif

clock_gettime(CLOCK_MONOTONIC, &inicio);

for (int passo = 1;
     passo <= passos;
     passo++) {
        double t_novo = passo*dt;
        monta_rhs_calor(u, f, n, h, dt, t_novo);

        if (modo == 0) {
            int itg = 0;
            if (!fmg_calor(u, &H, dt, t_novo, &itg)) {
                printf("Erro de memoria durante o FMG.\n");
                libera_hierarquia(&H); libera_matriz(u); libera_matriz(f); libera_matriz(r);
                return 1;
            }
            iter_grosso_total += itg;
        } else {
            aplica_contorno_calor(u, n, h, t_novo);
        }

        calcula_residuo_calor(u, f, r, n, h, dt);
        residuo = norma_residuo_calor(r, n, h);
#ifdef DEBUG
	fprintf(faioul,"passo %d, tempo %.9lf, resíduo: %.12e\n", passo, t_novo, residuo);
#endif
        int ciclo = 0;
        while (residuo > tol && ciclo < M) {
            v_cycle(u, f, n, h, dt, H.f, H.r, H.e, H.u_novo, 0, um_nivel);
            calcula_residuo_calor(u, f, r, n, h, dt);
            residuo = norma_residuo_calor(r, n, h);
#ifdef PRINT 
	    printf("Tempo: %.12e, Ciclo %d, resíduo: %.12e\n",t_novo, ciclo+1, residuo);
#endif 

#ifdef DEBUG
	    fprintf(faioul,"Ciclo %d, resíduo: %.12e\n", ciclo+1, residuo);
#endif
            ciclo++;
        }
        ciclos_ultimo = ciclo;
        ciclos_totais += ciclo;
        t = t_novo;
#ifdef PLOT

if (passo == passo_meio)
    salva_plot_calor(u,
                     n,
                     t,
                     "meio");
if (passo % intervalo_quadros == 0 || passo == passos) {
    snprintf(sufixo, sizeof(sufixo), "%05d", numero_quadro++);
    salva_plot_calor(u, n, t, sufixo);
}

#endif

    }
    clock_gettime(CLOCK_MONOTONIC, &fim);
#ifdef DEBUG
    fclose(faioul);
#endif

    #ifdef PLOT
salva_plot_calor(u,
                 n,
                 t,
                 "fim");
#endif
    

    double erro, erro_inf;
    calcula_erro_calor(u, n, h, t, &erro, &erro_inf);
    printf("\nRESULTADOS\n");
    printf("Tempo final atingido: %.12e\n", t);
    printf("Ciclos V-cycle no ultimo passo: %d\n", ciclos_ultimo);
    printf("Ciclos V-cycle totais: %lld\n", ciclos_totais);
    printf("Iteracoes Jacobi no nivel mais grosso: %lld\n", iter_grosso_total);
    printf("Residuo final: %.12e\n", residuo);
    printf("Norma 2 do erro: %.12e\n", erro);
    printf("Norma infinito do erro: %.12e\n", erro_inf);
    printf("Tempo de execucao: %.6f segundos\n", tempo_decorrido(inicio, fim));

    if (out_h) *out_h = h;
    if (out_dt) *out_dt = dt;
    if (out_erro_inf) *out_erro_inf = erro_inf;
    if (captura_solucao_temporal) copia_matriz(captura_solucao_temporal, u, n);
    if (captura_solucao_espacial) copia_matriz(captura_solucao_espacial, u, n);
    libera_hierarquia(&H); libera_matriz(u); libera_matriz(f); libera_matriz(r);
    return 0;
}

