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
int fmg_calor(Matriz u, Hierarquia *H, double dt, double t_novo, int *iter_grosso)
{
    restringe_rhs(H);
    int nivel = H->niveis - 1;
    int n_atual = H->N[nivel];
    double h_atual = 1.0/n_atual;
    Matriz u_atual = aloca_matriz(n_atual);
    Matriz u_novo = aloca_matriz(n_atual);
    Matriz r_atual = aloca_matriz(n_atual);
    if (!u_atual || !u_novo || !r_atual) {
        libera_matriz(u_atual); libera_matriz(u_novo); libera_matriz(r_atual);
        return 0;
    }
#ifdef DEBUG
    FILE *faioul = fopen("resultados.txt", "a");

    if (faioul == NULL) {
      perror("resultados.txt");
      libera_matriz(u_atual);
      libera_matriz(u_novo);
      libera_matriz(r_atual);
      return 0;
    }
#endif
    zera_matriz(u_atual, n_atual);
    aplica_contorno_calor(u_atual, n_atual, h_atual, t_novo);
    calcula_residuo_calor(u_atual, H->f[nivel], r_atual, n_atual, h_atual, dt);
    double res = norma_residuo_calor(r_atual, n_atual, h_atual);
    *iter_grosso = 0;
    while (res > T16 && *iter_grosso < M) {
        suavizacoes_calor(u_atual, u_novo, H->f[nivel], n_atual, h_atual, dt, 1);
        calcula_residuo_calor(u_atual, H->f[nivel], r_atual, n_atual, h_atual, dt);
        res = norma_residuo_calor(r_atual, n_atual, h_atual);
#ifdef DEBUG	
	fprintf(faioul,"Resolvendo em 16x16\t Iteração %d, Resíduo: %.12e\n",(*iter_grosso), res);
#endif
        (*iter_grosso)++;
    }
    libera_matriz(u_novo);
    libera_matriz(r_atual);

    while (nivel > 0) {
        int n_fino = H->N[nivel-1];
        double h_fino = 1.0/n_fino;
        Matriz u_fino = aloca_matriz(n_fino);
        if (!u_fino) {
#ifdef DEBUG
            fclose(faioul);
#endif
            libera_matriz(u_atual);
            return 0;
        }
#ifdef DEBUG
	fprintf(faioul,"Prolongado para %d pontos\n", n_fino);
#endif
        prolonga(u_atual, u_fino, n_fino);
        aplica_contorno_calor(u_fino, n_fino, h_fino, t_novo);
        v_cycle(u_fino, H->f[nivel-1], n_fino, h_fino, dt,
                H->f, H->r, H->e, H->u_novo, nivel-1, 0);
#ifdef DEBUG
Matriz r_debug = aloca_matriz(n_fino);

if (r_debug) {
    calcula_residuo_calor(
        u_fino, H->f[nivel-1],
        r_debug, n_fino, h_fino, dt);

    double ress =
        norma_residuo_calor(r_debug, n_fino, h_fino);

    fprintf(faioul,
            "nivel %d, resíduo: %.12e\n",
            nivel, ress);

    libera_matriz(r_debug);
}
#endif
        libera_matriz(u_atual);
        u_atual = u_fino;
        nivel--;
    }
    copia_matriz(u, u_atual, H->N[0]);
    libera_matriz(u_atual);
    
#ifdef DEBUG
fclose(faioul);
#endif

    return 1;
}

