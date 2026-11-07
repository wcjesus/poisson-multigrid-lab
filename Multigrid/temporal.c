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
/* Armazena a solucao final somente durante a validacao temporal por tres niveis. */
Matriz captura_solucao_temporal = NULL;

static double diferenca(Matriz u1, Matriz u2, int n)
{
    double resultado = 0.0;
#ifdef DIM2
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= n; j++) {
            double d = fabs(u1[i][j] - u2[i][j]);
            if (!isfinite(d)) return NAN;
            if (d > resultado) resultado = d;
        }
#elif defined(DIM3)
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= n; j++)
            for (int k = 0; k <= n; k++) {
                double d = fabs(u1[i][j][k] - u2[i][j][k]);
                if (!isfinite(d)) return NAN;
                if (d > resultado) resultado = d;
            }
#else
    for (int i = 0; i <= n; i++) {
        double d = fabs(u1[i] - u2[i]);
        if (!isfinite(d)) return NAN;
        if (d > resultado) resultado = d;
    }
#endif
    return resultado;
}

/* Usa duas matrizes auxiliares: apos comparar o primeiro par, reutiliza a
   memoria da solucao mais grossa para a solucao mais fina. */
int valida_ordem_temporal_solucoes(int n, int modo, int um_nivel,
                                         int passos1)
{
    if (!valida_malha(n, um_nivel)) return 1;
    if (passos1 <= 0 || passos1 > INT_MAX/4) {
        printf("Erro: numero de passos fora do intervalo permitido.\n");
        return 1;
    }

    Matriz u1 = aloca_matriz(n);
    Matriz u2 = aloca_matriz(n);
    if (!u1 || !u2) {
        printf("Erro de alocacao na comparacao temporal.\n");
        libera_matriz(u1); libera_matriz(u2);
        return 1;
    }

    int status;
    captura_solucao_temporal = u1;
    status = Multigrid(n, modo, um_nivel, DT_PADRAO, passos1,
                       -1.0, NULL, NULL, NULL);
    if (status) goto fim;

    captura_solucao_temporal = u2;
    status = Multigrid(n, modo, um_nivel, DT_PADRAO, 2*passos1,
                       -1.0, NULL, NULL, NULL);
    if (status) goto fim;

    double diferenca_grossa_media = diferenca(u1, u2, n);

    captura_solucao_temporal = u1;
    status = Multigrid(n, modo, um_nivel, DT_PADRAO, 4*passos1,
                       -1.0, NULL, NULL, NULL);
    if (status) goto fim;

    double diferenca_media_fina = diferenca(u2, u1, n);
    printf("\nDiferenca infinito entre solucoes (dt1, dt2): %.12e\n",
           diferenca_grossa_media);
    printf("Diferenca infinito entre solucoes (dt2, dt3): %.12e\n",
           diferenca_media_fina);
    if (isfinite(diferenca_grossa_media) && isfinite(diferenca_media_fina)
        && diferenca_grossa_media > 0.0 && diferenca_media_fina > 0.0) {
        double p = log(diferenca_grossa_media/diferenca_media_fina)/log(2.0);
        printf("Ordem de convergencia temporal (p): %.6f\n", p);
    } else {
        printf("Ordem de convergencia temporal (p): N/A (diferencas nulas ou invalidas)\n");
    }
    status = 0;

fim:
    captura_solucao_temporal = NULL;
    libera_matriz(u1);
    libera_matriz(u2);
    return status;
}
