#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "nucleo.h"
#include "calor.h"
#include "eqcalor.h"

Matriz captura_solucao_espacial = NULL;

static double diferenca_malhas(Matriz grossa, Matriz fina, int ng, int nf)
{
    int fator = nf/ng;
    double resultado = 0.0;
#ifdef DIM2
    for (int i = 0; i <= ng; i++)
        for (int j = 0; j <= ng; j++) {
            double d = fabs(grossa[i][j] - fina[fator*i][fator*j]);
            if (!isfinite(d)) return NAN;
            if (d > resultado) resultado = d;
        }
#elif defined(DIM3)
    for (int i = 0; i <= ng; i++)
        for (int j = 0; j <= ng; j++)
            for (int k = 0; k <= ng; k++) {
                double d = fabs(grossa[i][j][k] - fina[fator*i][fator*j][fator*k]);
                if (!isfinite(d)) return NAN;
                if (d > resultado) resultado = d;
            }
#else
    for (int i = 0; i <= ng; i++) {
        double d = fabs(grossa[i] - fina[fator*i]);
        if (!isfinite(d)) return NAN;
        if (d > resultado) resultado = d;
    }
#endif
    return resultado;
}

int valida_ordem_espacial_solucoes(int n1, int n2, int n3, int modo, int um_nivel, double dt)
{
    if (!valida_malha(n1, um_nivel) || !valida_malha(n2, um_nivel)
        || !valida_malha(n3, um_nivel)) return 1;
    if (n2 <= n1 || n3 <= n2 || n2 % n1 != 0 || n3 % n2 != 0
        || n2/n1 != n3/n2) {
        printf("Erro: no modo s com tres malhas, as malhas devem ter razao de refinamento constante.\n");
        return 1;
    }
    if (!isfinite(dt) || dt <= 0.0 || TEMPO_FINAL/dt > INT_MAX) {
        printf("Erro: dt deve ser positivo e gerar numero valido de passos.\n");
        return 1;
    }

    Matriz u1 = aloca_matriz(n1), u2 = aloca_matriz(n2), u3 = aloca_matriz(n3);
    if (!u1 || !u2 || !u3) {
        printf("Erro de alocacao na comparacao espacial.\n");
        libera_matriz(u1); libera_matriz(u2); libera_matriz(u3);
        return 1;
    }
    int status;
    captura_solucao_espacial = u1;
    status = Multigrid(n1, modo, um_nivel, dt, 0, -1.0, NULL, NULL, NULL);
    if (status) goto fim;
    captura_solucao_espacial = u2;
    status = Multigrid(n2, modo, um_nivel, dt, 0, -1.0, NULL, NULL, NULL);
    if (status) goto fim;
    captura_solucao_espacial = u3;
    status = Multigrid(n3, modo, um_nivel, dt, 0, -1.0, NULL, NULL, NULL);
    if (status) goto fim;

    double dif12 = diferenca_malhas(u1, u2, n1, n2);
    double dif23 = diferenca_malhas(u2, u3, n2, n3);
    printf("\nDiferenca infinito entre solucoes (h1, h2): %.12e\n", dif12);
    printf("Diferenca infinito entre solucoes (h2, h3): %.12e\n", dif23);
    if (isfinite(dif12) && isfinite(dif23) && dif12 > 0.0 && dif23 > 0.0)
        printf("Ordem de convergencia espacial por tres malhas (p): %.6f\n", log(dif12/dif23)/log((double)n2/n1));
    else
        printf("Ordem de convergencia espacial por tres malhas (p): N/A (diferencas nulas ou invalidas)\n");
    status = 0;

fim:
    captura_solucao_espacial = NULL;
    libera_matriz(u1); libera_matriz(u2); libera_matriz(u3);
    return status;
}
