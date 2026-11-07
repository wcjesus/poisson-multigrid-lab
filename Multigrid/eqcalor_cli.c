#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "nucleo.h"
#include "interface.h"
#include "eqcalor.h"
#include "calor.h"
void uso(const char *prog)
{
    printf("Equacao do calor - uso:\n");
    printf("  %s <malha> <dt>\n", prog);
    printf("  %s s <malha grossa> <malha fina> <dt>\n", prog);
    printf("  %s t <malha> <dt> <dt menor>\n", prog);
    printf("  %s n <malha> <numero de passos>\n", prog);
    printf("  %s n <malha> <passos 1> <passos 2>\n", prog);
    printf("  %s s <malha 1> <malha 2> <malha 3> <dt>\n", prog);
    printf("  %s t <malha> <dt 1> <dt 2> <dt 3>\n", prog);
    printf("  %s n <malha> <passos 1> <passos 2> <passos 3>\n", prog);
    printf("\nPrefixos opcionais:\n");
    printf("  %s v <...>       usa apenas V-cycle\n", prog);
    printf("  %s v 1 <...>     usa V-cycle limitado a um engrossamento\n", prog);
    printf("\nExemplos:\n");
    printf("  %s v s 128 256 1e-4\n", prog);
    printf("  %s v 1 t 256 1e-3 5e-4\n", prog);
#ifdef DIM2
    printf("\nExecutavel compilado para 2D.\n");
#elif defined(DIM3)
    printf("\nExecutavel compilado para 3D.\n");
#else
    printf("\nExecutavel compilado para 1D.\n");
#endif
}

int executa_argumentos(int argc, char **argv, int base, int modo, int submodo)
{
    int restantes = argc - base;

    if (restantes == 2) {
        int n = atoi(argv[base]);
        double dt = atof(argv[base+1]);
        return Multigrid(n, modo, submodo, dt, 0, -1.0, NULL, NULL, NULL);
    }

    if (restantes == 4 && strcmp(argv[base], "s") == 0) {
        int n1 = atoi(argv[base+1]);
        int n2 = atoi(argv[base+2]);
        double dt = atof(argv[base+3]);
        double h1, e1, h2, e2;
        if (Multigrid(n1, modo, submodo, dt, 0, -1.0, &h1, NULL, &e1) != 0) return 1;
        if (Multigrid(n2, modo, submodo, dt, 0, -1.0, &h2, NULL, &e2) != 0) return 1;
        imprime_ordem(e1, e2, h1, h2, "espacial");
        return 0;
    }

    if (restantes == 5 && strcmp(argv[base], "s") == 0) {
        int n1 = atoi(argv[base+1]);
        int n2 = atoi(argv[base+2]);
        int n3 = atoi(argv[base+3]);
        double dt = atof(argv[base+4]);
        return valida_ordem_espacial_solucoes(n1, n2, n3, modo, submodo, dt);
    }

    if ((restantes == 4 || restantes == 5) && strcmp(argv[base], "t") == 0) {
        int n = atoi(argv[base+1]);
        double dt_req1 = atof(argv[base+2]);
        double dt_req2 = atof(argv[base+3]);
        if (!isfinite(dt_req1) || !isfinite(dt_req2) || dt_req1 <= 0.0
            || dt_req2 <= 0.0 || TEMPO_FINAL/dt_req1 > INT_MAX/4
            || TEMPO_FINAL/dt_req2 > INT_MAX/2) {
            printf("Erro: intervalos de tempo invalidos.\n");
            return 1;
        }
        if (dt_req2 >= dt_req1) {
            printf("Erro: no modo t, o segundo dt deve ser menor que o primeiro.\n");
            return 1;
        }
        int passos1 = (int)ceil(TEMPO_FINAL/dt_req1);
        int passos2 = (int)ceil(TEMPO_FINAL/dt_req2);
        if (passos2 != 2*passos1) {
            printf("Erro: os dt efetivos devem ter razao 2, com o mesmo tempo final.\n");
            printf("Os dt pedidos resultariam em %d e %d passos; seriam necessarios %d e %d.\n",
                   passos1, passos2, passos1, 2*passos1);
            return 1;
        }
        if (restantes == 5) {
            double dt_req3 = atof(argv[base+4]);
            if (!isfinite(dt_req3) || dt_req3 <= 0.0 || dt_req3 >= dt_req2
                || TEMPO_FINAL/dt_req3 > INT_MAX) {
                printf("Erro: intervalos de tempo invalidos.\n");
                return 1;
            }
            int passos3 = (int)ceil(TEMPO_FINAL/dt_req3);
            if (passos3 != 4*passos1) {
                printf("Erro: os dt efetivos devem ter razao 2, com o mesmo tempo final.\n");
                return 1;
            }
        }
        return valida_ordem_temporal_solucoes(n, modo, submodo, passos1);
    }

    if ((restantes == 3 || restantes == 4 || restantes == 5) && strcmp(argv[base], "n") == 0) {
        int n = atoi(argv[base+1]);
        int p1 = atoi(argv[base+2]);
        if (p1 <= 0) {
            printf("Erro: o numero de passos deve ser positivo.\n");
            return 1;
        }
        if (restantes == 3)
            return Multigrid(n, modo, submodo, DT_PADRAO, p1, -1.0, NULL, NULL, NULL);

        int p2 = atoi(argv[base+3]);
        if (p2 <= 0) {
            printf("Erro: o numero de passos deve ser positivo.\n");
            return 1;
        }
        if (p1 > INT_MAX/4 || p2 != 2*p1) {
            printf("Erro: no modo n comparativo, passos 2 deve ser o dobro de passos 1.\n");
            return 1;
        }
        if (restantes == 5) {
            int p3 = atoi(argv[base+4]);
            if (p3 <= 0 || p3 != 4*p1) {
                printf("Erro: no modo n comparativo, passos 3 deve ser o quadruplo de passos 1.\n");
                return 1;
            }
        }
        return valida_ordem_temporal_solucoes(n, modo, submodo, p1);
    }

    uso(argv[0]);
    return 1;
}
