#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "nucleo.h"
#include "interface.h"
#include "poisson.h"
void uso(const char *prog)
{
    printf("Poisson - uso:\n");
    printf("  %s <malha>\n", prog);
    printf("  %s <malha> <malha fina>\n", prog);
    printf("  %s 0 <malha> <expoente da tolerancia>\n", prog);
    printf("\nPrefixos opcionais:\n");
    printf("  %s v <...>       usa apenas V-cycle\n", prog);
    printf("  %s v 1 <...>     usa V-cycle limitado a um engrossamento\n", prog);
    printf("\nExemplo: %s v 1 1024 2048\n", prog);
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
    if (restantes == 1) {
        int n = atoi(argv[base]);
        return Multigrid(n, modo, submodo, -1.0, NULL, NULL);
    }
    if (restantes == 2) {
        int n1 = atoi(argv[base]);
        int n2 = atoi(argv[base+1]);
        double h1, e1, h2, e2;
        if (Multigrid(n1, modo, submodo, -1.0, &h1, &e1) != 0) return 1;
        if (Multigrid(n2, modo, submodo, -1.0, &h2, &e2) != 0) return 1;
        imprime_ordem(e1, e2, h1, h2, "espacial");
        return 0;
    }
    if (restantes == 3 && strcmp(argv[base], "0") == 0) {
        int n = atoi(argv[base+1]);
        int exp = atoi(argv[base+2]);
        if (exp < 0) {
            printf("Erro: o expoente da tolerancia deve ser nao negativo.\n");
            return 1;
        }
        double tol = 1.0;
        for (int k = 0; k < exp; k++) tol /= 10.0;
        return Multigrid(n, modo, submodo, tol, NULL, NULL);
    }
    if (restantes == 4 && strcmp(argv[base], "0") == 0) {
	    int n1 = atoi(argv[base+1]);
	    int n2 = atoi(argv[base+2]);
	    int exp = atoi(argv[base+3]);
        if (exp < 0) {
            printf("Erro: o expoente da tolerancia deve ser nao negativo.\n");
            return 1;
        }
	    double tol = 1.0;
	    double h1,e1,h2,e2;
	    for (int k=0;k<exp;k++) tol /= 10.0;
	    if (Multigrid(n1,modo,submodo,tol,&h1,&e1) !=0 ) return 1;
	    if (Multigrid(n2,modo,submodo,tol,&h2,&e2) !=0 ) return 1;
	    imprime_ordem(e1,e2,h1,h2, "espacial");
        return 0;
    }
    uso(argv[0]);
    return 1;
}
