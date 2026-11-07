#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "nucleo.h"
double tempo_decorrido(struct timespec inicio, struct timespec fim)
{
    return (fim.tv_sec - inicio.tv_sec)
         + (fim.tv_nsec - inicio.tv_nsec)/1e9;
}

void imprime_ordem(double e1, double e2, double q1, double q2, const char *nome)
{
    if (e1 <= 0.0 || e2 <= 0.0 || q1 <= 0.0 || q2 <= 0.0 || q1 == q2) {
        printf("\nNao foi possivel calcular a ordem de convergencia %s.\n", nome);
        return;
    }
    double p = log(e1/e2)/log(q1/q2);
    printf("\nOrdem de convergencia %s (p): %.6f\n", nome, p);
}

