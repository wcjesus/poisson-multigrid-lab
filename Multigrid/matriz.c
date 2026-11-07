#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matriz.h"

Matriz aloca_matriz(int n)
{
    size_t lado = (size_t)n + 1;

#ifdef DIM2
    double **A = malloc(lado * sizeof(*A));
    double *dados = calloc(lado * lado, sizeof(*dados));
    if (!A || !dados) {
        free(A);
        free(dados);
        return NULL;
    }
    for (int i = 0; i <= n; i++)
        A[i] = dados + (size_t)i * lado;
    return A;
#elif defined(DIM3)
    double ***A = malloc(lado * sizeof(*A));
    double **planos = malloc(lado * lado * sizeof(*planos));
    double *dados = calloc(lado * lado * lado, sizeof(*dados));
    if (!A || !planos || !dados) {
        free(A);
        free(planos);
        free(dados);
        return NULL;
    }
    for (int i = 0; i <= n; i++) {
        A[i] = planos + (size_t)i * lado;
        for (int j = 0; j <= n; j++)
            A[i][j] = dados + ((size_t)i * lado + (size_t)j) * lado;
    }
    return A;
#else
    double *A = calloc (lado, sizeof(*A));
    if (!A) return NULL;

    return A;
    
#endif
}

void libera_matriz(Matriz A)
{
    if (!A) return;
#ifdef DIM2
    free(A[0]);
    free(A);
#elif defined(DIM3)
    free(A[0][0]);
    free(A[0]);
    free(A);
#else
    free (A);
#endif
}

void zera_matriz(Matriz A, int n)
{
    if (!A) return;
    size_t lado = (size_t)n + 1; //quantos bytes preencher? n+1 colunas
#ifdef DIM2
    memset(A[0], 0, lado * lado * sizeof(double)); //preenchendo (n+1)linhas*(n+1)colunas de A com 0
#elif defined(DIM3)
    memset(A[0][0], 0, lado * lado * lado * sizeof(double)); //preenchendo (n+1)linhas(n+1)colunas(n+1)cotas com 0
#else 
    memset(A, 0, lado*sizeof(double)); //preenchendo cada posição com 0
#endif
}

void copia_matriz(Matriz destino, Matriz origem, int n)
{
    size_t lado = (size_t)n + 1;
#ifdef DIM2
    memcpy(destino[0], origem[0], lado * lado * sizeof(double));
#elif defined(DIM3)
    memcpy(destino[0][0], origem[0][0], lado * lado * lado * sizeof(double));
#else
    memcpy(destino, origem, lado*sizeof(double));
#endif
}
