#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "nucleo.h"
int valida_malha(int n, int um_nivel)
{
    if (n < 16) {
        printf("Erro: a malha deve ser maior ou igual a 16.\n");
        return 0;
    }
    int temp = n;
    while (temp > 16) {
        if (temp % 2 != 0) {
            printf("Erro: a malha deve permitir sucessivas divisoes por 2 ate 16.\n");
            return 0;
        }
        temp /= 2;
    }
    if (temp != 16) {
        printf("Erro: a malha deve ser 16 vezes uma potencia de 2.\n");
        return 0;
    }
    if (um_nivel && n < 32) {
        printf("Erro: o modo 'v 1' exige malha pelo menos 32 para engrossar uma vez.\n");
        return 0;
    }
    return 1;
}

static int numero_niveis(int n, int um_nivel) //numero de pontos e um binario: 0, se nao é vcycle de 1 nivel; 1, c.c.
{
    int niveis = 1; //16 pontos corresponde é 1 nivel
    while (n > 16) { //cada vez que 16 puder ser divvidido por 2 é mais um nivel. ex. 32/2 = 16, logo 32 2 niveis 64/2/2 3 niveis
        n /= 2;
        niveis++;
    }
    if (um_nivel && niveis > 2) //se tiver pedido pra restringir vcycle a 1 nivel, nao pode engrossar alem disso
        niveis = 2;
    return niveis;
}

void libera_hierarquia(Hierarquia *H)
{
    if (!H) return;
    if (H->N) {
        for (int k = 0; k < H->niveis; k++) {
            if (H->r && H->r[k]) libera_matriz(H->r[k]);
            if (H->e && H->e[k]) libera_matriz(H->e[k]);
            if (H->u_novo && H->u_novo[k]) libera_matriz(H->u_novo[k]);
            if (k > 0 && H->f && H->f[k]) libera_matriz(H->f[k]);
        }
    }
    free(H->f);
    free(H->r);
    free(H->e);
    free(H->u_novo);
    free(H->N);
    memset(H, 0, sizeof(*H));
}

int aloca_hierarquia(Hierarquia *H, int n, Matriz f_fina, int um_nivel) //uma variavel pra preencher, numero de pontos, lado direito, é vcycle 1 nivel ou nao
{
    memset(H, 0, sizeof(*H)); //zera todo o H
    H->niveis = numero_niveis(n, um_nivel);  //numero de niveis (camada) da hierarquia
    H->N = malloc((size_t)H->niveis * sizeof(*H->N)); //os níveis em si. aloca memória para N niveis de H
    H->f = calloc((size_t)H->niveis, sizeof(*H->f)); //os niveis para o lado direito (comecam zerados)
    H->r = calloc((size_t)H->niveis, sizeof(*H->r)); //os niveis para o residuo
    H->e = calloc((size_t)H->niveis, sizeof(*H->e)); //para o erro
    H->u_novo = calloc((size_t)H->niveis, sizeof(*H->u_novo)); //para a nova aproximacao obtida
    if (!H->N || !H->f || !H->r || !H->e || !H->u_novo)
        goto erro;

    H->N[0] = n; //n pontos no nivel 1
    for (int k = 1; k < H->niveis; k++) //n/2 pontos no nivel 2; n/4 pontos no 3; n/8 no 4; n/16 ...
        H->N[k] = H->N[k-1]/2;

    H->f[0] = f_fina; //no nivel 1 ja foi preenchida a funçao, usa o que recebeu
    for (int k = 0; k < H->niveis; k++) {
        if (k > 0) H->f[k] = aloca_matriz(H->N[k]); //no primeiro nivel, f ja foi preenchida, a partir do segundo, aloca espaço com os pontos do nivel
        H->r[k] = aloca_matriz(H->N[k]); //aloca espaço para todos os pontos do residuo no nivel
        H->e[k] = aloca_matriz(H->N[k]); //para o erro
        H->u_novo[k] = aloca_matriz(H->N[k]); //para a solucao
        if ((k > 0 && !H->f[k]) || !H->r[k] || !H->e[k] || !H->u_novo[k])
            goto erro;
    }
    return 1;

erro:
    printf("Erro de alocacao de memoria na hierarquia multigrid.\n");
    libera_hierarquia(H);
    return 0;
}

void restringe_rhs(Hierarquia *H)
{
    for (int k = 1; k < H->niveis; k++)
        restringe(H->f[k-1], H->f[k], H->N[k-1]);//full weighting recebendo função no nivel fino, funcao agora, e numero de pontos do nivel fino (cf. "operadores.c")
}

