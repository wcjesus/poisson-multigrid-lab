#ifndef NUCLEO_H
#define NUCLEO_H
#include "matriz.h"
#include "vcycle.h"
#include "operadores.h"
#include <time.h>
#define M 100000
#define T 1e-7
#define T16 1e-12
#define CONSTANTE 0.01

typedef struct {
    int niveis;
    int *N;
    Matriz *f;
    Matriz *r;
    Matriz *e;
    Matriz *u_novo;
} Hierarquia;

int valida_malha(int n, int um_nivel);
void libera_hierarquia(Hierarquia *H);
int aloca_hierarquia(Hierarquia *H, int n, Matriz f_fina, int um_nivel);
void restringe_rhs(Hierarquia *H);
double tempo_decorrido(struct timespec inicio, struct timespec fim);
void imprime_ordem(double e1, double e2, double q1, double q2, const char *nome);
#endif
