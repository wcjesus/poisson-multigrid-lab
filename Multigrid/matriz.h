#ifndef MATRIZ_H
#define MATRIZ_H

#ifdef DIM2
typedef double **Matriz;
#elif DIM3
typedef double ***Matriz;
#else 
typedef double *Matriz;
#endif

Matriz aloca_matriz(int n);
void libera_matriz(Matriz A);
void zera_matriz(Matriz A, int n);
void copia_matriz(Matriz destino, Matriz origem, int n);

#endif

