#include <stdio.h>
#include <stdlib.h>
#include "vcycle.h"
#include "operadores.h"
#include "matriz.h"
#include "erro.h" //necessário??

#include "residuo_calor.h"
#include "suavizacao_calor.h"

/*
V-CYCLE MULTIGRID.

u = aproximacao atual
f = lado direito
nx = numero de intervalos da malha
h = espacamento
profundidade = quantas vezes o grid vai ser engrossado

   suavizacao ->
   calcula residuo ->
   restringe ->
   resolve problema coarse ->
   prolonga correcao
   corrige u ->
   suaviza

*/

void v_cycle(Matriz u, Matriz f, int n, double h, double dt,
             Matriz *f_niveis, Matriz *r_niveis,
             Matriz *e_niveis, Matriz *u_novo_niveis,
             int nivel_atual, int um_nivel)
{
    Matriz u_novo = u_novo_niveis[nivel_atual];


    /*
       Caso base normal: malha 16.
       Em "v 1", a chamada inicial esta no nivel 0 e a primeira
       malha grossa esta no nivel 1; portanto paramos ali.
    */
    if (n <= 16 || (um_nivel && nivel_atual >= 1)) {
        suavizacoes_calor(u, u_novo, f, n, h, dt, COARSE_ITER);
        return;
    }

    Matriz r = r_niveis[nivel_atual];
    int n_grossa = n/2;
    Matriz r_grossa = f_niveis[nivel_atual+1];
    Matriz e_grossa = e_niveis[nivel_atual+1];
    Matriz e_fina = u_novo;
#ifdef DEBUG
    FILE *quiv = fopen ("resultados.txt","a");
    if (!quiv){
	    perror("resultados.txt");
	    exit(EXIT_FAILURE);
	}
    setvbuf(quiv, NULL, _IOLBF, 0);
	calcula_residuo_calor(u,f,r,n,h,dt);
	fprintf(quiv, "nível %d, n = %d, resíduo atualmente (antes de NU1): %.12e \n",nivel_atual,n,norma_residuo_calor(r,n,h) );
#endif
    // Pre-smoothing - NU1 SUAVIZACOES
    suavizacoes_calor(u, u_novo, f, n, h, dt, NU1);
#ifdef DEBUG
    calcula_residuo_calor(u,f,r,n,h,dt);
    fprintf(quiv,"nivel %d n= %d residuo apoś NU1: %.12e\n",nivel_atual,n, norma_residuo_calor(r,n,h));
#endif
    // Residuo
    calcula_residuo_calor(u, f, r, n, h, dt);

#ifdef DEBUG
    fprintf(quiv,"resíduo em %d pontos: %.12e\n", n, norma_residuo_calor(r,n,h));
#endif
    // Restricao do residuo para a malha mais grossa
    restringe(r, r_grossa, n);
#ifdef DEBUG
    fprintf(quiv,"vou à malha com %d pontos resolver A_H e_H = Rr\n", n/2);
#endif
    zera_matriz(e_grossa, n_grossa);

    /*
       4/5. resolve recursivamente o problema do erro na malha grossa.
       Na malha grossa H = 2h. O dt permanece o mesmo.
    */
#ifdef DEBUG
    
    fprintf(quiv,"chamada recursiva de vciclo pra resolver correcao\n");
    
#endif

    v_cycle(e_grossa, r_grossa, n_grossa, 2.0*h, dt,
            f_niveis, r_niveis, e_niveis, u_novo_niveis,
            nivel_atual+1, um_nivel);
#ifdef DEBUG
    calcula_residuo_calor(e_grossa, r_grossa, r_niveis[nivel_atual + 1], n_grossa, 2.0*h, dt);
    fprintf(quiv, "[nivel %d, n=%d] residuo A_H e_H = Rr: %.12e\n",
            nivel_atual + 1, n_grossa, norma_residuo_calor(r_niveis[nivel_atual + 1], n_grossa, 2.0*h));
    fprintf(quiv,"Resolvi a correcao já, vou prolongar para a fina\n");
#endif
    //Prolongacao da correcao para a malha fina 
    prolonga(e_grossa, e_fina, n);

    //Correcao apenas no interior. o erro tem contorno zero. */
#ifdef DIM2
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            u[i][j] += e_fina[i][j];
#elif defined(DIM3)
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            for (int k = 1; k < n; k++)
                u[i][j][k] += e_fina[i][j][k];
#else
    for (int i = 1; i < n; i++)
	    u[i] += e_fina[i];
#endif
#ifdef DEBUG
    calcula_residuo_calor(u, f, r, n, h, dt);
    fprintf(quiv, "[nivel %d, n=%d] residuo apos correcao: %.12e\n",
            nivel_atual, n, norma_residuo_calor(r, n, h));
#endif
    //Possmoothing 
    suavizacoes_calor(u, u_novo, f, n, h, dt, NU2);
#ifdef DEBUG
    calcula_residuo_calor(u, f, r, n, h, dt);
    fprintf(quiv, "[nivel %d, n=%d] residuo apos NU2: %.12e\n",
            nivel_atual, n, norma_residuo_calor(r, n, h));
    fclose(quiv);
#endif
}
