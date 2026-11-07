#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "nucleo.h"
#include "problema.h"
#include "residuo.h"
#include "erro.h"
#include "suavizacao.h"
#include "poisson.h"
#include "poisson_plot.h"
static void aplica_contorno_poisson(Matriz u, int n, double h)  //matriz da aproximação, pontos, passo
{
	//preenche as paredes (bordo) com dirichlet (contorno está em "problema.c")
#ifdef DIM2
    for (int i = 0; i <= n; i++) {
        u[i][0] = contorno(i*h, 0.0); //parede de baixo
        u[i][n] = contorno(i*h, 1.0); //parede de cima
    }
    for (int j = 0; j <= n; j++) {
        u[0][j] = contorno(0.0, j*h); //parede esq.
        u[n][j] = contorno(1.0, j*h); //parede dir.
    }
#elif defined(DIM3)
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= n; j++) {
            u[i][j][0] = contorno(i*h, j*h, 0.0); //plano de baixo (fundo da caixa)
            u[i][j][n] = contorno(i*h, j*h, 1.0); //plano de cima (tampa da caixa)
        }
    for (int i = 0; i <= n; i++)
        for (int k = 0; k <= n; k++) {
            u[i][0][k] = contorno(i*h, 0.0, k*h); //plano esquerdo
            u[i][n][k] = contorno(i*h, 1.0, k*h); //plano direito
        }
    for (int j = 0; j <= n; j++)
        for (int k = 0; k <= n; k++) {
            u[0][j][k] = contorno(0.0, j*h, k*h); //plano de trás
            u[n][j][k] = contorno(1.0, j*h, k*h); //plano da frente
        }
#else 
    for (int i=0;i<=n;i++) {
	u[0] = contorno(0.0*h); //inicio da linha
	u[n] = contorno(1.0); //fim da linha
    }
#endif
}

static void monta_poisson(Matriz u, Matriz f, int n, double h) //recebe matriz u da aproximação, f do lado direito, pontos da malha, passo da malha
{
    zera_matriz(u, n); /*zera as n posições de u para usar o chute inicial (está em "matriz")*/
#ifdef DIM2
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= n; j++)
            f[i][j] = func(i*h, j*h); //preenche o lado direito com os valores da função manufaturada em cada ponto (cada ponto é índice do ponto*passo)
#elif defined(DIM3)
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= n; j++)
            for (int k = 0; k <= n; k++)
                f[i][j][k] = func(i*h, j*h, k*h);
#else 
    for (int i=0; i<=n; i++)
	    f[i] = func(i*h);
#endif
    aplica_contorno_poisson(u, n, h); //preenche bordos de u com poisson
}

static int fmg_poisson(Matriz u, Hierarquia *H, int *iter_grosso)
{
    restringe_rhs(H); //passa f no nivel atual para preencher f restringida a n/2 pontos
    int nivel = H->niveis - 1;
    int n_atual = H->N[nivel]; //numero de pontos atual é o numero de pontos no nivel atual
    double h_atual = 1.0/n_atual; //passo atual
    Matriz u_atual = aloca_matriz(n_atual); //espaco para n_atual pontos
    Matriz u_novo = aloca_matriz(n_atual);
    Matriz r_atual = aloca_matriz(n_atual);
    if (!u_atual || !u_novo || !r_atual) {
        libera_matriz(u_atual); libera_matriz(u_novo); libera_matriz(r_atual);
        return 0;
    }
#ifdef DEBUG
    FILE* faioul = fopen ("resultados.txt", "a");
    if (!faioul) {
        perror("resultados.txt");
        libera_matriz(u_atual);
        libera_matriz(u_novo);
        libera_matriz(r_atual);
        return 0;
    }
#endif
    zera_matriz(u_atual, n_atual); //u_atual ainda vai ser preenchido, zeramos
    aplica_contorno_poisson(u_atual, n_atual, h_atual); //preenche os bordos
    calcula_residuo(u_atual, H->f[nivel], r_atual, n_atual, h_atual); /*passa o u pra preencher, a f no nivel pra calcular Au-f, o residuo para ser atualizado recebendo Au-f, pontos no nível atual, passo*/
    double res = norma_residuo(r_atual, n_atual, h_atual); //agora que o residuo foi preenchido, passa ele para calcular norma euclidiana (cf. residuo.c)
    *iter_grosso = 0;
#ifdef PRINT
    	printf("Resolvendo em 16 pontos\n");
#endif
    while (res > T16 && *iter_grosso < M) {
        suavizacoes(u_atual, u_novo, H->f[nivel], n_atual, h_atual, 1); //UMA suavização (jacobi pond.) e o valor novo fica em u_atual
        calcula_residuo(u_atual, H->f[nivel], r_atual, n_atual, h_atual);
        res = norma_residuo(r_atual, n_atual, h_atual);
        (*iter_grosso)++;
    }
#ifdef PRINT
    printf("%d iterações em 16 pontos\n", (*iter_grosso));
#endif
    libera_matriz(u_novo); //nao precisa mais, é u_atual que vira o chute inicial do vcycle
    libera_matriz(r_atual);

    while (nivel > 0) {
        int n_fino = H->N[nivel-1]; //usa o numero de pontos do nivel anterior para ser o numero de pontos da malha fina
        double h_fino = 1.0/n_fino;
        Matriz u_fino = aloca_matriz(n_fino);
        if (!u_fino) {
#ifdef DEBUG
            fclose(faioul);
#endif
            libera_matriz(u_atual);
            return 0;
        }
#ifdef DEBUG
	fprintf(faioul,"Estamos no nível %d com %d pontos\n", nivel,H->N[nivel]);
#endif
        prolonga(u_atual, u_fino, n_fino);//interpolação (cf. operadores.h) para obbter erro na fina
#ifdef DEBUG
	fprintf(faioul,"Prolongado u para %d pontos\n", n_fino);
#endif	
        aplica_contorno_poisson(u_fino, n_fino, h_fino);
        v_cycle(u_fino, H->f[nivel-1], n_fino, h_fino, 0.0,
                H->f, H->r, H->e, H->u_novo, nivel-1, 0); /*/uu na fina, f nafina, pontos da fina, passo fino, dt =0, f aqui, r aqui, e aqui, u_novo, nivel atual, resposta: não é vcycle de um nivel*/
#ifdef DEBUG
	
	fprintf(faioul,"concluí 1 vciclo\n");
	Matriz r_debug = aloca_matriz(n_fino); //espaço pro residuo
    if (r_debug) {
        calcula_residuo(u_fino, H->f[nivel-1], r_debug, n_fino, h_fino);
        fprintf(faioul, "residuo apos V-cycle (n=%d): %.12e\n",
                n_fino, norma_residuo(r_debug, n_fino, h_fino));//como está o residuo atualmente
        libera_matriz(r_debug);
    } else {
        fprintf(faioul, "falha ao alocar residuo para DEBUG\n");
    }
							  
#endif
        libera_matriz(u_atual);
        u_atual = u_fino; //atualiza u para apos o v ciclo
        nivel--;
    }

    copia_matriz(u, u_atual, H->N[0]); //passa u_atual pra u em pontos no nivel mais fino
    libera_matriz(u_atual);
#ifdef DEBUG
    fprintf(faioul,"A partir de agora é vciclo apenas\n");
    fclose(faioul);
#endif
    return 1;
}

int Multigrid(int n, int modo, int um_nivel, double custom_tol,
              double *out_h, double *out_erro_inf)
{
    if (!valida_malha(n, um_nivel)) return 1; //verifica se foi passado um número de pontos com o qual se pode trabalhar

    double h = 1.0/n; //calcula o passo da malha
    double tol = custom_tol > 0.0 ? custom_tol : T; //se passou tolerancia personalizada, tol é esse valor passado. c.c. tol é o valor base definido

    Matriz u = aloca_matriz(n), f = aloca_matriz(n), r = aloca_matriz(n); //aloca espaço para a aproximação, o lado direito (f) e o resíduo  r = Av - f
    if (!u || !f || !r) {
        printf("Erro de alocacao de memoria.\n");
        libera_matriz(u); libera_matriz(f); libera_matriz(r);
        return 1;
    }
    monta_poisson(u, f, n, h);//acima. zera as posições de u e depois preenche os bordos com dirichlet. preenche as posiçoes de f com a sol.manufaturada

    Hierarquia H; //struct com numero de niveis, pontos por nivel, lado direito em cada nivel, residuo emc cada nivel, erro em cada nivel e aproximacao emc ada nivel
    if (!aloca_hierarquia(&H, n, f, um_nivel)) {
        libera_matriz(u); libera_matriz(f); libera_matriz(r);
        return 1;
    }

    printf("\nPROBLEMA DE POISSON\n\n%s \n", um_nivel ? "TWO GRID" :  modo ? "Multigrid V-Cycle"  : "Full Multigrid");
#ifdef DIM2
    printf("Dimensao = 2D\nnx = %d\nny = %d\n", n, n);
#elif defined(DIM3)
    printf("Dimensao = 3D\nnx = %d\nny = %d\nnz = %d\n", n, n, n);
#else 
    printf("Dimensao = 1D\nnx = %d\n", n);
#endif
    printf("h = %.12e\nOmega = %.12f\n", h, OMEGA);
    printf("Pre-smoothing = %d\nPos-smoothing = %d\nTolerancia = %.12e\n", NU1, NU2, tol);
    if (um_nivel) printf("V-cycle limitado a um engrossamento.\n");
    printf("\n");

    struct timespec inicio, fim;
    clock_gettime(CLOCK_MONOTONIC, &inicio); //nao contamos o tempo pra alocar espaço, começamos agora a contar 
    int iter_grosso = 0;
    if (modo == 0) { //se for full multigrid...
        if (!fmg_poisson(u, &H, &iter_grosso)) {
            printf("Erro de memoria durante o FMG.\n");
            libera_hierarquia(&H); libera_matriz(u); libera_matriz(f); libera_matriz(r);
            return 1;
        }
    }

#ifdef DEBUG
    FILE* faioul = fopen("resultados.txt", "a");
    if (!faioul) {
        perror("resultados.txt");
        libera_hierarquia(&H);
        libera_matriz(u); libera_matriz(f); libera_matriz(r);
        return 1;
    }
#endif
    calcula_residuo(u, f, r, n, h);//um ciclo o if fez, olha o residuo pra ver se faz mais
    double residuo = norma_residuo(r, n, h);
    int ciclo = 0;
    while (residuo > tol && ciclo < M) {
        v_cycle(u, f, n, h, 0.0, H.f, H.r, H.e, H.u_novo, 0, um_nivel);
        calcula_residuo(u, f, r, n, h);
        residuo = norma_residuo(r, n, h);
#ifdef PRINT
	printf("Ciclo %d: Resíduo: %.12e\n", ciclo +1, residuo);
#endif
#ifdef DEBUG	
	fprintf(faioul, "Ciclo %d: Resíduo: %.12e\n", ciclo+1,residuo);
#endif
        ciclo++;
    }
    clock_gettime(CLOCK_MONOTONIC, &fim); //nao precisa contar o tempo de calular o erro
#ifdef DEBUG
    fclose(faioul);
#endif
#ifdef PLOT
    salva_plot_poisson(u, n);
#endif
    double erro = norma_erro(u, n, h);
    double erro_inf = norma_infinito_erro(u, n, h);
    printf("\nRESULTADOS\n");
    printf("Ciclos V-cycle no nivel mais fino: %d\n", ciclo);
    printf("Iteracoes Jacobi no nivel mais grosso: %d\n", iter_grosso);
    printf("Residuo final: %.12e\n", residuo);
    printf("Norma 2 do erro: %.12e\n", erro);
    printf("Norma infinito do erro: %.12e\n", erro_inf);
    printf("Tempo de execucao: %.6f segundos\n", tempo_decorrido(inicio, fim));
    printf("%s\n", residuo < tol ? "Convergencia atingida." : "Maximo de ciclos atingido.");

    if (out_h) *out_h = h; //se foram passadas duas malhas ele vai trocar os valores para ccalcular o p
    if (out_erro_inf) *out_erro_inf = erro_inf;
    libera_hierarquia(&H); libera_matriz(u); libera_matriz(f); libera_matriz(r);
    return 0;
}
