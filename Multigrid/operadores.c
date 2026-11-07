
#include "matriz.h"
#include "operadores.h"

/*
Full weighting

passamos da malha fina para a malha grossa.

   r_H[I][J] =
   1/16 * (
       4*r_h[2I][2J]
       + 2*(vizinhos
       + (cantos)
   )

*/

void restringe(Matriz r_fina, Matriz r_grossa, int n_fina)
{
    int n_grossa = n_fina / 2;
    zera_matriz(r_grossa, n_grossa);

#ifdef DIM2
    for (int I = 1; I < n_grossa; I++) {
        for (int J = 1; J < n_grossa; J++) {
            int i = 2 * I, j = 2 * J;
            r_grossa[I][J] = (
                4.0 * r_fina[i][j]
                + 2.0 * (r_fina[i-1][j] + r_fina[i+1][j]
                       + r_fina[i][j-1] + r_fina[i][j+1])
                + (r_fina[i-1][j-1] + r_fina[i-1][j+1]
                 + r_fina[i+1][j-1] + r_fina[i+1][j+1])
            ) / 16.0;
        }
    }
#elif defined(DIM3)
    for (int I = 1; I < n_grossa; I++) {
        for (int J = 1; J < n_grossa; J++) {
            for (int K = 1; K < n_grossa; K++) {
                int i = 2 * I, j = 2 * J, k = 2 * K;
                double soma = 0.0;
                for (int di = -1; di <= 1; di++) {
                    int wi = (di == 0) ? 2 : 1;
                    for (int dj = -1; dj <= 1; dj++) {
                        int wj = (dj == 0) ? 2 : 1;
                        for (int dk = -1; dk <= 1; dk++) {
                            int wk = (dk == 0) ? 2 : 1;
                            soma += (double)(wi * wj * wk)
                                  * r_fina[i+di][j+dj][k+dk];
                        }
                    }
                }
                r_grossa[I][J][K] = soma / 64.0;
            }
        }
    }
#else 
    for (int I=1; I < n_grossa; I++) {
	    int i = 2*I;
	    r_grossa[I] = (
	        2.0*r_fina[i] 
		+
		r_fina[i-1]
		+r_fina[i+1]
		)/4.0;
    }
#endif
}




/*
correcao da malha grossa e interpolada
para a malha fina.

pontos comum recebem diretamente
o valor da malha grossa.

Os outros recebem medias
*/

void prolonga(Matriz e_grossa, Matriz e_fina, int n_fina)
{
    int n_grossa = n_fina / 2;
    zera_matriz(e_fina, n_fina);

#ifdef DIM2
    for (int i = 0; i <= n_fina; i++) {
        int I = i / 2;
        int px = i % 2;
        for (int j = 0; j <= n_fina; j++) {
            int J = j / 2;
            int py = j % 2;
            double soma = 0.0;
            int nx = px ? 2 : 1;
            int ny = py ? 2 : 1;
            for (int a = 0; a < nx; a++)
                for (int b = 0; b < ny; b++)
                    soma += e_grossa[I+a][J+b];
            e_fina[i][j] = soma / (double)(nx * ny);
        }
    }
#elif defined(DIM3)
    for (int i = 0; i <= n_fina; i++) {
        int I = i / 2;
        int px = i % 2;
        for (int j = 0; j <= n_fina; j++) {
            int J = j / 2;
            int py = j % 2;
            for (int k = 0; k <= n_fina; k++) {
                int K = k / 2;
                int pz = k % 2;
                int nx = px ? 2 : 1;
                int ny = py ? 2 : 1;
                int nz = pz ? 2 : 1;
                double soma = 0.0;
                for (int a = 0; a < nx; a++)
                    for (int b = 0; b < ny; b++)
                        for (int c = 0; c < nz; c++)
                            soma += e_grossa[I+a][J+b][K+c];
                e_fina[i][j][k] = soma / (double)(nx * ny * nz);
            }
        }
    }
#else 
    for (int i=0; i<=n_fina; i++) {
	int I = i/2;
	int px = i%2;
	int nx = px ? 2:1;
	double soma = 0.0;
	for (int a=0; a<nx; a++)
		soma += e_grossa[I+a];
	e_fina[i] = soma/(double)nx;
    }
#endif
    (void)n_grossa;
}
