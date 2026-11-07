#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "nucleo.h"
#include "calor.h"
#include "residuo_calor.h"
#include "suavizacao_calor.h"
#include "calor_interno.h"
#ifdef PLOT
void salva_plot_calor(Matriz u, int n, double t, const char *instante)
{
    char nome_dat[64];
    snprintf(nome_dat, sizeof(nome_dat),
             "grafico_%s.dat", instante);

    FILE *gnu = fopen(nome_dat, "w");

    if (gnu == NULL) {
        perror(nome_dat);
        exit(EXIT_FAILURE);
    }

    double h = 1.0/n;

#ifdef DIM2

    /*
      Formato:
     
      x y numerica exata erro
     */
    for (int j = 0; j <= n; j++) {
        for (int i = 0; i <= n; i++) {

            double x = i*h;
            double y = j*h;

            double numerica = u[i][j];
            double exata = uex_calor(x, y, t);

            fprintf(gnu,
                    "%.15e %.15e %.15e %.15e %.15e\n",
                    x,
                    y,
                    numerica,
                    exata,
                    fabs(numerica - exata));
        }

        fprintf(gnu, "\n");
    }

#elif defined(DIM3)

    /*
     Formato:
     
     x y z numerica exata erro
     */
    for (int k = 0; k <= n; k++) {
        for (int j = 0; j <= n; j++) {
            for (int i = 0; i <= n; i++) {

                double x = i*h;
                double y = j*h;
                double z = k*h;

                double numerica = u[i][j][k];
                double exata = uex_calor(x, y, z, t);

                fprintf(gnu,
                        "%.15e %.15e %.15e %.15e %.15e %.15e\n",
                        x,
                        y,
                        z,
                        numerica,
                        exata,
                        fabs(numerica - exata));
            }

            fprintf(gnu, "\n");
        }

        fprintf(gnu, "\n");
    }

#else

    /*
      Formato:
     
      x numerica exata erro
     */
    for (int i = 0; i <= n; i++) {

        double x = i*h;

        double numerica = u[i];
        double exata = uex_calor(x, t);

        fprintf(gnu,
                "%.15e %.15e %.15e %.15e\n",
                x,
                numerica,
                exata,
                fabs(numerica - exata));
    }

#endif

    fclose(gnu);


#ifdef DIM3

    /*
      Em 3D também gera um VTK para visualização
      completa no ParaView.
     */

    char nome_vtk[64];

    snprintf(nome_vtk, sizeof(nome_vtk),
             "grafico_%s.vtk", instante);

    FILE *vtk = fopen(nome_vtk, "w");

    if (vtk == NULL) {
        perror(nome_vtk);
        exit(EXIT_FAILURE);
    }

    /*
      Existem n+1 pontos:
     
      0, h, 2h, ..., 1
     */
    int N = n + 1;

    long long pontos =
        (long long)N * N * N;

    fprintf(vtk,
        "# vtk DataFile Version 3.0\n"
        "Multigrid 3D Calor\n"
        "ASCII\n"
        "DATASET STRUCTURED_POINTS\n"
        "DIMENSIONS %d %d %d\n"
        "ORIGIN 0 0 0\n"
        "SPACING %.15e %.15e %.15e\n"
        "POINT_DATA %lld\n"
        "SCALARS numerica double 1\n"
        "LOOKUP_TABLE default\n",
        N, N, N,
        h, h, h,
        pontos);


    /*
     NUMÉRICA
     
      x deve variar mais rápido no formato VTK.
     */
    for (int k = 0; k <= n; k++)
        for (int j = 0; j <= n; j++)
            for (int i = 0; i <= n; i++)
                fprintf(vtk,
                        "%.15e\n",
                        u[i][j][k]);


    /*
     * EXATA
     */
    fprintf(vtk,
            "SCALARS exata double 1\n"
            "LOOKUP_TABLE default\n");

    for (int k = 0; k <= n; k++)
        for (int j = 0; j <= n; j++)
            for (int i = 0; i <= n; i++)
                fprintf(vtk,
                        "%.15e\n",
                        uex_calor(i*h,
                                  j*h,
                                  k*h,
                                  t));


    /*
     * ERRO ABSOLUTO
     */
    fprintf(vtk,
            "SCALARS erro double 1\n"
            "LOOKUP_TABLE default\n");

    for (int k = 0; k <= n; k++)
        for (int j = 0; j <= n; j++)
            for (int i = 0; i <= n; i++) {

                double exata =
                    uex_calor(i*h,
                              j*h,
                              k*h,
                              t);

                fprintf(vtk,
                        "%.15e\n",
                        fabs(u[i][j][k] - exata));
            }

    fclose(vtk);

#endif
}
#endif
