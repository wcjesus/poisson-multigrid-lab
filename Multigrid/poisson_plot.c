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
#ifdef PLOT
void salva_plot_poisson(Matriz u, int n)
{
    FILE *gnu = fopen("grafico.dat", "w");
#ifdef DIM2 
    for (int inteiro = 0; inteiro<n; inteiro++){
    	for (int jota =0; jota<n; jota++){
		double x = (double)inteiro/n;
		double y = (double)jota/n;
		fprintf(gnu,"%lf %lf %lf %lf\n",x,y,u[inteiro][jota], fabs(u[inteiro][jota]-uex(x,y)));
		
	}
	fprintf(gnu,"\n");
}fclose(gnu);
#elif defined(DIM1)
    
    for (int inteiro =0; inteiro<n; inteiro++){
		double x = (double)inteiro/n;
    		fprintf(gnu,"%lf %lf %lf\n", x, u[inteiro], fabs(u[inteiro]-uex(x)));}fprintf(gnu,"\n");fclose(gnu);		
#elif defined(DIM3)

    int passo = 1;
    int N = n / passo;

    FILE *vtk = fopen("grafico.vtk", "w");

    if (vtk == NULL) {
        perror("Erro ao abrir grafico.vtk");
        exit(1);
    }

    /*
      Arquivo VTK para o ParaView.
     
      A malha possui N x N x N pontos:
     
      0, 4/n, 8/n, ..., (n-4)/n
     
      Para n = 128:
      N = 32
     */

    fprintf(vtk,
        "# vtk DataFile Version 3.0\n"
        "Multigrid 3D\n"
        "ASCII\n"
        "DATASET STRUCTURED_POINTS\n"
        "DIMENSIONS %d %d %d\n"
        "ORIGIN 0 0 0\n"
        "SPACING %.15e %.15e %.15e\n"
        "\n"
        "POINT_DATA %d\n"
        "SCALARS numerica double 1\n"
        "LOOKUP_TABLE default\n",
        N, N, N,
        (double)passo / n,
        (double)passo / n,
        (double)passo / n,
        N * N * N);


    //grafico.dat continua sendo usado pelo gnuplot.
     
    

    if (gnu == NULL) {
        perror("Erro ao abrir grafico.dat");
        fclose(vtk);
        exit(1);
    }


    //NUMÉRICA
     
    for (int inteiro = 0; inteiro < n; inteiro += passo)
    {
        for (int jota = 0; jota < n; jota += passo)
        {
            for (int ka = 0; ka < n; ka += passo)
            {
                double x = (double)inteiro / n;
                double y = (double)jota / n;
                double z = (double)ka / n;

                double exata = uex(x, y, z);

                fprintf(gnu,
                    "%lf %lf %lf %lf %lf %lf\n",
                    x, y, z,
                    u[inteiro][jota][ka],
                    fabs(u[inteiro][jota][ka] - exata),
                    exata);

                fprintf(vtk,
                    "%.15e\n",
                    u[inteiro][jota][ka]);
            }

            fprintf(gnu, "\n");
        }

        fprintf(gnu, "\n");
    }


    //SOLUÇÃO EXATA
     
    fprintf(vtk,
        "\n"
        "SCALARS exata double 1\n"
        "LOOKUP_TABLE default\n");

    for (int inteiro = 0; inteiro < n; inteiro += passo)
    {
        for (int jota = 0; jota < n; jota += passo)
        {
            for (int ka = 0; ka < n; ka += passo)
            {
                double x = (double)inteiro / n;
                double y = (double)jota / n;
                double z = (double)ka / n;

                fprintf(vtk,
                    "%.15e\n",
                    uex(x, y, z));
            }
        }
    }


    //ERRO ABSOLUTO
     
    fprintf(vtk,
        "\n"
        "SCALARS erro double 1\n"
        "LOOKUP_TABLE default\n");

    for (int inteiro = 0; inteiro < n; inteiro += passo)
    {
        for (int jota = 0; jota < n; jota += passo)
        {
            for (int ka = 0; ka < n; ka += passo)
            {
                double x = (double)inteiro / n;
                double y = (double)jota / n;
                double z = (double)ka / n;

                double exata = uex(x, y, z);

                fprintf(vtk,
                    "%.15e\n",
                    fabs(u[inteiro][jota][ka] - exata));
            }
        }
    }

    fclose(gnu);
    fclose(vtk);
#endif
}
#endif
