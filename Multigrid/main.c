#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "nucleo.h"
#include "interface.h"
int main(int argc, char **argv)
{
    if (argc < 2) {
        uso(argv[0]);
        return 0;
    }

    int modo = 0;
    int submodo = 0;

    // strcmp exige que o argumento seja exatamente v
    if (argc > 1 && strcmp(argv[1], "v") == 0)
        modo = 1;

    //* submodo so existe depois de v: "v 1". 
    if (modo && argc > 2 && strcmp(argv[2], "1") == 0)
        submodo = 1;

    int base = 1 + modo + submodo;
    if (argc <= base) {
        uso(argv[0]);
        return 1;
    }

    return executa_argumentos(argc, argv, base, modo, submodo);
}
