#include "../include/fem_other.h"

double femMin(double *x, int n)
{
    double myMin = x[0];
    for (int i = 1; i < n; i++) { myMin = fmin(myMin, x[i]); }
    return myMin;
}

double femMax(double *x, int n)
{
    double myMax = x[0];
    for (int i = 1; i < n; i++) { myMax = fmax(myMax, x[i]); }
    return myMax;
}

void femError(char *text, int line, char *file)
{
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s:%d at line %d : \n  %s\n", file, line, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);
}

void femErrorScan(int test, int line, char *file)
{
    if (test >= 0) { return; }

    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in fscanf or fgets in %s:%d at line %d : \n", file, line, line);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);
}

void femWarning(char *text, int line, char *file)
{
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s:%d at line %d : \n  %s\n", file, line, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
}