#include <stdio.h>
#include <math.h>

#include "fem.h"
#include "glfem.h"



double integrate(double x[3],double y[3],double(*f)(double,double));
double integrateRecursive(double x[3],double y[3],double(*f)(double,double),int n);

double fun(double x,double y)      { return cos(x) + y * y; }
double stupid(double x,double y)   { return 1.0; }


int main(int argc, char* argv[])
{

    char   theMessage[256];
   
    double x[3] = { 0, 1, 0};
    double y[3] = { 0, 0, 1};
 
    glfemWindowCreate("EPL1110 : Integrate",480,480,3,x,y);
    do {
        glfemReshape(x,y,3);

        double I = integrateRecursive(x,y,fun,2);
        sprintf(theMessage, "Integral = %14.7e",I); 
        glfemDrawMessage(theMessage,(double[2]){16.0, 30.0});

        glfemWindowUpdate();
    } while(!glfemWindowShouldClose());
    
    glfemWindowFree();
    exit(EXIT_SUCCESS);
    return 0;
}


