#include <stdio.h>
#include <math.h>
#include "glfem.h"


double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];

    // TODO: BEGIN

    // Pour dessiner l'element, les sommets du triangle :-)
    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);

    // Decommenter la ligne pour dessiner aussi les points d'integration
    // glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);

    // TODO: END

    return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{
    // TODO: BEGIN

    // Modifier cette ligne aussi !
    double I = integrate(x,y,f);
    
    // TODO: END

    return I;
}