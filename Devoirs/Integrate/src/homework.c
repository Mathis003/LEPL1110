#include <math.h>
#include "glfem.h"

// Strip BEGIN
double interpolate(double u[3], double xi, double eta) { return u[0] * xi + u[1] * eta + u[2] * (1 - xi - eta); }

double getJacobien(double x[3], double y[3]) { return (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]); }
// Strip END

double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0.0;
    double xLoc[3];
    double yLoc[3];

    // Strip BEGIN
    const double value = 1.0 / 6.0;
    const double xi[3]     = { value, 4 * value, value };
    const double eta[3]    = { value, value, 4 * value };
    const double weight[3] = { value, value, value };

    double jacobien = getJacobien(x, y);

    for (int i = 0; i < 3; i++)
    {
        xLoc[i] = interpolate(x, xi[i], eta[i]);
        yLoc[i] = interpolate(y, xi[i], eta[i]);
        I += weight[i] * f(xLoc[i], yLoc[i]);
    }

    I *= fabs(jacobien);
    // Strip END

    // Pour dessiner l'element, les sommets du triangle :-)
    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x, y, 3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x, y, 3);

    // Decommenter la ligne pour dessiner aussi les points d'integration
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc, yLoc, 3);

    return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{
    // Strip BEGIN
    if (n <= 0) { return integrate(x, y, f); }
    
    // Calcul des nouveaux points (milieux des côtés du triangle)
    double a[2] = { fabs((x[1] + x[2]) / 2), fabs((y[1] + y[2]) / 2) };
    double b[2] = { fabs((x[2] + x[0]) / 2), fabs((y[2] + y[0]) / 2) };
    double c[2] = { fabs((x[1] + x[0]) / 2), fabs((y[1] + y[0]) / 2) };
    
    // Calcul des coordonnées des quatres triangle
    double newPtsX[4][3] = {{x[0], c[0], b[0]}, {c[0], x[1], a[0]}, {b[0], a[0], x[2]}, {b[0], a[0], c[0]}};
    double newPtsY[4][3] = {{y[0], c[1], b[1]}, {c[1], y[1], a[1]}, {b[1], a[1], y[2]}, {b[1], a[1], c[1]}};

    double I = 0.0;
    for (int i = 0; i < 4; i++) { I += integrateRecursive(newPtsX[i], newPtsY[i], f, n - 1); }
    // Strip END

    return I;
}
