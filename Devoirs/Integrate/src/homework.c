#include <math.h>
#include "glfem.h"

// Strip BEGIN
double interpolate(double v[3], double xi, double eta) { return v[0] * (1 - xi - eta) + v[1] * eta + v[2] * xi; }

double getJacobien(double x[3], double y[3]) { return (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]); }
// Strip END

double integrate(double x[3], double y[3], double (*f) (double, double))
{
    // Strip BEGIN
    // Const values for the integration points (see guidelines)
    const double value = 1.0 / 6.0;
    const double xi[3]  = { value, 4 * value, value };
    const double eta[3] = { value, value, 4 * value };
    const double weight      = value; // All the weights are the same

    const double jacobien = getJacobien(x, y);

    double I = 0.0;
    double xLoc[3];
    double yLoc[3];

    // Calculate the approximate integral
    for (int i = 0; i < 3; i++)
    {
        xLoc[i] = interpolate(x, xi[i], eta[i]);
        yLoc[i] = interpolate(y, xi[i], eta[i]);
        I += f(xLoc[i], yLoc[i]);
    }

    // Multiply by the jacobien for the change of variable
    I *= weight * fabs(jacobien);
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
    
    // Calculate the new points (= middle of the edges of the triangle)
    const double a_x = (x[1] + x[2]) / 2;
    const double a_y = (y[1] + y[2]) / 2;
    const double b_x = (x[2] + x[0]) / 2;
    const double b_y = (y[2] + y[0]) / 2;
    const double c_x = (x[1] + x[0]) / 2;
    const double c_y = (y[1] + y[0]) / 2;
    
    // Calculate new coordinates for the 4 new triangles
    double newPtsX[4][3] = {{x[0], c_x, b_x}, {c_x, x[1], a_x}, {b_x, a_x, x[2]}, {b_x, a_x, c_x}};
    double newPtsY[4][3] = {{y[0], c_y, b_y}, {c_y, y[1], a_y}, {b_y, a_y, y[2]}, {b_y, a_y, c_y}};

    double I = 0.0;

    // Calculate the integral for each of the 4 new triangles recursively
    for (int i = 0; i < 4; i++) { I += integrateRecursive(newPtsX[i], newPtsY[i], f, n - 1); }
    // Strip END

    return I;
}
