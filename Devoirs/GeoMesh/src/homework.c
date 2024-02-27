#include "fem.h"

// Strip : BEGIN
/*
This function is a 1D hermite interpolation function.
It is used to interpolate the value of a function at a given point, given the value of the
function at two other points and the distance at which the function reaches its maximum value.

More explicitely, the Hermite interpolation is mathematically defined as follows:

    Conditions:
        1: h(0) = h_min
        2: h(dist_interp) = h_max
        3: h'(0) = 0
        4: h'(dist_interp) = 0

    Resolution:
        double a = - 2 * (h_max - h_min) / pow(dist_interp, 3);
        double b =   3 * (h_max - h_min) / pow(dist_interp, 2);
        double c = 0.0;
        double d = h_min;

        double h = a * pow(x, 3) + b * pow(x, 2) + c * x + d;

@params:
    - x:           the distance at which we want to interpolate the value of the function.
    - h_max:       the maximum value that the function can reach.
    - h_min:       the minimum value that the function can reach.
    - dist_interp: the distance at which the function reaches its maximum value.
                   Between [0, dist_interp], the function is interpolated using the hermite interpolation, defined above.

@returns:
    - the value of the hermite interpolation at the point x.
*/
double hermiteInterpolation(double x, double h_max, double h_min, double dist_interp)
{
    if (x >= dist_interp) return h_max;
    if (x <= 0)           return h_min;

    return (1 / (dist_interp * dist_interp)) * (h_max - h_min) * x * x * (- 2 * x / dist_interp + 3) + h_min;
}

/*
This function is used to calculate the euclidian distance between two points (x, y) and (x_ref, y_ref).

@params:
    - x1:  the x coordinate of the first point.
    - y1:  the y coordinate of the first point.
    - x2:  the x coordinate of the second point.
    - y2:  the y coordinate of the second point.

@returns:
    - the euclidian distance between (x1, y1) and (x2, y2).
*/
double getEuclidianDistance(double x1, double y1, double x2, double y2) { return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)); }
// Strip : END

/*
This function is used to get the geometry size of the plate at a given point.

@params:
    - x: the x coordinate of the point at which we want to get the geometry size.
    - y: the y coordinate of the point at which we want to get the geometry size.

@returns:
    - the geometry size of the plate at the point (x, y).
*/
double geoSize(double x, double y)
{
    femGeo *theGeometry = geoGetGeometry();
    
    double h = theGeometry->h;

    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    double h0 = theGeometry->hNotch;
    double d0 = theGeometry->dNotch;
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
    double h1 = theGeometry->hHole;
    double d1 = theGeometry->dHole;

    // Strip : BEGIN
    double d_EuclNotch = getEuclidianDistance(x, y, x0, y0);
    double d_EuclHole  = getEuclidianDistance(x, y, x1, y1);

    double h_Notch = hermiteInterpolation(d_EuclNotch - r0, h, h0, d0);
    double h_Hole  = hermiteInterpolation(d_EuclHole - r1, h, h1, d1);

    return fmin(h_Notch, h_Hole);
    // Strip : END
}


// #define ___ 0

/*
Generate the mesh of the plate.

@params:
    - None

@returns:
    - the mesh of the plate.
*/
void geoMeshGenerate()
{
    femGeo *theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
     
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
 
    //
    //  -1- Construction de la geometrie avec OpenCascade
    //      On cree le rectangle
    //      On cree les deux cercles
    //      On soustrait les cercles du rectangle :-)
    //

    // Strip : BEGIN
    double xPlate = theGeometry->xPlate;
    double yPlate = theGeometry->yPlate;
    double zPlate = 0.0; // 2D
 
    int ierr;

    int idPlate = gmshModelOccAddRectangle(x0, y0, 0.0, w, h, -1, 0, &ierr);
    ErrorGmsh(ierr);

    int idNotch  = gmshModelOccAddDisk(x0, y0, 0.0, r0, r0, -1, NULL, 0, NULL, 0, &ierr);    
    ErrorGmsh(ierr);

    int idHole = gmshModelOccAddDisk(x1, y1, 0.0, r1, r1, -1, NULL, 0, NULL, 0, &ierr);
    ErrorGmsh(ierr);

    // First parameter is the dimension : 2 because 2D
    // Second parameter is the id of the object (the tag of the object in Gmsh's memory)
    int plate[] = {2, idPlate};
    int notch[] = {2, idNotch};
    int hole[]  = {2, idHole};

    gmshModelOccCut(plate, 2, notch, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);

    gmshModelOccCut(plate, 2, hole, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);
    // Strip : END

    //
    //  -2- Definition de la fonction callback pour la taille de reference
    //      Synchronisation de OpenCascade avec gmsh
    //      Generation du maillage (avec l'option Mesh.SaveAll :-)
   
    geoSetSizeCallback(geoSize);
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  
       
    //
    //  Generation de quads :-)
    //
    
    // gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    // gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    // gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  chk(ierr);
    // gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  chk(ierr);
    // gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  chk(ierr);
    // gmshModelMeshGenerate(2, &ierr);  

    
    //
    //  Plot of Fltk
    //

    // gmshFltkInitialize(&ierr);
    // gmshFltkRun(&ierr);  chk(ierr);
}
