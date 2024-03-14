#include "fem.h"
#include <math.h>

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
// double geoSize(double x, double y)
// {
//     femGeo *theGeometry = geoGetGeometry();
    
//     double h = theGeometry->h;

//     double x0 = theGeometry->xNotch;
//     double y0 = theGeometry->yNotch;
//     double r0 = theGeometry->rNotch;
//     double h0 = theGeometry->hNotch;
//     double d0 = theGeometry->dNotch;
    
//     double x1 = theGeometry->xHole;
//     double y1 = theGeometry->yHole;
//     double r1 = theGeometry->rHole;
//     double h1 = theGeometry->hHole;
//     double d1 = theGeometry->dHole;

//     // Strip : BEGIN
//     double d_EuclNotch = getEuclidianDistance(x, y, x0, y0);
//     double d_EuclHole  = getEuclidianDistance(x, y, x1, y1);

//     double h_Notch = hermiteInterpolation(d_EuclNotch - r0, h, h0, d0);
//     double h_Hole  = hermiteInterpolation(d_EuclHole - r1, h, h1, d1);

//     return fmin(h_Notch, h_Hole);
//     // Strip : END
// }


// #define ___ 0

/*
Generate the mesh of the plate.

@params:
    - None

@returns:
    - the mesh of the plate.
*/
// void geoMeshGenerate()
// {
//     femGeo *theGeometry = geoGetGeometry();

//     double w = theGeometry->LxPlate;
//     double h = theGeometry->LyPlate;
     
//     double x0 = theGeometry->xNotch;
//     double y0 = theGeometry->yNotch;
//     double r0 = theGeometry->rNotch;
    
//     double x1 = theGeometry->xHole;
//     double y1 = theGeometry->yHole;
//     double r1 = theGeometry->rHole;
 
//     //
//     //  -1- Construction de la geometrie avec OpenCascade
//     //      On cree le rectangle
//     //      On cree les deux cercles
//     //      On soustrait les cercles du rectangle :-)
//     //

//     // Strip : BEGIN
//     double xPlate = theGeometry->xPlate;
//     double yPlate = theGeometry->yPlate;
//     double zPlate = 0.0; // 2D
 
//     int ierr;

//     int idPlate = gmshModelOccAddRectangle(x0, y0, 0.0, w, h, -1, 0, &ierr);
//     // ErrorGmsh(ierr);

//     int idNotch  = gmshModelOccAddDisk(x0, y0, 0.0, r0, r0, -1, NULL, 0, NULL, 0, &ierr);    
//     // ErrorGmsh(ierr);

//     int idHole = gmshModelOccAddDisk(x1, y1, 0.0, r1, r1, -1, NULL, 0, NULL, 0, &ierr);
//     // ErrorGmsh(ierr);

//     // First parameter is the dimension : 2 because 2D
//     // Second parameter is the id of the object (the tag of the object in Gmsh's memory)
//     int plate[] = {2, idPlate};
//     int notch[] = {2, idNotch};
//     int hole[]  = {2, idHole};

//     gmshModelOccCut(plate, 2, notch, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
//     // ErrorGmsh(ierr);

//     gmshModelOccCut(plate, 2, hole, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
//     // ErrorGmsh(ierr);
//     // Strip : END

//     //
//     //  -2- Definition de la fonction callback pour la taille de reference
//     //      Synchronisation de OpenCascade avec gmsh
//     //      Generation du maillage (avec l'option Mesh.SaveAll :-)
   
//     geoSetSizeCallback(geoSize);
//     gmshModelOccSynchronize(&ierr);       
//     gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//     gmshModelMeshGenerate(2, &ierr);  

//     //
//     //  Generation de quads :-)
//     //
    
//     // gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//     // gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
//     // gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  chk(ierr);
//     // gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  chk(ierr);
//     // gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  chk(ierr);
//     // gmshModelMeshGenerate(2, &ierr);  

    
//     //
//     //  Plot of Fltk
//     //

//     // gmshFltkInitialize(&ierr);
//     // gmshFltkRun(&ierr);  chk(ierr);
// }


double geoSize(double x, double y)
{
    return 0.7;
}


void geoMeshGenerate()
{
    double rxLittleArcs = 5.0;
    double rxBigArcs = 2 * rxLittleArcs;
    double ryLittleArcs = 5.0;
    double ryBigArcs = 4.0;

    double widthPillier = 3.0;
    double offsetPont = 4.0;

    double Lx = 2 * rxBigArcs + 6 * rxLittleArcs + 4 * widthPillier;
    double Ly = ryLittleArcs + offsetPont;
 
    int ierr;

    // Ajout du plateau de base

    int idPlate = gmshModelOccAddRectangle(- Lx / 2, 0.0, 0.0, Lx, Ly, -1, 0, &ierr);

    // Ajout des arcs

    double xLittleArc1 = - Lx / 2;
    double yLittleArc1 = 0.0;

    double xLittleArc2 = - widthPillier - rxBigArcs - rxLittleArcs;
    double yLittleArc2 = 0.0;

    double xBigArc = 0.0;
    double yBigArc = 0.0;

    double xLittleArc3 = widthPillier + rxBigArcs + rxLittleArcs;
    double yLittleArc3 = 0.0;

    double xLittleArc4 = Lx / 2;
    double yLittleArc4 = 0.0;

    int idLittleArc1 = gmshModelOccAddDisk(xLittleArc1, yLittleArc1, 0.0, rxLittleArcs, ryLittleArcs, -1, NULL, 0, NULL, 0, &ierr);
    int idLittleArc2 = gmshModelOccAddDisk(xLittleArc2, yLittleArc2, 0.0, rxLittleArcs, ryLittleArcs, -1, NULL, 0, NULL, 0, &ierr);
    int idBigArc     = gmshModelOccAddDisk(xBigArc, yBigArc, 0.0, rxBigArcs, ryBigArcs, -1, NULL, 0, NULL, 0, &ierr);
    int idLittleArc3 = gmshModelOccAddDisk(xLittleArc3, yLittleArc3, 0.0, rxLittleArcs, ryLittleArcs, -1, NULL, 0, NULL, 0, &ierr);
    int idLittleArc4 = gmshModelOccAddDisk(xLittleArc4, yLittleArc4, 0.0, rxLittleArcs, ryLittleArcs, -1, NULL, 0, NULL, 0, &ierr);

    // Ajout des piliers

    double depthPillier = 5.0;

    double xPillier1 = - Lx / 2 + rxLittleArcs;
    double yPillier1 = - depthPillier;

    double xPillier2 = - rxBigArcs - widthPillier;
    double yPillier2 = - depthPillier;

    double xPillier3 = rxBigArcs;
    double yPillier3 = - depthPillier;

    double xPillier4 = Lx / 2 - rxLittleArcs - widthPillier;
    double yPillier4 = - depthPillier;

    int idPillier1 = gmshModelOccAddRectangle(xPillier1, yPillier1, 0.0, widthPillier, depthPillier, -1, 0, &ierr);
    int idPillier2 = gmshModelOccAddRectangle(xPillier2, yPillier2, 0.0, widthPillier, depthPillier, -1, 0, &ierr);
    int idPillier3 = gmshModelOccAddRectangle(xPillier3, yPillier3, 0.0, widthPillier, depthPillier, -1, 0, &ierr);
    int idPillier4 = gmshModelOccAddRectangle(xPillier4, yPillier4, 0.0, widthPillier, depthPillier, -1, 0, &ierr);

    // Ajout des colonnes

    double widthBigColumn = 2.5;
    double heightBigColumn = 6.0;

    double xBigColumn1 = - rxBigArcs - widthPillier - rxLittleArcs - widthBigColumn / 2;
    double yBigColumn1 = Ly;

    double xBigColumn2 = - rxBigArcs - widthPillier - rxLittleArcs - widthBigColumn / 4;
    double yBigColumn2 = Ly + heightBigColumn;

    double xBigColumn3 = - rxBigArcs - widthPillier - rxLittleArcs - widthBigColumn / 8;
    double yBigColumn3 = Ly + 5 * heightBigColumn / 3;

    double rBigDisk1 = 0.8;
    double xBigDisk1 = - rxBigArcs - widthPillier - rxLittleArcs;
    double yBigDisk1 = Ly + 2 * heightBigColumn + rBigDisk1 / 3;

    double xBigColumn4 = rxBigArcs + widthPillier + rxLittleArcs - widthBigColumn / 2;
    double yBigColumn4 = Ly;

    double xBigColumn5 = rxBigArcs + widthPillier + rxLittleArcs - widthBigColumn / 4;
    double yBigColumn5 = Ly + heightBigColumn;

    double xBigColumn6 = rxBigArcs + widthPillier + rxLittleArcs - widthBigColumn / 8;
    double yBigColumn6 = Ly + 5 * heightBigColumn / 3;

    double rBigDisk2 = 0.8;
    double xBigDisk2 = rxBigArcs + widthPillier + rxLittleArcs;
    double yBigDisk2 = Ly + 2 * heightBigColumn + rBigDisk2 / 3;

    int idBigColumn1 = gmshModelOccAddRectangle(xBigColumn1, yBigColumn1, 0.0, widthBigColumn, heightBigColumn, -1, 0, &ierr);
    int idBigColumn2 = gmshModelOccAddRectangle(xBigColumn2, yBigColumn2, 0.0, widthBigColumn / 2, heightBigColumn / 1.5, -1, 0, &ierr);
    int idBigColumn3 = gmshModelOccAddRectangle(xBigColumn3, yBigColumn3, 0.0, widthBigColumn / 4, heightBigColumn / 3, -1, 0, &ierr);
    int idBigDisk1   = gmshModelOccAddDisk(xBigDisk1, yBigDisk1, 0.0, rBigDisk1, rBigDisk1, -1, NULL, 0, NULL, 0, &ierr);
    int idBigColumn4 = gmshModelOccAddRectangle(xBigColumn4, yBigColumn4, 0.0, widthBigColumn, heightBigColumn, -1, 0, &ierr);
    int idBigColumn5 = gmshModelOccAddRectangle(xBigColumn5, yBigColumn5, 0.0, widthBigColumn / 2, heightBigColumn / 1.5, -1, 0, &ierr);
    int idBigColumn6 = gmshModelOccAddRectangle(xBigColumn6, yBigColumn6, 0.0, widthBigColumn / 4, heightBigColumn / 3, -1, 0, &ierr);
    int idBigDisk2   = gmshModelOccAddDisk(xBigDisk2, yBigDisk2, 0.0, rBigDisk2, rBigDisk2, -1, NULL, 0, NULL, 0, &ierr);

    // Ajout du second plateau

    double widthPlate2 = Lx;
    double heightPlate2 = 1.0;

    double xPlate2 = - Lx / 2;
    double yPlate2 = 0.0;

    int idPlate2 = gmshModelOccAddRectangle(xPlate2, yPlate2, 0.0, widthPlate2, heightPlate2, -1, 0, &ierr);

    // Ajout des colonnes du second plateau

    double widthLittleColumn = 0.5;
    double heightLittleColumn = 4;

    double xLittleColumn1 = - Lx / 2 + rxLittleArcs / 2 - widthLittleColumn / 2;
    double yLittleColumn1 = heightPlate2;

    double xLittleColumn2 = - rxBigArcs - widthPillier - 3 * rxLittleArcs / 2 - widthLittleColumn / 2;
    double yLittleColumn2 = heightPlate2;

    double xLittleColumn3 = - rxBigArcs - widthPillier - rxLittleArcs / 2 - widthLittleColumn / 2;
    double yLittleColumn3 = heightPlate2;

    double xLittleColumn4 = - rxBigArcs / 4 - widthLittleColumn / 2;
    double yLittleColumn4 = heightPlate2;

    double xLittleColumn5 = - 3 * rxBigArcs / 4 - widthLittleColumn / 2;
    double yLittleColumn5 = heightPlate2;

    double xLittleColumn6 =  rxBigArcs / 4 - widthLittleColumn / 2;
    double yLittleColumn6 = heightPlate2;

    double xLittleColumn7 = 3 * rxBigArcs / 4 - widthLittleColumn / 2;
    double yLittleColumn7 = heightPlate2;

    double xLittleColumn8 = rxBigArcs + widthPillier + rxLittleArcs / 2 - widthLittleColumn / 2;
    double yLittleColumn8 = heightPlate2;

    double xLittleColumn9 = rxBigArcs + widthPillier + 3 *rxLittleArcs / 2 - widthLittleColumn / 2;
    double yLittleColumn9 = heightPlate2;

    double xLittleColumn10 = Lx / 2 - rxLittleArcs / 2 - widthLittleColumn / 2;
    double yLittleColumn10 = heightPlate2;

    int idLittleColumn1  = gmshModelOccAddRectangle(xLittleColumn1, yLittleColumn1, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn2  = gmshModelOccAddRectangle(xLittleColumn2, yLittleColumn2, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn3  = gmshModelOccAddRectangle(xLittleColumn3, yLittleColumn3, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn4  = gmshModelOccAddRectangle(xLittleColumn4, yLittleColumn4, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn5  = gmshModelOccAddRectangle(xLittleColumn5, yLittleColumn5, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn6  = gmshModelOccAddRectangle(xLittleColumn6, yLittleColumn6, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn7  = gmshModelOccAddRectangle(xLittleColumn7, yLittleColumn7, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn8  = gmshModelOccAddRectangle(xLittleColumn8, yLittleColumn8, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn9  = gmshModelOccAddRectangle(xLittleColumn9, yLittleColumn9, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn10 = gmshModelOccAddRectangle(xLittleColumn10, yLittleColumn10, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);

    double widthWindow = widthPillier / 2;
    double heightWindow = 0.6;

    double xWindow1 = - rxBigArcs - widthPillier + widthPillier / 4;
    double yWindow1 = heightPlate2 + 0.5;

    double xWindow2 = - rxBigArcs - 2 * widthPillier - 2 * rxLittleArcs + widthPillier / 4;
    double yWindow2 = heightPlate2 + 0.5;

    double xWindow3 = rxBigArcs + widthPillier / 4;
    double yWindow3 = heightPlate2 + 0.5;

    double xWindow4 = rxBigArcs + widthPillier + 2 * rxLittleArcs + widthPillier / 4;
    double yWindow4 = heightPlate2 + 0.5;

    int idWindow1 = gmshModelOccAddRectangle(xWindow1, yWindow1, 0.0, widthWindow, heightWindow, -1, 0, &ierr);
    int idWindow2 = gmshModelOccAddRectangle(xWindow2, yWindow2, 0.0, widthWindow, heightWindow, -1, 0, &ierr);
    int idWindow3 = gmshModelOccAddRectangle(xWindow3, yWindow3, 0.0, widthWindow, heightWindow, -1, 0, &ierr);
    int idWindow4 = gmshModelOccAddRectangle(xWindow4, yWindow4, 0.0, widthWindow, heightWindow, -1, 0, &ierr);

    // Ajout des cables

    double widthCable = 0.2;
    double heightCable = 14.5;
    double distanceCable = 0.8;

    double xCable1 = xBigColumn2 + widthBigColumn / 2;
    double yCable1 = yBigColumn3;

    double xCable2 = xBigColumn2 + widthBigColumn / 2;
    double yCable2 = yCable1 - distanceCable;

    double xCable3 = xBigColumn2 + widthBigColumn / 2;
    double yCable3 = yCable2 - distanceCable;

    double xCable4 = xBigColumn2 + widthBigColumn / 2;
    double yCable4 = yCable3 - distanceCable;

    double xCable5 = xBigColumn1 + widthBigColumn;
    double yCable5 = yCable4 - 2 * distanceCable;

    double xCable6 = xBigColumn1 + widthBigColumn;
    double yCable6 = yCable5 - distanceCable;

    double xCable7 = xBigColumn1 + widthBigColumn;
    double yCable7 = yCable6 - distanceCable;

    double xCable8 = xBigColumn1 + widthBigColumn;
    double yCable8 = yCable7 - distanceCable;

    double xCable9 = xBigColumn1 + widthBigColumn;
    double yCable9 = yCable8 - distanceCable;

    int idCable1 = gmshModelOccAddRectangle(xCable1, yCable1, 0.0, widthCable, heightCable, -1, 0, &ierr);
    int idCable2 = gmshModelOccAddRectangle(xCable2, yCable2, 0.0, widthCable, heightCable - 1, -1, 0, &ierr);
    int idCable3 = gmshModelOccAddRectangle(xCable3, yCable3, 0.0, widthCable, heightCable - 2, -1, 0, &ierr);
    int idCable4 = gmshModelOccAddRectangle(xCable4, yCable4, 0.0, widthCable, heightCable - 3, -1, 0, &ierr);
    int idCable5 = gmshModelOccAddRectangle(xCable5, yCable5, 0.0, widthCable, heightCable - 5, -1, 0, &ierr);
    int idCable6 = gmshModelOccAddRectangle(xCable6, yCable6, 0.0, widthCable, heightCable - 6, -1, 0, &ierr);
    int idCable7 = gmshModelOccAddRectangle(xCable7, yCable7, 0.0, widthCable, heightCable - 7, -1, 0, &ierr);
    int idCable8 = gmshModelOccAddRectangle(xCable8, yCable8, 0.0, widthCable, heightCable - 8, -1, 0, &ierr);
    int idCable9 = gmshModelOccAddRectangle(xCable9, yCable9, 0.0, widthCable, heightCable - 9, -1, 0, &ierr);

    heightCable = 14.5;

    double xCable10 = xBigColumn5 + widthBigColumn / 2;
    double yCable10 = yBigColumn6;

    double xCable11 = xBigColumn5 + widthBigColumn / 2;
    double yCable11 = yCable10 - distanceCable;

    double xCable12 = xBigColumn5 + widthBigColumn / 2;
    double yCable12 = yCable11 - distanceCable;

    double xCable13 = xBigColumn5 + widthBigColumn / 2;
    double yCable13 = yCable12 - distanceCable;

    double xCable14 = xBigColumn4 + widthBigColumn;
    double yCable14 = yCable13 - 2 * distanceCable;

    double xCable15 = xBigColumn4 + widthBigColumn;
    double yCable15 = yCable14 - distanceCable;

    double xCable16 = xBigColumn4 + widthBigColumn;
    double yCable16 = yCable15 - distanceCable;

    double xCable17 = xBigColumn4 + widthBigColumn;
    double yCable17 = yCable16 - distanceCable;

    double xCable18 = xBigColumn4 + widthBigColumn;
    double yCable18 = yCable17 - distanceCable;

    int idCable10 = gmshModelOccAddRectangle(xCable10, yCable10, 0.0, widthCable, heightCable, -1, 0, &ierr);
    int idCable11 = gmshModelOccAddRectangle(xCable11, yCable11, 0.0, widthCable, heightCable - 1, -1, 0, &ierr);
    int idCable12 = gmshModelOccAddRectangle(xCable12, yCable12, 0.0, widthCable, heightCable - 2, -1, 0, &ierr);
    int idCable13 = gmshModelOccAddRectangle(xCable13, yCable13, 0.0, widthCable, heightCable - 3, -1, 0, &ierr);
    int idCable14 = gmshModelOccAddRectangle(xCable14, yCable14, 0.0, widthCable, heightCable - 5, -1, 0, &ierr);
    int idCable15 = gmshModelOccAddRectangle(xCable15, yCable15, 0.0, widthCable, heightCable - 6, -1, 0, &ierr);
    int idCable16 = gmshModelOccAddRectangle(xCable16, yCable16, 0.0, widthCable, heightCable - 7, -1, 0, &ierr);
    int idCable17 = gmshModelOccAddRectangle(xCable17, yCable17, 0.0, widthCable, heightCable - 8, -1, 0, &ierr);
    int idCable18 = gmshModelOccAddRectangle(xCable18, yCable18, 0.0, widthCable, heightCable - 9, -1, 0, &ierr);


    int plate[]       = {2, idPlate};
    int plate2[]      = {2, idPlate2};

    int littleArc1[]  = {2, idLittleArc1};
    int littleArc2[]  = {2, idLittleArc2};
    int bigArc[]      = {2, idBigArc};
    int littleArc3[]  = {2, idLittleArc3};
    int littleArc4[]  = {2, idLittleArc4};

    int pillier1[]    = {2, idPillier1};
    int pillier2[]    = {2, idPillier2};
    int pillier3[]    = {2, idPillier3};
    int pillier4[]    = {2, idPillier4};

    int littleColumn1[]  = {2, idLittleColumn1};
    int littleColumn2[]  = {2, idLittleColumn2};
    int littleColumn3[]  = {2, idLittleColumn3};
    int littleColumn4[]  = {2, idLittleColumn4};
    int littleColumn5[]  = {2, idLittleColumn5};
    int littleColumn6[]  = {2, idLittleColumn6};
    int littleColumn7[]  = {2, idLittleColumn7};
    int littleColumn8[]  = {2, idLittleColumn8};
    int littleColumn9[]  = {2, idLittleColumn9};
    int littleColumn10[] = {2, idLittleColumn10};

    int bigColumn1[]  = {2, idBigColumn1};
    int bigColumn2[]  = {2, idBigColumn2};
    int bigColumn3[]  = {2, idBigColumn3};
    int bigColumn4[]  = {2, idBigColumn4};
    int bigColumn5[]  = {2, idBigColumn5};
    int bigColumn6[]  = {2, idBigColumn6};

    int bigDisk1[]  = {2, idBigDisk1};
    int bigDisk2[]  = {2, idBigDisk2};

    int window1[]  = {2, idWindow1};
    int window2[]  = {2, idWindow2};
    int window3[]  = {2, idWindow3};
    int window4[]  = {2, idWindow4};

    int cable1[]  = {2, idCable1};
    int cable2[]  = {2, idCable2};
    int cable3[]  = {2, idCable3};
    int cable4[]  = {2, idCable4};
    int cable5[]  = {2, idCable5};
    int cable6[]  = {2, idCable6};
    int cable7[]  = {2, idCable7};
    int cable8[]  = {2, idCable8};
    int cable9[]  = {2, idCable9};

    int cable10[]  = {2, idCable10};
    int cable11[]  = {2, idCable11};
    int cable12[]  = {2, idCable12};
    int cable13[]  = {2, idCable13};
    int cable14[]  = {2, idCable14};
    int cable15[]  = {2, idCable15};
    int cable16[]  = {2, idCable16};
    int cable17[]  = {2, idCable17};
    int cable18[]  = {2, idCable18};

    const double PI = 3.14159265358979323846;
    double angle_cable = -135 * PI / 180;
    
    gmshModelOccRotate(cable1, 2, xCable1, yCable1, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable2, 2, xCable2, yCable2, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable3, 2, xCable3, yCable3, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable4, 2, xCable4, yCable4, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable5, 2, xCable5, yCable5, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable6, 2, xCable6, yCable6, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable7, 2, xCable7, yCable7, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable8, 2, xCable8, yCable8, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable9, 2, xCable9, yCable9, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);

    gmshModelOccRotate(cable10, 2, xCable10, yCable10, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable11, 2, xCable11, yCable11, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable12, 2, xCable12, yCable12, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable13, 2, xCable13, yCable13, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable14, 2, xCable14, yCable14, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable15, 2, xCable15, yCable15, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable16, 2, xCable16, yCable16, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable17, 2, xCable17, yCable17, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    gmshModelOccRotate(cable18, 2, xCable18, yCable18, 0.0, 0.0, 0.0, 1.0, angle_cable, &ierr);
    
    // On soustrait les arcs du plateau de base
    gmshModelOccCut(plate, 2, littleArc1, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, littleArc2, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, bigArc, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, littleArc3, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, littleArc4, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    // Fusionner la structure de base ensemble
    gmshModelOccFuse(plate, 2, plate2, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    gmshModelOccFuse(plate, 2, pillier1, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(plate, 2, pillier2, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(plate, 2, pillier3, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(plate, 2, pillier4, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    gmshModelOccFuse(plate, 2, littleColumn1, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(plate, 2, littleColumn2, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(plate, 2, littleColumn3, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(plate, 2, littleColumn4, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(plate, 2, littleColumn5, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(plate, 2, littleColumn6, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(plate, 2, littleColumn7, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(plate, 2, littleColumn8, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(plate, 2, littleColumn9, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(plate, 2, littleColumn10, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    gmshModelOccFuse(bigColumn1, 2, bigColumn2, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn1, 2, bigColumn3, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn1, 2, bigDisk1, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    gmshModelOccFuse(bigColumn4, 2, bigColumn5, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn4, 2, bigColumn6, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn4, 2, bigDisk2, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    gmshModelOccCut(plate, 2, window1, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, window2, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, window3, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, window4, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    gmshModelOccFuse(bigColumn1, 2, cable1, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn1, 2, cable2, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn1, 2, cable3, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn1, 2, cable4, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn1, 2, cable5, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn1, 2, cable6, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn1, 2, cable7, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn1, 2, cable8, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn1, 2, cable9, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    gmshModelOccFuse(bigColumn4, 2, cable10, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn4, 2, cable11, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn4, 2, cable12, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn4, 2, cable13, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn4, 2, cable14, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn4, 2, cable15, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn4, 2, cable16, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn4, 2, cable17, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn4, 2, cable18, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    gmshModelOccFuse(plate, 2, bigColumn1, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(plate, 2, bigColumn4, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    geoSetSizeCallback(geoSize);
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);

    //  Plot of Fltk
    gmshFltkInitialize(&ierr);
}