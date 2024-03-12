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
    return 0.8;
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

    int idPlate = gmshModelOccAddRectangle(- Lx / 2, - Ly, 0.0, Lx, Ly, -1, 0, &ierr);

    // Ajout des arcs

    double xLittleArc1 = - Lx / 2;
    double yLittleArc1 = - Ly;

    double xLittleArc2 = - widthPillier - rxBigArcs - rxLittleArcs;
    double yLittleArc2 = - Ly;

    double xBigArc = 0;
    double yBigArc = - Ly;

    double xLittleArc3 = widthPillier + rxBigArcs + rxLittleArcs;
    double yLittleArc3 = - Ly;

    double xLittleArc4 = Lx / 2;
    double yLittleArc4 = - Ly;

    int idLittleArc1 = gmshModelOccAddDisk(xLittleArc1, yLittleArc1, 0.0, rxLittleArcs, ryLittleArcs, -1, NULL, 0, NULL, 0, &ierr);
    int idLittleArc2 = gmshModelOccAddDisk(xLittleArc2, yLittleArc2, 0.0, rxLittleArcs, ryLittleArcs, -1, NULL, 0, NULL, 0, &ierr);
    int idBigArc = gmshModelOccAddDisk(xBigArc, yBigArc, 0.0, rxBigArcs, ryBigArcs, -1, NULL, 0, NULL, 0, &ierr);
    int idLittleArc3 = gmshModelOccAddDisk(xLittleArc3, yLittleArc3, 0.0, rxLittleArcs, ryLittleArcs, -1, NULL, 0, NULL, 0, &ierr);
    int idLittleArc4 = gmshModelOccAddDisk(xLittleArc4, yLittleArc4, 0.0, rxLittleArcs, ryLittleArcs, -1, NULL, 0, NULL, 0, &ierr);

    // Ajout des piliers

    double depthPillier = 5.0;

    double xPillier1 = - Lx / 2 + rxLittleArcs;
    double yPillier1 = - Ly - depthPillier;

    double xPillier2 = - rxBigArcs - widthPillier;
    double yPillier2 = - Ly - depthPillier;

    double xPillier3 = rxBigArcs;
    double yPillier3 = - Ly - depthPillier;

    double xPillier4 = Lx / 2 - rxLittleArcs - widthPillier;
    double yPillier4 = - Ly - depthPillier;

    int idPillier1 = gmshModelOccAddRectangle(xPillier1, yPillier1, 0.0, widthPillier, depthPillier, -1, 0, &ierr);
    int idPillier2 = gmshModelOccAddRectangle(xPillier2, yPillier2, 0.0, widthPillier, depthPillier, -1, 0, &ierr);
    int idPillier3 = gmshModelOccAddRectangle(xPillier3, yPillier3, 0.0, widthPillier, depthPillier, -1, 0, &ierr);
    int idPillier4 = gmshModelOccAddRectangle(xPillier4, yPillier4, 0.0, widthPillier, depthPillier, -1, 0, &ierr);

    // Ajout des colonnes

    double widthBigColumn = 1.5;
    double heightBigColumn = 15.0;

    double xBigColumn1 = - rxBigArcs - widthPillier - rxLittleArcs - widthBigColumn / 2;
    double yBigColumn1 = 0.0;

    double xBigColumn2 = rxBigArcs + widthPillier + rxLittleArcs - widthBigColumn / 2;
    double yBigColumn2 = 0.0;

    int idBigColumn1 = gmshModelOccAddRectangle(xBigColumn1, yBigColumn1, 0.0, widthBigColumn, heightBigColumn, -1, 0, &ierr);
    int idBigColumn2 = gmshModelOccAddRectangle(xBigColumn2, yBigColumn2, 0.0, widthBigColumn, heightBigColumn, -1, 0, &ierr);

    // Ajout du second plateau

    double widthPlate2 = Lx;
    double heightPlate2 = 1.0;

    double xPlate2 = - Lx / 2;
    double yPlate2 = - Ly;

    int idPlate2 = gmshModelOccAddRectangle(xPlate2, yPlate2, 0.0, widthPlate2, heightPlate2, -1, 0, &ierr);

    // Ajout des colonnes du second plateau

    double widthLittleColumn = 0.5;
    double heightLittleColumn = 4;

    double xLittleColumn1 = - Lx / 2 + rxLittleArcs / 2 - widthLittleColumn / 2;
    double yLittleColumn1 = - Ly + heightPlate2;

    double xLittleColumn2 = - rxBigArcs - widthPillier - 3 * rxLittleArcs / 2 - widthLittleColumn / 2;
    double yLittleColumn2 = - Ly + heightPlate2;

    double xLittleColumn3 = - rxBigArcs - widthPillier - rxLittleArcs / 2 - widthLittleColumn / 2;
    double yLittleColumn3 = - Ly + heightPlate2;

    double xLittleColumn4 = - rxBigArcs / 4 - widthLittleColumn / 2;
    double yLittleColumn4 = - Ly + heightPlate2;

    double xLittleColumn5 = - 3 * rxBigArcs / 4 - widthLittleColumn / 2;
    double yLittleColumn5 = - Ly + heightPlate2;

    double xLittleColumn6 =  rxBigArcs / 4 - widthLittleColumn / 2;
    double yLittleColumn6 = - Ly + heightPlate2;

    double xLittleColumn7 = 3 * rxBigArcs / 4 - widthLittleColumn / 2;
    double yLittleColumn7 = - Ly + heightPlate2;

    double xLittleColumn8 = rxBigArcs + widthPillier + rxLittleArcs / 2 - widthLittleColumn / 2;
    double yLittleColumn8 = - Ly + heightPlate2;

    double xLittleColumn9 = rxBigArcs + widthPillier + 3 *rxLittleArcs / 2 - widthLittleColumn / 2;
    double yLittleColumn9 = - Ly + heightPlate2;

    double xLittleColumn10 = Lx / 2 - rxLittleArcs / 2 - widthLittleColumn / 2;
    double yLittleColumn10 = - Ly + heightPlate2;

    int idLittleColumn1 = gmshModelOccAddRectangle(xLittleColumn1, yLittleColumn1, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn2 = gmshModelOccAddRectangle(xLittleColumn2, yLittleColumn2, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn3 = gmshModelOccAddRectangle(xLittleColumn3, yLittleColumn3, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn4 = gmshModelOccAddRectangle(xLittleColumn4, yLittleColumn4, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn5 = gmshModelOccAddRectangle(xLittleColumn5, yLittleColumn5, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn6 = gmshModelOccAddRectangle(xLittleColumn6, yLittleColumn6, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn7 = gmshModelOccAddRectangle(xLittleColumn7, yLittleColumn7, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn8 = gmshModelOccAddRectangle(xLittleColumn8, yLittleColumn8, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn9 = gmshModelOccAddRectangle(xLittleColumn9, yLittleColumn9, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);
    int idLittleColumn10 = gmshModelOccAddRectangle(xLittleColumn10, yLittleColumn10, 0.0, widthLittleColumn, heightLittleColumn, -1, 0, &ierr);

    // TODO : Ajouter les fils entre les deux grosses colonnes ici !

    
    // TODO : Ajouter aussi les fils verticaux entre les fils et le plateau de base ici !


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
    

    // TODO : Faire l'union entre les petites colonnes et les arcs ici !

    // On soustrait les arcs du plateau de base
    gmshModelOccCut(plate, 2, littleArc1, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, littleArc2, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, bigArc, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, littleArc3, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, littleArc4, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    // Garder qu'une des deux parties qui sont intersectées par les deux plateaux
    gmshModelOccFragment(plate2, 2, plate, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    

    geoSetSizeCallback(geoSize);
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
}