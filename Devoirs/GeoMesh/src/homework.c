#include "fem.h"

// Strip : BEGIN
double hermiteInterpolation(double x, double h_max, double h_min, double dist_interp)
{
    // if (x >= dist_interp) return h_max;
    // if (x <= 0)           return h_min;

    // double dist_interp_squared = dist_interp * dist_interp;
    // double dist_interp_cubed = dist_interp_squared * dist_interp;

    // double a = - 2 * (h_max - h_min) / dist_interp_cubed;
    // double b = 3 * (h_max - h_min) / dist_interp_squared;
    // double d = h_min;
    // // int c = 0

    // double x_squared = x * x;
    // double x_cubed = x_squared * x;

    // double h = a * x_cubed + b * x_squared + d;
    // return h;

    if (x >= dist_interp) return h_max;
    if (x <= 0)           return h_min;
    return (1 / (dist_interp * dist_interp)) * (h_max - h_min) * x * x * (- 2 * x / dist_interp + 3) + h_min;
}
// Strip : END

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
    double d_Notch = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
    double d_Hole  = sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1));

    double h_Notch = hermiteInterpolation(d_Notch - r0, h, h0, d0);
    double h_Hole  = hermiteInterpolation(d_Hole - r1, h, h1, d1);

    h = fmin(h_Notch, h_Hole);
    // Strip : END
     
    return h;
}


#define ___ 0

void geoMeshGenerate()
{
    femGeo* theGeometry = geoGetGeometry();

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
