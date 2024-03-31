#ifndef _FEM_GMSH_H_
#define _FEM_GMSH_H_

#include "gmshc.h"
#include "fem_geometry.h"
#include "fem_elasticity.h"

double geoSizeDefault(double x, double y);
double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data);
void geoInitialize(void);
void geoFinalize(void);
void geoMeshImport(void);
void femErrorGmsh(int ierr, int line, char *file);

// Defined in Project/PreProcessing/homework.c
void geoMeshGenerate(void);
int createRectangle(double x, double y, double width, double height);
int createDisk(double xc, double yc, double rx, double ry);
void cutElement(int *mainElement, int *cutElement);
void fuseElement(int *mainElement, int *fuseElement);
void rotateElement(int *element, double posX, double posY, double angle);

void setDomainsName(void);
void createBoundaryConditions(femProblem *theProblem);
// int geoMeshRefine(int nRefine); // TODO

#endif // _FEM_GMSH_H_