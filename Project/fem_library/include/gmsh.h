#ifndef _GMSH_H_
#define _GMSH_H_

#include "fem.h"
#include "gmshc.h"

double geoSizeDefault(double x, double y);
double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data);
void geoInitialize(void);
void geoFinalize(void);
void geoMeshImport(void);
void femErrorGmsh(int ierr, int line, char *file);

// Defined in Project/PreProcessing/homework.c
void geoMeshGenerate(void);
void geoMeshGenerateExample(void);
void setDomainsName(void);
void createBoundaryConditions(femProblem *theProblem);

#endif // _GMSH_H_