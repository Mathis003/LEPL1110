#ifndef _GMSH_H_
#define _GMSH_H_

#include "fem.h"
#include "gmshc.h"

double geoSizeDefault(double x, double y);
double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data);
void geoInitialize(void);
void geoFinalize(void);
void geoMeshImport(femDiscreteType discreteType);
void femErrorGmsh(int ierr, int line, char *file);

// Defined in Project/PreProcessing/homework.c
void geoMeshGenerate(femDiscreteType discreteType);
void geoMeshGenerateExample(femDiscreteType discreteType);
void setDomainsName(void);
void createBoundaryConditions(femProblem *theProblem);

#endif // _GMSH_H_