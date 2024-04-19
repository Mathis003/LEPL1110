#ifndef _GMSH_H_
#define _GMSH_H_

#include "../../Project/src/fem.h"

#include "gmshc.h"

double geoSizeDefault(double x, double y);
double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data);
void geoInitialize(void);
void geoFinalize(void);
void geoMeshImport(femDiscreteType discreteType);
void femErrorGmsh(int ierr, int line, char *file);

void geoMeshGenerate(femDiscreteType discreteType, int bridgeSimplified);
void geoMeshGenerateExample(femDiscreteType discreteType, int beam_example);
void geoMeshGenerate_UForm(femDiscreteType discreteType);
void geoMeshGenerate_Beam(femDiscreteType discreteType);
void setDomainsName(int bridgeSimplified);
void createBoundaryConditions(femProblem *theProblem, int bridgeSimplified);

#endif // _GMSH_H_