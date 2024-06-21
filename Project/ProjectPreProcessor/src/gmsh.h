#ifndef _GMSH_H_
#define _GMSH_H_

/*
Note:
The gmsh.h file is used to define the functions that are used to generate the mesh using Gmsh.
The functions are used to :
    * generate the mesh for the bridge, U example, and beam example.
    * set the domains' names and create the boundary conditions.
    * initialize and finalize the geometry.
    * import the mesh and generate the mesh.
*/

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