#ifndef _FEM_GMSH_H_
#define _FEM_GMSH_H_

#include "gmshc.h"
#include "fem_geometry.h"

double geoSizeDefault(double x, double y);
double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data);
void geoInitialize(void);
void geoFinalize(void);
void geoMeshImport(void);

void femErrorGmsh(int ierr, int line, char *file);

void geoMeshGenerate(void);
// void geoMeshGenerateGeo(void);
// void geoMeshGenerateGeoFile(const char *filename);
// void geoMeshGenerateMshFile(const char *filename);

#endif // _FEM_GMSH_H_