#ifndef _FEM_GEOMETRY_H_
#define _FEM_GEOMETRY_H_

#define ErrorGmsh(a) femErrorGmsh(a, __LINE__, __FILE__)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gmshc.h"

#include "fem_other.h"

void geoInitialize(void);
femGeo *geoGetGeometry(void);
double geoSize(double x, double y);
double geoSizeDefault(double x, double y);
double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data);
void geoInitialize(void);
void geoFree(void);
void geoFinalize(void);
void geoSetSizeCallback(double (*geoSize)(double x, double y));
void geoMeshImport(void);
void geoMeshGenerate(void);
void geoMeshPrint(void);
void geoMeshWrite(const char *filename);
void geoMeshRead(const char *filename);
void geoSetDomainName(int iDomain, char *name);
int geoGetDomain(char *name);

void femErrorGmsh(int test, int line, char *file);

#endif // _FEM_GEOMETRY_H_