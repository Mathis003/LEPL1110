#ifndef _FEM_GEOMETRY_H_
#define _FEM_GEOMETRY_H_

#define ErrorGmsh(a) femErrorGmsh(a, __LINE__, __FILE__)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fem_other.h"

extern femGeometry theGeometry;

femGeometry *geoGetGeometry(void);
double geoSize(double x, double y);
void geoFree(void);
void geoSetSizeCallback(double (*geoSize)(double x, double y));
void geoMeshPrint(void);
void geoMeshWrite(const char *filename);
void geoMeshRead(const char *filename);
void geoSetDomainName(int iDomain, char *name);
int geoGetDomain(char *name);

#endif // _FEM_GEOMETRY_H_