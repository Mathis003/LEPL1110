#ifndef _FEM_INTEGRATE_H_
#define _FEM_INTEGRATE_H_

#include <stdio.h>
#include <stdlib.h>

#include "fem_other.h"

static const double _gaussQuad4Xsi[4]    = {-0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Eta[4]    = {0.577350269189626, -0.577350269189626, -0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Weight[4] = {1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000};
static const double _gaussTri3Xsi[3]     = {0.166666666666667, 0.666666666666667, 0.166666666666667};
static const double _gaussTri3Eta[3]     = {0.166666666666667, 0.166666666666667, 0.666666666666667};
static const double _gaussTri3Weight[3]  = {0.166666666666667, 0.166666666666667, 0.166666666666667};

typedef struct {
    int n;
    femElementType type;
    void (*x2) (double *xsi, double *eta);
    void (*phi2) (double xsi, double eta, double *phi);
    void (*dphi2dx) (double xsi, double eta, double *dphidxsi, double *dphideta);
    void (*x) (double *xsi);
    void (*phi) (double xsi, double *phi);
    void (*dphidx) (double xsi, double *dphidxsi);
} femDiscrete;

typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

femIntegration *femIntegrationCreate(int n, femElementType type);
void femIntegrationFree(femIntegration *theRule);

femDiscrete *femDiscreteCreate(int n, femElementType type);
void femDiscreteFree(femDiscrete *mySpace);
void femDiscretePrint(femDiscrete *mySpace);
void femDiscreteXsi2(femDiscrete *mySpace, double *xsi, double *eta);
void femDiscretePhi2(femDiscrete *mySpace, double xsi, double eta, double *phi);
void femDiscreteDphi2(femDiscrete *mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);
void femDiscreteXsi(femDiscrete *mySpace, double *xsi);
void femDiscretePhi(femDiscrete *mySpace, double xsi, double *phi);
void femDiscreteDphi(femDiscrete *mySpace, double xsi, double *dphidxsi);

#endif // _FEM_INTEGRATE_H_