#include "../include/fem_integrate.h"

femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (theRule == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    if (type == FEM_QUAD && n == 4)
    {
        theRule->n      = 4;
        theRule->xsi    = _gaussQuad4Xsi;
        theRule->eta    = _gaussQuad4Eta;
        theRule->weight = _gaussQuad4Weight;
    }
    else if (type == FEM_TRIANGLE && n == 3)
    {
        theRule->n      = 3;
        theRule->xsi    = _gaussTri3Xsi;
        theRule->eta    = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight;
    }
    else { Error("Cannot create such an integration rule !"); }
    return theRule;
}

void femIntegrationFree(femIntegration *theRule) { free(theRule); }

void _q1c0_x(double *xsi, double *eta)
{
    xsi[0] = 1.0;  eta[0] = 1.0;
    xsi[1] = -1.0; eta[1] = 1.0;
    xsi[2] = -1.0; eta[2] = -1.0;
    xsi[3] = 1.0;  eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;
    phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
    phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
    phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] =  (1.0 + eta) / 4.0; dphidxsi[1] = -(1.0 + eta) / 4.0;
    dphidxsi[2] = -(1.0 - eta) / 4.0; dphidxsi[3] =  (1.0 - eta) / 4.0;
    dphideta[0] =  (1.0 + xsi) / 4.0; dphideta[1] =  (1.0 - xsi) / 4.0;
    dphideta[2] = -(1.0 - xsi) / 4.0; dphideta[3] = -(1.0 + xsi) / 4.0;
}

void _p1c0_x(double *xsi, double *eta)
{
    xsi[0] = 0.0; eta[0] = 0.0;
    xsi[1] = 1.0; eta[1] = 0.0;
    xsi[2] = 0.0; eta[2] = 1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;
    phi[1] = xsi;
    phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = -1.0;
    dphidxsi[1] =  1.0;
    dphidxsi[2] =  0.0;
    dphideta[0] = -1.0;
    dphideta[1] =  0.0;
    dphideta[2] =  1.0;
}

femDiscrete *femDiscreteCreate(int n, femElementType type)
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete));
    if (theSpace == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    if (type == FEM_QUAD && n == 4)
    {
        theSpace->n       = 4;
        theSpace->x2      = _q1c0_x;
        theSpace->phi2    = _q1c0_phi;
        theSpace->dphi2dx = _q1c0_dphidx;
    }
    else if (type == FEM_TRIANGLE && n == 3)
    {
        theSpace->n       = 3;
        theSpace->x2      = _p1c0_x;
        theSpace->phi2    = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx;
    }
    else { Error("Cannot create such a discrete space !"); }
    return theSpace;
}

void femDiscreteFree(femDiscrete *theSpace) { free(theSpace); }

void femDiscreteXsi2(femDiscrete *mySpace, double *xsi, double *eta) { mySpace->x2(xsi, eta); }

void femDiscretePhi2(femDiscrete *mySpace, double xsi, double eta, double *phi) { mySpace->phi2(xsi, eta, phi); }

void femDiscreteDphi2(femDiscrete *mySpace, double xsi, double eta, double *dphidxsi, double *dphideta) { mySpace->dphi2dx(xsi, eta, dphidxsi, dphideta); }

void femDiscreteXsi(femDiscrete *mySpace, double *xsi) { mySpace->x(xsi); }

void femDiscretePhi(femDiscrete *mySpace, double xsi, double *phi) { mySpace->phi(xsi, phi); }

void femDiscreteDphi(femDiscrete *mySpace, double xsi, double *dphidxsi) { mySpace->dphidx(xsi, dphidxsi); }

void femDiscretePrint(femDiscrete *mySpace)
{
    double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];
    int n = mySpace->n;

    femDiscreteXsi2(mySpace, xsi, eta);
    for (int i = 0; i < n; i++)
    {
        femDiscretePhi2(mySpace, xsi[i], eta[i], phi);
        femDiscreteDphi2(mySpace, xsi[i], eta[i], dphidxsi, dphideta);

        for (int j = 0; j < n; j++)
        {
            printf("(xsi = %+.1f, eta = %+.1f) : ", xsi[i], eta[i]);
            printf(" phi(%d) = %+.1f", j, phi[j]);
            printf("   dphidxsi(%d) = %+.1f", j, dphidxsi[j]);
            printf("   dphideta(%d) = %+.1f \n", j, dphideta[j]);
        }
        printf(" \n");
    }
}