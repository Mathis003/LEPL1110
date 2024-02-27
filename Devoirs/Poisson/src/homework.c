#include "fem.h"

# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femGeo *theGeometry = geoMeshCreate(filename);
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));

    theProblem->geo  = theGeometry;
    femMesh *theMesh = theGeometry->theElements;

    if (theMesh->nLocalNode == 4)
    {
        theProblem->space = femDiscreteCreate(4, FEM_QUAD);
        theProblem->rule  = femIntegrationCreate(4, FEM_QUAD);
    }
    else if (theMesh->nLocalNode == 3)
    {
        theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
        theProblem->rule  = femIntegrationCreate(3, FEM_TRIANGLE);
    }

    theProblem->system = femFullSystemCreate(theMesh->nodes->nNodes);
    return theProblem;
}

# endif
# ifndef NOPOISSONBOUNDARY

void femPoissonFindBoundaryNodes(femPoissonProblem *theProblem)
{
    femGeo *theGeometry = theProblem->geo;  
    femMesh *theEdges = theGeometry->theEdges; 
    int nBoundary = 0;
    
    // Strip : BEGIN
    nBoundary = theEdges->nElem * theEdges->nLocalNode;
    // Strip : END

    femDomain *theBoundary = malloc(sizeof(femDomain));
    theGeometry->nDomains++;
    theGeometry->theDomains = realloc(theGeometry->theDomains, theGeometry->nDomains * sizeof(femDomain*));
    theGeometry->theDomains[theGeometry->nDomains - 1] = theBoundary;
    theBoundary->nElem = nBoundary;
    theBoundary->elem = malloc(nBoundary * sizeof(int));

    // Strip : BEGIN
    for (int i = 0; i < nBoundary; i++) { theBoundary->elem[i] = theEdges->elem[i]; }
    // Strip : END

    theBoundary->mesh = NULL;
    sprintf(theBoundary->name, "Boundary");
}
    
# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
    // Strip : BEGIN
    geoMeshFree(theProblem->geo);
    femDiscreteFree(theProblem->space);
    femIntegrationFree(theProblem->rule);
    femFullSystemFree(theProblem->system);
    free(theProblem);
    // Strip : END
}
    
# endif
# ifndef NOPOISSONLOCAL

void femPoissonLocal(femPoissonProblem *theProblem, const int iElem, int *map, double *x, double *y)
{
    femMesh *theMesh = theProblem->geo->theElements;
    
    // Strip : BEGIN
    for (int i = 0; i < theMesh->nLocalNode; i++)
    {
        map[i] = theMesh->elem[i + iElem * theMesh->nLocalNode];
        x[i]   = theMesh->nodes->X[map[i]];
        y[i]   = theMesh->nodes->Y[map[i]];
    }
    // Strip : END
}

# endif
# ifndef NOPOISSONSOLVE


void femPoissonSolve(femPoissonProblem *theProblem)
{
    // Extract necessary components from the problem struct
    femMesh *theMesh         = theProblem->geo->theElements;
    femDomain *theBoundary   = geoGetDomain(theProblem->geo,"Boundary");
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule  = theProblem->rule;
    femDiscrete *theSpace    = theProblem->space;
    
    // Check if the discrete space size is as expected
    if (theSpace->n > 4) { Error("Unexpected discrete space size !"); }

    // Arrays to store local coordinates, shape functions, and derivatives
    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int iElem, iInteg, iEdge, i, j, map[4];
    int nLocal = theMesh->nLocalNode;

    // Strip : BEGIN

    // Iterate over the elements
    for (iElem = 0; iElem < theMesh->nElem; iElem++)
    {
        // Get the local nodes (global -> local mapping)
        femPoissonLocal(theProblem, iElem, map, x, y);

        // Iterate over the integration points
        for (iInteg = 0; iInteg < theRule->n; iInteg++)
        {
            // Get the integration point coordinates and weight
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            // Compute the shape functions and their derivatives
            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
            
            // Compute the Jacobian components
            double dxdxsi = 0; double dxdeta = 0;
            double dydxsi = 0; double dydeta = 0;
            for (i = 0; i < theSpace->n; i++) 
            {
                dxdxsi += x[i] * dphidxsi[i];  dxdeta += x[i] * dphideta[i];   
                dydxsi += y[i] * dphidxsi[i];  dydeta += y[i] * dphideta[i];
            }

            // Compute the Jacobian
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

            // Compute the derivatives of the shape functions
            for (i = 0; i < theSpace->n; i++)
            { 
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            // Compute the weighted Jacobian
            double weightedJac = jac * weight;
            
            // local -> global mapping for the matrix and vector
            for (i = 0; i < theSpace->n; i++)
            { 
                for (j = 0; j < theSpace->n; j++)
                {
                    // Assemble the local matrix A_e with the global matrix A
                    theSystem->A[map[i]][map[j]] += (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * weightedJac;
                }
                // Assemble the local vector B_e with the global vector B
                theSystem->B[map[i]] += phi[i] * weightedJac;
            }
            // For Debugging
            // femFullSystemPrint(theSystem);
        }
    }

    // Apply boundary conditions (0.0 on the boundary)
    for (int i = 0; i < theBoundary->nElem; i++) { femFullSystemConstrain(theSystem, theBoundary->elem[i], 0.0); }

    // For Debugging
    // femFullSystemPrint(theSystem);

    // Solve the linear system
    femFullSystemEliminate(theSystem);

    // Strip : END
}

# endif