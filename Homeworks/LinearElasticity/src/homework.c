#include "fem.h"

void geoMeshGenerate()
{
    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    
    int ierr;
    double r = w / 4;
    int idRect = gmshModelOccAddRectangle(0.0, 0.0, 0.0, w, h, -1, 0.0, &ierr); 
    int idDisk = gmshModelOccAddDisk(w / 2.0, h / 2.0, 0.0, r, r, -1, NULL, 0, NULL, 0, &ierr); 
    int idSlit = gmshModelOccAddRectangle(w / 2.0, h / 2.0 - r, 0.0, w, 2.0 * r, -1, 0.0, &ierr); 
    int rect[] = {2, idRect};
    int disk[] = {2, idDisk};
    int slit[] = {2, idSlit};

    gmshModelOccCut(rect, 2, disk, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); 
    gmshModelOccCut(rect, 2, slit, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); 
    gmshModelOccSynchronize(&ierr); 

    if (theGeometry->elementType == FEM_QUAD)
    {
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.Algorithm", 11, &ierr);  
        gmshOptionSetNumber("Mesh.SmoothRatio", 21.5, &ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr); 
        gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);  
        gmshModelMeshGenerate(2, &ierr);
    }
  
    if (theGeometry->elementType == FEM_TRIANGLE)
    {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);
    }

    return;
}

double *femElasticitySolve(femProblem *theProblem)
{
    femFullSystem  *theSystem   = theProblem->system;
    femIntegration *theRule     = theProblem->rule;
    femDiscrete    *theSpace    = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes    = theGeometry->theNodes;
    femMesh        *theMesh     = theGeometry->theElements;

    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int iElem, iInteg, i, j, map[4], mapX[4], mapY[4];

    int nLocal = theMesh->nLocalNode;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;

    // Strip : BEGIN

    // Iterate over the elements
    for (iElem = 0; iElem < theMesh->nElem; iElem++)
    {
        // Get mapping local nodes
        for (i = 0; i < nLocal; i++)
        {
            map[i]  = theMesh->elem[i + iElem * nLocal];
            x[i]    = theNodes->X[map[i]];
            y[i]    = theNodes->Y[map[i]];
            mapX[i] = 2 * map[i];
            mapY[i] = 2 * map[i] + 1;
        }

        // Iterate over the integration points
        for (iInteg = 0; iInteg < theRule->n; iInteg++)
        {   
            // Get the integration point coordinates and weight
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            // Compute the shape functions and their derivatives
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);

            // Compute the Jacobian components
            double dxdxsi = 0.0; double dxdeta = 0.0;
            double dydxsi = 0.0; double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++)
            {
                dxdxsi += x[i] * dphidxsi[i]; dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i]; dydeta += y[i] * dphideta[i];
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
                // Assemble the local matrix A_e with the global matrix A
                for(j = 0; j < theSpace->n; j++)
                {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * weightedJac;
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * weightedJac;
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * weightedJac;
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * weightedJac;
                }
                // Assemble the local vector B_e with the global vector B
                B[mapY[i]] -= g * rho * phi[i] * weightedJac;
            }
        }
    }
    // Strip : End                
  
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (i = 0; i < theSystem->size; i++)
    {
        if (theConstrainedNodes[i] != -1)
        {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem, i, value);
        }
    }
                            
    return femFullSystemEliminate(theSystem);
}