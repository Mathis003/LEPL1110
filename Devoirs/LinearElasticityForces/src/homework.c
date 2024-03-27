#include "fem.h"

double **A_copy = NULL;
double *B_copy  = NULL;

void femElasticityAssembleElements(femProblem *theProblem)
{
    femFullSystem  *theSystem   = theProblem->system;
    femIntegration *theRule     = theProblem->rule;
    femDiscrete    *theSpace    = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes    = theGeometry->theNodes;
    femMesh        *theMesh     = theGeometry->theElements;
    femMesh        *theEdges    = theGeometry->theEdges;

    double x[4], y[4], phi[4], dphidxsi[4] ,dphideta[4], dphidx[4], dphidy[4];
    int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];
    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++)
    {
        for (j = 0; j < nLocal; j++)
        {
            map[j]  = theMesh->elem[iElem * nLocal + j];
            mapX[j] = 2 * map[j];
            mapY[j] = 2 * map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
        } 
        
        for (iInteg = 0; iInteg < theRule->n; iInteg++)
        {   
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
            
            double dxdxsi = 0.0; double dxdeta = 0.0;
            double dydxsi = 0.0; double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++)
            {  
                dxdxsi += x[i] * dphidxsi[i];       
                dxdeta += x[i] * dphideta[i];   
                dydxsi += y[i] * dphidxsi[i];   
                dydeta += y[i] * dphideta[i];
            }

            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++)
            {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            double weightedJac = jac * weight;

            for (i = 0; i < theSpace->n; i++)
            { 
                for (j = 0; j < theSpace->n; j++)
                {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * weightedJac;
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * weightedJac;
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * weightedJac;
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * weightedJac;
                }
                B[mapY[i]] -= phi[i] * g * rho * weightedJac;
            }
        }
    }
}


void femElasticityAssembleNeumann(femProblem *theProblem)
{
    femFullSystem  *theSystem   = theProblem->system;
    femIntegration *theRule     = theProblem->ruleEdge;
    femDiscrete    *theSpace    = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes    = theGeometry->theNodes;
    femMesh        *theEdges    = theGeometry->theEdges;

    double x[2], y[2], phi[2];
    int iBnd, iElem, iInteg, iEdge, i, j, d, map[2], mapU[2];
    int nLocal = 2;
    double *B  = theSystem->B;

    for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++)
    {
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        double value = theCondition->value;
    
        // Strip : BEGIN
        if (type == DIRICHLET_X || type == DIRICHLET_Y) { continue; }
        
        for (iEdge = 0; iEdge < theEdges->nElem; iEdge++)
        {
            if (theEdges->elem[iEdge * nLocal] != type) { continue; }

            for (j = 0; j < nLocal; j++)
            {
                map[j] = theEdges->elem[iEdge * nLocal + j];
                mapU[j] = 2 * map[j] + 1;
                x[j] = theNodes->X[map[j]];
                y[j] = theNodes->Y[map[j]];
            }

            for (iInteg = 0; iInteg < theRule->n; iInteg++)
            {
                double xsi    = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];

                femDiscretePhi(theSpace, xsi, phi);

                double dx = x[1] - x[0];
                double dy = y[1] - y[0];
                double jac = sqrt(dx * dx + dy * dy);

                for (i = 0; i < theSpace->n; i++)
                {
                    B[mapU[i]] += phi[i] * value * jac * weight;
                }
            }
        }
        // Strip : END
    }
}


// Strip : BEGIN
double *femElasticitySolve(femProblem *theProblem)
{
    // Assembly of stiffness matrix and load vector
    femElasticityAssembleElements(theProblem);

    // Assembly of Neumann boundary conditions
    femElasticityAssembleNeumann(theProblem);

    femFullSystem *theSystem = theProblem->system;
    int size = theSystem->size;

    // Allocate memory for the copy of the stiffness matrix A and the load vector B
    A_copy = (double **) malloc(sizeof(double *) * size);
    for (int i = 0; i < size; i++) { A_copy[i] = (double *) malloc(sizeof(double) * size); }
    B_copy = (double *) malloc(sizeof(double) * size);

    // Copy the stiffness matrix A and the load vector B
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++) { A_copy[i][j] = theSystem->A[i][j]; }
        B_copy[i] = theSystem->B[i];
    }

    // Apply Dirichlet boundary conditions (costraints the nodes)
    int *theConstrainedNodes = theProblem->constrainedNodes;
    for (int i = 0; i < size; i++)
    {
        if (theConstrainedNodes[i] != -1)
        {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem, i, value);
        }
    }

    // Solve the system and return the solution
    theProblem->soluce = femFullSystemEliminate(theSystem);
    
    return theProblem->soluce;
}
// Strip : END


// Strip : BEGIN
double *femElasticityForces(femProblem *theProblem)
{
    double *residuals = theProblem->residuals;
    double *soluce    = theProblem->soluce;
    int size = theProblem->system->size;

    // Allocate memory for residuals if not already done
    if (residuals == NULL) { residuals = (double *) malloc(sizeof(double) * size); }

    // Initialize residuals to zero
    for (int i = 0; i < size; i++) { residuals[i] = 0.0; }

    /*
    Compute residuals: R = A * U - B where A and B are the system matrix
    and load vector before applying Dirichlet boundary conditions
    */
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++) { residuals[i] += A_copy[i][j] * soluce[j]; }
        residuals[i] -= B_copy[i];

        // Invert residuals to get forces (action-reaction principle)
        residuals[i] *= -1;
    }

    // Free memory
    free(A_copy);
    free(B_copy);

    // Return the forces
    return residuals;
}
// Strip : END