#include "../../../fem_library/include/fem.h"

/*
TODO
Il faut generaliser ce code :
    - Ajouter l'axisymetrique (obligatoire)
    - Ajouter les conditions de Neumann (obligatoire)
    - Ajouter les conditions en normal et tangentiel (fortement conseille)
    - Remplacer le solveur plein par un solveur plus subtil (obligatoire)
        => Solveur bande + renumerotation / solveur frontal / solveur multigrid / ...
*/

double **A_copy = NULL;
double *B_copy  = NULL;

void femElasticityAssembleElements(femProblem *theProblem)
{
    femSolver *theSolver     = theProblem->solver;
    femIntegration *theRule  = theProblem->rule;
    femDiscrete *theSpace    = theProblem->space;
    femGeometry *theGeometry = theProblem->geometry;
    femNodes *theNodes       = theGeometry->theNodes;
    femMesh *theMesh         = theGeometry->theElements;

    if (theSpace->n > 4) Error("Unexpected discrete space size !");

    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int size, iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];
    void *theSystem;
    double **A, *B;

    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;
    double rho = theProblem->rho;
    double gx  = theProblem->gx;
    double gy  = theProblem->gy;

    A     = getMatrixA(theSolver);
    B     = getVectorB(theSolver);
    size = getSizeMatrix(theSolver);

    for (iElem = 0; iElem < theMesh->nElem; iElem++)
    {
        for (j = 0; j < nLocal; j++)
        {
            map[j] = theMesh->elem[iElem * nLocal + j];
            mapX[j] = 2 * map[j];
            mapY[j] = 2 * map[j] + 1;
            x[j] = theNodes->X[map[j]];
            y[j] = theNodes->Y[map[j]];
        }

        for (iInteg = 0; iInteg < theRule->n; iInteg++)
        {
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            double dxdxsi = 0.0; double dydxsi = 0.0;
            double dxdeta = 0.0; double dydeta = 0.0;
            double xLoc = 0.0;
            for (i = 0; i < theSpace->n; i++)
            {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
                xLoc   += x[i] * phi[i];
            }

            double jac = dxdxsi * dydeta - dxdeta * dydxsi;
            if (jac < 0.0) { printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n"); }
            jac = fabs(jac);

            for (i = 0; i < theSpace->n; i++)
            {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            double weightedJac = jac * weight;

            femSolverAssemble(theSolver, theProblem, mapX, mapY, phi, dphidx, dphidy, weightedJac, xLoc, theSpace->n);
        }
    }
}

void femElasticityAssembleNeumann(femProblem *theProblem)
{
    femSolver *theSolver     = theProblem->solver;
    femIntegration *theRule  = theProblem->ruleEdge;
    femDiscrete *theSpace    = theProblem->spaceEdge;
    femGeometry *theGeometry = theProblem->geometry;
    femNodes *theNodes       = theGeometry->theNodes;
    femMesh *theEdges        = theGeometry->theEdges;

    double x[2], y[2], phi[2];
    int iBnd, iElem, iInteg, iEdge, i, j, d, map[2], mapU[2];

    int nLocal = 2;
    double *B  = getVectorB(theSolver);

    for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++)
    {
        // Attention, pour le normal tangent on calcule la normale (sortante) au SEGMENT, surtout PAS celle de constrainedNodes
        
        // Une petite aide pour le calcul de la normale :
        // double tx = theNodes->X[node1] - theNodes->X[node0];
        // double ty = theNodes->Y[node1] - theNodes->Y[node0];
        // double nx = ty;
        // double ny = -tx;

        // A completer (COMPRENDS RIEN)
        // double tx = theNodes->X[node + 1] - theNodes->X[node];
        // double ty = theNodes->Y[node + 1] - theNodes->Y[node];
        // double nx = ty;
        // double ny = -tx;

        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        femDomain *domain = theCondition->domain;
        double value1 = theCondition->value1;
        // double value2 = theCondition->value2;

        if (type == DIRICHLET_X || type == DIRICHLET_Y || type == DIRICHLET_XY || type == DIRICHLET_N || type == DIRICHLET_T || type == DIRICHLET_NT) { continue; }
        if (type == NEUMANN_N || type == NEUMANN_T) { continue; } // TODO : A PRENDRE EN COMPTE APRES

        int shift;
        if (type == NEUMANN_X)      { shift = 0; }
        else if (type == NEUMANN_Y) { shift = 1; }

        for (iEdge = 0; iEdge < domain->nElem; iEdge++)
        {
            iElem = domain->elem[iEdge];

            for (j = 0; j < nLocal; j++)
            {
                map[j] = theEdges->elem[iElem * nLocal + j];
                mapU[j] = 2 * map[j] + shift;
                x[j] = theNodes->X[map[j]];
                y[j] = theNodes->Y[map[j]];
            }

            double dx = x[1] - x[0];
            double dy = y[1] - y[0];
            double length = sqrt(dx * dx + dy * dy);
            double jac = length / 2;

            for (iInteg = 0; iInteg < theRule->n; iInteg++)
            {
                double xsi    = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];

                femDiscretePhi(theSpace, xsi, phi);

                if (theProblem->planarStrainStress == PLANAR_STRAIN || theProblem->planarStrainStress == PLANAR_STRESS)
                {
                    for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value1 * jac * weight; }
                }
                else if (theProblem->planarStrainStress == AXISYM)
                {
                    for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value1 * jac * weight * x[i]; }
                }
                else { Error("Unexpected planarStrainStress value !"); }
            }
        }
    }
}

void femElasticityApplyDirichlet(femProblem *theProblem)
{
    femSolver *theSolver     = theProblem->solver;
    femGeometry *theGeometry = theProblem->geometry;
    femNodes *theNodes       = theGeometry->theNodes;

    for (int node = 0; node < theNodes->nNodes; node++)
    {
        femConstrainedNode *theConstrainedNode = &theProblem->constrainedNodes[node];
        if (theConstrainedNode->type == UNDEFINED) { continue; }
        femBoundaryType type = theConstrainedNode->type;

        if (type == DIRICHLET_X)
        {
            double value = theConstrainedNode->value1;
            femSolverSystemConstrain(theSolver, 2 * node + 0, value);
        }
        else if (type == DIRICHLET_Y)
        {
            double value = theConstrainedNode->value1;
            femSolverSystemConstrain(theSolver, 2 * node + 1, value);
        }
        else if (type == DIRICHLET_XY)
        {
            double value_x = theConstrainedNode->value1;
            double value_y = theConstrainedNode->value2;
            
            femSolverSystemConstrain(theSolver, 2 * node + 0, value_x);
            femSolverSystemConstrain(theSolver, 2 * node + 1, value_y);
        }
        else if (type == DIRICHLET_N)
        {
            double value = theConstrainedNode->value1;
            double nx = theConstrainedNode->nx;
            double ny = theConstrainedNode->ny;

            //We normalise the normal vector
            double norm = sqrt(nx * nx + ny * ny);
            nx /= norm;
            ny /= norm;

            femSolverSystemConstrain(theSolver, 2 * node + 0, value * nx);
            femSolverSystemConstrain(theSolver, 2 * node + 1, value * ny);
        }
        else if (type == DIRICHLET_T)
        {
            double value = theConstrainedNode->value1;
            double nx = theConstrainedNode->nx;
            double ny = theConstrainedNode->ny;
            double tx = ny;
            double ty = -nx;

            double norm = sqrt(tx * tx + ty * ty);
            tx /= norm;
            ty /= norm;

            femSolverSystemConstrain(theSolver, 2 * node + 0, value * tx);
            femSolverSystemConstrain(theSolver, 2 * node + 1, value * ty);
        }
        else if (type == DIRICHLET_NT)
        {
            double value_n = theConstrainedNode->value1;
            double value_t = theConstrainedNode->value2;
            double nx = theConstrainedNode->nx;
            double ny = theConstrainedNode->ny;
            double tx = ny;
            double ty = -nx;
            
            double norm_n = sqrt(nx * nx + ny * ny);
            double norm_t = sqrt(tx * tx + ty * ty);

            //We normalise the normal and tangent vectors
            nx /= norm_n;
            ny /= norm_n;
            tx /= norm_t;
            ty /= norm_t;

            femSolverSystemConstrain(theSolver, 2 * node + 0, value_n * nx + value_t * tx);
            femSolverSystemConstrain(theSolver, 2 * node + 1, value_n * ny + value_t * ty);
        }
    }
}

double *femElasticitySolve(femProblem *theProblem)
{
    femSolver *theSolver = theProblem->solver;
    femSolverInit(theSolver);
    femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem);

    int size;
    double **A, *B;
    if (theSolver->type == FEM_FULL)
    {
        femFullSystem *theSystem = (femFullSystem *) theSolver;
        A = theSystem->A;
        B = theSystem->B;
        size = theSystem->size;
    }
    else if (theSolver->type == FEM_BAND)
    {
        femBandSystem *theSystem = (femBandSystem *) theSolver;
        A = theSystem->A;
        B = theSystem->B;
        size = theSystem->size;
    }
    else { Error("Unexpected solver type !"); }

    if (A_copy == NULL)
    {
        A_copy = (double **) malloc(sizeof(double *) * size);
        for (int i = 0; i < size; i++) { A_copy[i] = (double *) malloc(sizeof(double) * size); }
    }
    if (B_copy == NULL) { B_copy = (double *) malloc(sizeof(double) * size); }

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++) { A_copy[i][j] = A[i][j]; }
        B_copy[i] = B[i];
    }

    femElasticityApplyDirichlet(theProblem);

    double *soluce = femSolverEliminate(theProblem->solver);
    memcpy(theProblem->soluce, soluce, getSizeMatrix(theProblem->solver) * sizeof(double));
    return theProblem->soluce;
}

double *femElasticityForces(femProblem *theProblem)
{
    double *residuals = theProblem->residuals;
    double *soluce    = theProblem->soluce;
    int size          = getSizeMatrix(theProblem->solver);

    if (residuals == NULL) { residuals = (double *) malloc(sizeof(double) * size); }
    for (int i = 0; i < size; i++) { residuals[i] = 0.0; }

    /*
    Compute residuals: R = A * U - B where A and B are the system matrix
    and load vector before applying Dirichlet boundary conditions.
    */
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++) { residuals[i] += A_copy[i][j] * soluce[j]; }
        residuals[i] -= B_copy[i];
    }

    for (int i = 0; i < size; i++) { free(A_copy[i]); A_copy[i] = NULL;}
    free(A_copy); free(B_copy);
    A_copy = NULL; B_copy = NULL;

    return residuals;
}