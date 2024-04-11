#include "../../../fem_library/include/fem.h"

/*
TODO :
    - Faire un solveur bande avec renumerotation RCMK ET solveur frontal creux
    - Ajouter les conditions de Dirichlets (N, T, NT)

DONE :
    - Ajouter l'axisymetrique
    - Ajouter les conditions de Dirichlet (x, y, xy)
    - Ajouter les conditions de Neumann (x, y, n, t)
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

    double x[2], y[2], phi[2], tx, ty, nx, ny, norm_n, norm_t, a, b;
    int iBnd, iElem, iInteg, iEdge, i, j, d, map[2], mapU[2];
    double oldValue, oldType;
    int changeType;

    int nLocal = 2;
    double **A = getMatrixA(theSolver);
    double *B  = getVectorB(theSolver);

    for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++)
    {
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        femDomain *domain = theCondition->domain;
        double value1 = theCondition->value1;
        
        int shift = 0;
        if (type == NEUMANN_Y)  { shift = 1; }
        else if (type == DIRICHLET_X || type == DIRICHLET_Y || type == DIRICHLET_N || type == DIRICHLET_T || type == DIRICHLET_NT) { continue; }

        for (iEdge = 0; iEdge < domain->nElem; iEdge++)
        {
            iElem = domain->elem[iEdge];

            double xLoc = 0.0;
            for (j = 0; j < nLocal; j++)
            {
                map[j] = theEdges->elem[iElem * nLocal + j];
                mapU[j] = 2 * map[j] + shift;
                x[j] = theNodes->X[map[j]];
                y[j] = theNodes->Y[map[j]];
                xLoc += phi[j] * x[j];
            }

            double dx = x[1] - x[0];
            double dy = y[1] - y[0];
            double length = sqrt(dx * dx + dy * dy);
            double jac = length / 2.0;

            for (iInteg = 0; iInteg < theRule->n; iInteg++)
            {
                double xsi    = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];

                double weightedJac = jac * weight;

                femDiscretePhi(theSpace, xsi, phi);

                if (theProblem->planarStrainStress == PLANAR_STRAIN || theProblem->planarStrainStress == PLANAR_STRESS)
                {
                    if (type == NEUMANN_X || type == NEUMANN_Y)
                    {
                        for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value1 * weightedJac; }
                    }
                    else if (type == NEUMANN_N)
                    {
                        tx = dx; ty = dy;
                        nx = ty; ny = -tx;
                        norm_n = sqrt(nx * nx + ny * ny);
                        nx /= norm_n; ny /= norm_n;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value1 * nx * weightedJac;
                            B[mapU[i] + 1] += phi[i] * value1 * ny * weightedJac;
                        }
                    }
                    else if (type == NEUMANN_T)
                    {
                        tx = dx; ty = dy;
                        norm_t = sqrt(tx * tx + ty * ty);
                        tx /= norm_t; ty /= norm_t;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value1 * tx * weightedJac;
                            B[mapU[i] + 1] += phi[i] * value1 * ty * weightedJac;
                        }
                    }
                }
                else if (theProblem->planarStrainStress == AXISYM)
                {
                    if (type == NEUMANN_X || type == NEUMANN_Y)
                    {
                        for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value1 * weightedJac * xLoc; }
                    }
                    else if (type == NEUMANN_N)
                    {
                        tx = dx; ty = dy;
                        nx = ty; ny = -tx;
                        norm_n = sqrt(nx * nx + ny * ny);
                        nx /= norm_n; ny /= norm_n;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value1 * nx * weightedJac * xLoc;
                            B[mapU[i] + 1] += phi[i] * value1 * ny * weightedJac * xLoc;
                        }
                    }
                    else if (type == NEUMANN_T)
                    {
                        tx = dx; ty = dy;
                        norm_t = sqrt(tx * tx + ty * ty);
                        tx /= norm_t; ty /= norm_t;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value1 * tx * weightedJac * xLoc;
                            B[mapU[i] + 1] += phi[i] * value1 * ty * weightedJac * xLoc;
                        }
                    }
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

    double **A = getMatrixA(theSolver);
    double *B  = getVectorB(theSolver);
    int size   = getSizeMatrix(theSolver);

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

            double norm_n = sqrt(nx * nx + ny * ny);
            nx /= norm_n; ny /= norm_n;

            int Ux = 2 * node;
            int Uy = 2 * node + 1;

            double a_n, b_n;
            int node1, node2;

            // nx * Ux + ny * Uy = value_n
            // <=> Uy = (value_n/ny) - (nx/ny) * Ux    (if fabs(ny) >= fabs(nx)
            // <=> Ux = (value_n/nx) - (ny/nx) * Uy    (if fabs(nx) >= fabs(ny)
            if (fabs(nx) >= fabs(ny))
            {
                a_n = - ny / nx;
                b_n = value / nx;
                node1 = Ux;
                node2 = Uy;
            }
            else
            {
                a_n = - nx / ny;
                b_n = value / ny;
                node1 = Uy;
                node2 = Ux;
            }

            for (int i = 0; i < size; i++)
            {
                if (i == node2) { A[node2][i] = 1.0; }
                else if (i == node1) { A[node2][i] = - a_n; }
                else { A[node2][i] = 0.0; }
            }
            B[node2] = b_n;
        }
        else if (type == DIRICHLET_T)
        {
            double value = theConstrainedNode->value1;
            double nx = theConstrainedNode->nx;
            double ny = theConstrainedNode->ny;
            double tx = ny;
            double ty = -nx;

            double norm_t = sqrt(tx * tx + ty * ty);
            tx /= norm_t; ty /= norm_t;

            int Ux = 2 * node;
            int Uy = 2 * node + 1;

            double a_t, b_t;
            int node1, node2;

            // tx * Ux + ty * Uy = value_t
            // <=> Uy = (value_t/ty) - (tx/ty) * Ux    (if fabs(ty) >= fabs(tx)
            // <=> Ux = (value_t/tx) - (ty/tx) * Uy    (if fabs(tx) >= fabs(ty)
            if (fabs(tx) >= fabs(ty))
            {
                a_t = - ty / tx;
                b_t = value / tx;
                node1 = Ux;
                node2 = Uy;
            }
            else
            {
                a_t = - tx / ty;
                b_t = value / ty;
                node1 = Uy;
                node2 = Ux;
            }

            for (int i = 0; i < size; i++)
            {
                if (i == node2) { A[node2][i] = 1.0; }
                else if (i == node1) { A[node2][i] = - a_t; }
                else { A[node2][i] = 0.0; }
            }
            B[node2] = b_t;
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

            nx /= norm_n; ny /= norm_n;
            tx /= norm_t; ty /= norm_t;

            int Ux = 2 * node;
            int Uy = 2 * node + 1;

            // nx * Ux + ny * Uy = value_n
            // <=> Uy = (value_n/ny) - (nx/ny) * Ux    (if fabs(ny) >= fabs(nx)
            // <=> Ux = (value_n/nx) - (ny/nx) * Uy    (if fabs(nx) >= fabs(ny)
            double a_n, b_n;
            int node1, node2;
            if (fabs(nx) >= fabs(ny))
            {
                a_n = - ny / nx;
                b_n = value_n / nx;
                node1 = Ux;
                node2 = Uy;
            }
            else
            {
                a_n = - nx / ny;
                b_n = value_n / ny;
                node1 = Uy;
                node2 = Ux;
            }

            for (int i = 0; i < size; i++)
            {
                if (i == node2) { A[node2][i] = 1.0; }
                else if (i == node1) { A[node2][i] = - a_n; }
                else { A[node2][i] = 0.0; }
            }
            B[node2] = b_n;


            // tx * Ux + ty * Uy = value_t
            // <=> Uy = (value_t/ty) - (tx/ty) * Ux    (if fabs(ty) >= fabs(tx)
            // <=> Ux = (value_t/tx) - (ty/tx) * Uy    (if fabs(tx) >= fabs(ty)
            double a_t, b_t;
            if (fabs(tx) >= fabs(ty))
            {
                a_t = - ty / tx;
                b_t = value_t / tx;
                node1 = Ux;
                node2 = Uy;
            }
            else
            {
                a_t = - tx / ty;
                b_t = value_t / ty;
                node1 = Uy;
                node2 = Ux;
            }

            for (int i = 0; i < size; i++)
            {
                if (i == node2) { A[node2][i] = 1.0; }
                else if (i == node1) { A[node2][i] = - a_t; }
                else { A[node2][i] = 0.0; }
            }
            B[node2] = b_t;
        }
    }
}

double *femElasticitySolve(femProblem *theProblem)
{
    femSolver *theSolver = theProblem->solver;
    femSolverInit(theSolver);
    femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem);

    double **A = getMatrixA(theSolver);
    double *B  = getVectorB(theSolver);
    int size   = getSizeMatrix(theSolver);

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

    double *soluce = femSolverEliminate(theSolver);
    memcpy(theProblem->soluce, soluce, size * sizeof(double));
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
