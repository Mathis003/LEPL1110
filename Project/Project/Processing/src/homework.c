#include "../../../fem_library/include/fem.h"

/*
TODO :
    - Faire un solveur bande                  (IN PROGRESS)
    - Renumerotation des noeuds X/Y           (IN PROGRESS)
    - Verifier la validite de l'axisymetrique (IN PROGRESS)
    - Renumerotation des noeuds RCMK
    - Solveur frontal creux
    - Renumerotation des elements (pour solveur frontal)
    - Elements billineaire (fonction de forme du deuxieme degre)
    - Finir les conditions sur la geometrie du pont
    - Calcul des tensions aux noeuds (PAS COMPRIS => + FORCES ?)
    - Ecrire un script permettant de faire une animation des résultats

DONE :
    - Axisymetrique / tension plane / deformation plane
    - Conditions de Dirichlet (X, Y, XY, N, T, NT)
    - Conditions de Neumann (X, Y, N, T)
    - Résidus pour les forces (par noeud + globale)
    - Visualisation des résidus (Post-Processing avec 'X' et 'Y')
    - Visualisation de la matrice (Post-Processing avec 'S')
*/

void femElasticityAssembleElements(femProblem *theProblem, double FACTOR)
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

            femSolverAssemble(theSolver, theProblem, mapX, mapY, phi, dphidx, dphidy, weightedJac, xLoc, theSpace->n, FACTOR);
        }
    }
}

void femElasticityAssembleNeumann(femProblem *theProblem, double FACTOR)
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
                        for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value1 * weightedJac * FACTOR; }
                    }
                    else if (type == NEUMANN_N)
                    {
                        tx = dx; ty = dy;
                        nx = ty; ny = -tx;
                        norm_n = sqrt(nx * nx + ny * ny);
                        nx /= norm_n; ny /= norm_n;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value1 * nx * weightedJac * FACTOR;
                            B[mapU[i] + 1] += phi[i] * value1 * ny * weightedJac * FACTOR;
                        }
                    }
                    else if (type == NEUMANN_T)
                    {
                        tx = dx; ty = dy;
                        norm_t = sqrt(tx * tx + ty * ty);
                        tx /= norm_t; ty /= norm_t;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value1 * tx * weightedJac * FACTOR;
                            B[mapU[i] + 1] += phi[i] * value1 * ty * weightedJac * FACTOR;
                        }
                    }
                }
                else if (theProblem->planarStrainStress == AXISYM)
                {
                    if (type == NEUMANN_X || type == NEUMANN_Y)
                    {
                        for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value1 * weightedJac * xLoc * FACTOR; }
                    }
                    else if (type == NEUMANN_N)
                    {
                        tx = dx; ty = dy;
                        nx = ty; ny = -tx;
                        norm_n = sqrt(nx * nx + ny * ny);
                        nx /= norm_n; ny /= norm_n;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value1 * nx * weightedJac * xLoc * FACTOR;
                            B[mapU[i] + 1] += phi[i] * value1 * ny * weightedJac * xLoc * FACTOR;
                        }
                    }
                    else if (type == NEUMANN_T)
                    {
                        tx = dx; ty = dy;
                        norm_t = sqrt(tx * tx + ty * ty);
                        tx /= norm_t; ty /= norm_t;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value1 * tx * weightedJac * xLoc * FACTOR;
                            B[mapU[i] + 1] += phi[i] * value1 * ty * weightedJac * xLoc * FACTOR;
                        }
                    }
                }
                else { Error("Unexpected planarStrainStress value !"); }
            }
        }
    }
}

void getValueConditionsNT(double vect_x, double vect_y, double value, int Ux, int Uy, int *node1, int *node2, double *a, double *b)
{
    if (fabs(vect_x) >= fabs(vect_y))
    {
        *a = - vect_y / vect_x;
        *b = value / vect_x;
        *node1 = Ux; *node2 = Uy;
    }
    else
    {
        *a = - vect_x / vect_y;
        *b = value / vect_y;
        *node1 = Uy; *node2 = Ux;
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

        int Ux = 2 * node;
        int Uy = 2 * node + 1;

        if (type == DIRICHLET_X)
        {
            double value = theConstrainedNode->value1;
            femSolverSystemConstrainXY(theSolver, Ux, value);
        }
        else if (type == DIRICHLET_Y)
        {
            double value = theConstrainedNode->value1;
            femSolverSystemConstrainXY(theSolver, Uy, value);
        }
        else if (type == DIRICHLET_XY)
        {
            double value_x = theConstrainedNode->value1;
            double value_y = theConstrainedNode->value2;
            
            femSolverSystemConstrainXY(theSolver, Ux, value_x);
            femSolverSystemConstrainXY(theSolver, Uy, value_y);
        }
        else
        {
            double nx, ny, norm_n, tx, ty, norm_t, a_n, b_n, value_n, a_t, b_t, value_t;
            int node1, node2;

            value_n = theConstrainedNode->value1;
            value_t = theConstrainedNode->value2;
            nx = theConstrainedNode->nx;
            ny = theConstrainedNode->ny;

            if (type == DIRICHLET_N)
            {
                norm_n = sqrt(nx * nx + ny * ny);
                nx /= norm_n; ny /= norm_n;
                getValueConditionsNT(nx, ny, value_n, Ux, Uy, &node1, &node2, &a_n, &b_n);
                femSolverSystemConstrainNT(theSolver, node1, node2, a_n, b_n);
            }
            else if (type == DIRICHLET_T)
            {
                tx = ny; ty = -nx;
                norm_t = sqrt(tx * tx + ty * ty);
                tx /= norm_t; ty /= norm_t;
                getValueConditionsNT(tx, ty, value_t, Ux, Uy, &node1, &node2, &a_t, &b_t);
                femSolverSystemConstrainNT(theSolver, node1, node2, a_t, b_t);
            }
            else if (type == DIRICHLET_NT)
            {
                tx = ny; ty = -nx;
                norm_n = sqrt(nx * nx + ny * ny);
                norm_t = sqrt(tx * tx + ty * ty);
                nx /= norm_n; ny /= norm_n;
                tx /= norm_t; ty /= norm_t;
                getValueConditionsNT(nx, ny, value_n, Ux, Uy, &node1, &node2, &a_n, &b_n);
                femSolverSystemConstrainNT(theSolver, node1, node2, a_n, b_n);
                getValueConditionsNT(tx, ty, value_t, Ux, Uy, &node1, &node2, &a_t, &b_t);
                femSolverSystemConstrainNT(theSolver, node1, node2, a_t, b_t);
            }
        }
    }
}

double *femElasticitySolve(femProblem *theProblem, double FACTOR)
{
    femSolver *theSolver = theProblem->solver;
    femSolverInit(theSolver);
    femElasticityAssembleElements(theProblem, FACTOR);
    femElasticityAssembleNeumann(theProblem, FACTOR);

    double **A = getMatrixA(theSolver);
    double *B  = getVectorB(theSolver);
    int size   = getSizeMatrix(theSolver);
    
    // Copy the Dirichlet unconstrained system
    femSystemWrite(A, B, size, "../data/dirichletUnconstrainedSystem.txt");

    femElasticityApplyDirichlet(theProblem);

    A = getMatrixA(theSolver);
    B = getVectorB(theSolver);
    size = getSizeMatrix(theSolver);

    // Copy the final system
    femSystemWrite(A, B, size, "../data/finalSystem.txt");

    double *soluce = femSolverEliminate(theSolver);
    memcpy(theProblem->soluce, soluce, size * sizeof(double));
    return theProblem->soluce;
}