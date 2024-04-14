#include "../../../fem_library/include/fem.h"

/*
TODO :
    - Faire un solveur bande                  (IN PROGRESS)
    - Verifier la validite de l'axisymetrique (IN PROGRESS)
    - Elements billineaire (fonction de forme du deuxieme degre)   (IN PROGRESS)
    - Renumerotation des noeuds RCMK
    - Solveur frontal creux
    - Renumerotation des elements (pour solveur frontal)
    - Finir les conditions sur la geometrie du pont
    - Calcul des tensions aux noeuds (PAS COMPRIS => + FORCES ?)

DONE :
    - Axisymetrique / tension plane / deformation plane
    - Conditions de Dirichlet (X, Y, XY, N, T, NT)
    - Conditions de Neumann (X, Y, N, T)
    - Résidus pour les forces (par noeud + globale)
    - Visualisation des résidus (Post-Processing avec 'X' et 'Y')
    - Visualisation de la matrice (Post-Processing avec 'S')
    - Script d'animation
    - Renumerotation des noeuds X/Y
*/

void femElasticityAssembleElements(femProblem *theProblem, double FACTOR)
{
    femSolver *theSolver     = theProblem->solver;
    femIntegration *theRule  = theProblem->rule;
    femDiscrete *theSpace    = theProblem->space;
    femGeometry *theGeometry = theProblem->geometry;
    femNodes *theNodes       = theGeometry->theNodes;
    femMesh *theMesh         = theGeometry->theElements;

    int nLocal = theSpace->n;
    int *number = theMesh->nodes->number;

    double xLoc, x[nLocal], y[nLocal], phi[nLocal], dphidxsi[nLocal], dphideta[nLocal], dphidx[nLocal], dphidy[nLocal];
    int iElem, iInteg, iEdge, i, map[nLocal], mapX[nLocal], mapY[nLocal];

    for (iElem = 0; iElem < theMesh->nElem; iElem++)
    {
        for (i = 0; i < theSpace->n; i++)
        {
            map[i] = theMesh->elem[iElem * nLocal + i];
            x[i] = theNodes->X[map[i]];
            y[i] = theNodes->Y[map[i]];
            map[i] = number[map[i]];
            mapX[i] = 2 * map[i];
            mapY[i] = 2 * map[i] + 1;
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
            xLoc = 0.0;
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

    int nLocal  = theSpace->n;
    int *number = theNodes->number;

    double xLoc, x[nLocal], y[nLocal], phi[nLocal], tx, ty, nx, ny, norm_n, norm_t, a, b;
    int iBnd, iElem, iInteg, iEdge, i, j, d, map[nLocal], mapU[nLocal];
    double oldValue, oldType;
    int changeType;

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

            for (j = 0; j < nLocal; j++)
            {
                map[j] = theEdges->elem[iElem * nLocal + j];
                x[j] = theNodes->X[map[j]];
                y[j] = theNodes->Y[map[j]];
                map[j] = number[map[j]];
                mapU[j] = nLocal * map[j] + shift;
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

                xLoc = 0.0;
                for (j = 0; i < theSpace->n; i++) { xLoc += phi[j] * x[j]; }

                if (theProblem->planarStrainStress == PLANAR_STRAIN || theProblem->planarStrainStress == PLANAR_STRESS)
                {
                    if (type == NEUMANN_X || type == NEUMANN_Y)
                    {
                        for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value1 * FACTOR * weightedJac; }
                    }
                    else if (type == NEUMANN_N)
                    {
                        tx = dx; ty = dy;
                        nx = ty; ny = -tx;
                        norm_n = sqrt(nx * nx + ny * ny);
                        nx /= norm_n; ny /= norm_n;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value1 * FACTOR * nx * weightedJac;
                            B[mapU[i] + 1] += phi[i] * value1 * FACTOR * ny * weightedJac;
                        }
                    }
                    else if (type == NEUMANN_T)
                    {
                        tx = dx; ty = dy;
                        norm_t = sqrt(tx * tx + ty * ty);
                        tx /= norm_t; ty /= norm_t;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value1 * FACTOR * tx * weightedJac;
                            B[mapU[i] + 1] += phi[i] * value1 * FACTOR * ty * weightedJac;
                        }
                    }
                }
                else if (theProblem->planarStrainStress == AXISYM)
                {
                    if (type == NEUMANN_X || type == NEUMANN_Y)
                    {
                        for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value1 * FACTOR * weightedJac * xLoc; }
                    }
                    else if (type == NEUMANN_N)
                    {
                        tx = dx; ty = dy;
                        nx = ty; ny = -tx;
                        norm_n = sqrt(nx * nx + ny * ny);
                        nx /= norm_n; ny /= norm_n;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value1 * FACTOR * nx * weightedJac * xLoc;
                            B[mapU[i] + 1] += phi[i] * value1 * FACTOR * ny * weightedJac * xLoc;
                        }
                    }
                    else if (type == NEUMANN_T)
                    {
                        tx = dx; ty = dy;
                        norm_t = sqrt(tx * tx + ty * ty);
                        tx /= norm_t; ty /= norm_t;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value1 * FACTOR * tx * weightedJac * xLoc;
                            B[mapU[i] + 1] += phi[i] * value1 * FACTOR * ty * weightedJac * xLoc;
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

void femElasticityApplyDirichlet(femProblem *theProblem, double FACTOR)
{
    femSolver *theSolver     = theProblem->solver;
    femGeometry *theGeometry = theProblem->geometry;
    femNodes *theNodes       = theGeometry->theNodes;

    int *number = theNodes->number;

    for (int iNode = 0; iNode < theNodes->nNodes; iNode++)
    {
        femConstrainedNode *theConstrainedNode = &theProblem->constrainedNodes[iNode];
        if (theConstrainedNode->type == UNDEFINED) { continue; }
        femBoundaryType type = theConstrainedNode->type;

        iNode = number[iNode];

        int Ux = 2 * iNode;
        int Uy = 2 * iNode + 1;

        if (type == DIRICHLET_X)
        {
            double value = theConstrainedNode->value1;
            femSolverSystemConstrainXY(theSolver, Ux, value * FACTOR);
        }
        else if (type == DIRICHLET_Y)
        {
            double value = theConstrainedNode->value1;
            femSolverSystemConstrainXY(theSolver, Uy, value * FACTOR);
        }
        else if (type == DIRICHLET_XY)
        {
            double value_x = theConstrainedNode->value1;
            double value_y = theConstrainedNode->value2;
            
            femSolverSystemConstrainXY(theSolver, Ux, value_x * FACTOR);
            femSolverSystemConstrainXY(theSolver, Uy, value_y * FACTOR);
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
                getValueConditionsNT(nx, ny, value_n * FACTOR, Ux, Uy, &node1, &node2, &a_n, &b_n);
                femSolverSystemConstrainNT(theSolver, node1, node2, a_n, b_n);
            }
            else if (type == DIRICHLET_T)
            {
                tx = ny; ty = -nx;
                norm_t = sqrt(tx * tx + ty * ty);
                tx /= norm_t; ty /= norm_t;
                getValueConditionsNT(tx, ty, value_t * FACTOR, Ux, Uy, &node1, &node2, &a_t, &b_t);
                femSolverSystemConstrainNT(theSolver, node1, node2, a_t, b_t);
            }
            else if (type == DIRICHLET_NT)
            {
                tx = ny; ty = -nx;
                norm_n = sqrt(nx * nx + ny * ny);
                norm_t = sqrt(tx * tx + ty * ty);
                nx /= norm_n; ny /= norm_n;
                tx /= norm_t; ty /= norm_t;
                getValueConditionsNT(nx, ny, value_n * FACTOR, Ux, Uy, &node1, &node2, &a_n, &b_n);
                femSolverSystemConstrainNT(theSolver, node1, node2, a_n, b_n);
                getValueConditionsNT(tx, ty, value_t * FACTOR, Ux, Uy, &node1, &node2, &a_t, &b_t);
                femSolverSystemConstrainNT(theSolver, node1, node2, a_t, b_t);
            }
        }
    }
}

double *femElasticitySolve(femProblem *theProblem, femRenumType renumType, double FACTOR)
{
    femSolver *theSolver = theProblem->solver;
    femSolverInit(theSolver);

    femElasticityAssembleElements(theProblem, FACTOR);

    femElasticityAssembleNeumann(theProblem, FACTOR);

    // Copy the Dirichlet unconstrained system
    // femSystemWrite(theSolver, "../data/dirichletUnconstrainedSystem.txt");

    femElasticityApplyDirichlet(theProblem, FACTOR);

    for (int i = 0; i < 10; i++)
    {
        int start = i;
        int end = i + 10;
        for (int j = start; j < end; j++)
        {
            printf("elem (%d, %d) : %f\n", i, j, femSolverGet(theSolver, i, j));
        }
    }

    // Copy the final system
    // femSystemWrite(theSolver, "../data/finalSystem.txt");

    double *soluce = femSolverEliminate(theSolver);
    femNodes *theNodes = theProblem->geometry->theNodes;
    for (int i = 0; i < theNodes->nNodes; i++)
    {
        theProblem->soluce[2 * i] = soluce[2 * theNodes->number[i]];
        theProblem->soluce[2 * i + 1] = soluce[2 * theNodes->number[i] + 1];
    }

    return theProblem->soluce;
}