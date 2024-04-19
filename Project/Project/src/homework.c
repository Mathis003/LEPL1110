#include "fem.h"

void femElasticityAssembleElements(femProblem *theProblem)
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

            femSolverAssemble(theSolver, theProblem, mapX, mapY, phi, dphidx, dphidy, weightedJac, xLoc, theSpace->n);
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

    double xLoc, x[nLocal], y[nLocal], phi[nLocal], tx, ty, nx, ny, norm_n, norm_t;
    int iBnd, iElem, iInteg, iEdge, i, j, map[nLocal], mapU[nLocal];

    double **A = femSolverGetA(theSolver);
    double *B  = femSolverGetB(theSolver);

    for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++)
    {
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        femDomain *domain = theCondition->domain;
        double value = theCondition->value1 * FACTOR;
        
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
                        for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value * weightedJac; }
                    }
                    else if (type == NEUMANN_N)
                    {
                        tx = dx; ty = dy;
                        nx = ty; ny = -tx;
                        norm_n = sqrt(nx * nx + ny * ny);
                        nx /= norm_n; ny /= norm_n;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value * nx * weightedJac;
                            B[mapU[i] + 1] += phi[i] * value * ny * weightedJac;
                        }
                    }
                    else if (type == NEUMANN_T)
                    {
                        tx = dx; ty = dy;
                        norm_t = sqrt(tx * tx + ty * ty);
                        tx /= norm_t; ty /= norm_t;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value * tx * weightedJac;
                            B[mapU[i] + 1] += phi[i] * value * ty * weightedJac;
                        }
                    }
                }
                else if (theProblem->planarStrainStress == AXISYM)
                {
                    if (type == NEUMANN_X || type == NEUMANN_Y)
                    {
                        for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value * weightedJac * xLoc; }
                    }
                    else if (type == NEUMANN_N)
                    {
                        tx = dx; ty = dy;
                        nx = ty; ny = -tx;
                        norm_n = sqrt(nx * nx + ny * ny);
                        nx /= norm_n; ny /= norm_n;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value * nx * weightedJac * xLoc;
                            B[mapU[i] + 1] += phi[i] * value * ny * weightedJac * xLoc;
                        }
                    }
                    else if (type == NEUMANN_T)
                    {
                        tx = dx; ty = dy;
                        norm_t = sqrt(tx * tx + ty * ty);
                        tx /= norm_t; ty /= norm_t;
                        for (i = 0; i < theSpace->n; i++)
                        {
                            B[mapU[i]] += phi[i] * value * tx * weightedJac * xLoc;
                            B[mapU[i] + 1] += phi[i] * value * ty * weightedJac * xLoc;
                        }
                    }
                }
                else { Error("Unexpected planarStrainStress value !"); }
            }
        }
    }
}

void femElasticityApplyDirichlet(femProblem *theProblem, double FACTOR)
{
    femSolver *theSolver     = theProblem->solver;
    femGeometry *theGeometry = theProblem->geometry;
    femNodes *theNodes       = theGeometry->theNodes;

    int renumberNode;
    int *number = theNodes->number;
    int size = theSolver->size;

    for (int iNode = 0; iNode < theNodes->nNodes; iNode++)
    {
        femConstrainedNode *theConstrainedNode = &theProblem->constrainedNodes[iNode];
        if (theConstrainedNode->type == UNDEFINED) { continue; }
        femBoundaryType type = theConstrainedNode->type;

        renumberNode = number[iNode];

        int Ux = 2 * renumberNode;
        int Uy = 2 * renumberNode + 1;

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
            double nx, ny, tx, ty, myValue_n, myValue_t;

            // Already normalized
            nx = theConstrainedNode->nx;
            ny = theConstrainedNode->ny;
            tx = ny;
            ty = -nx;

            if (type == DIRICHLET_N)
            {
                myValue_n = theConstrainedNode->value1;
                femSolverSystemConstrainNT(theSolver, nx, ny, tx, ty, myValue_n * FACTOR, Ux, Uy);
            }
            else if (type == DIRICHLET_T)
            {
                myValue_t = theConstrainedNode->value1;
                femSolverSystemConstrainNT(theSolver, tx, ty, nx, ny, myValue_t * FACTOR, Ux, Uy);
            }
            else if (type == DIRICHLET_NT)
            {
                myValue_n = theConstrainedNode->value1;
                myValue_t = theConstrainedNode->value2;
                femSolverSystemConstrainNT(theSolver, nx, ny, tx, ty, myValue_n * FACTOR, Ux, Uy);
                femSolverSystemConstrainNT(theSolver, tx, ty, nx, ny, myValue_t * FACTOR, Ux, Uy);
            }
        }
    }
}

double *femElasticitySolve(femProblem *theProblem, femRenumType renumType, double FACTOR)
{
    femSolver *theSolver = theProblem->solver;
    femSolverInit(theSolver);
    femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem, FACTOR);
    femElasticityApplyDirichlet(theProblem, FACTOR);

    double *soluce = femSolverEliminate(theSolver);
    femNodes *theNodes = theProblem->geometry->theNodes;
    for (int i = 0; i < theNodes->nNodes; i++)
    {
        theProblem->soluce[2 * i] = soluce[2 * theNodes->number[i]];
        theProblem->soluce[2 * i + 1] = soluce[2 * theNodes->number[i] + 1];
    }
    return theProblem->soluce;
}