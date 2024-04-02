#include "../../../fem_library/include/fem_elasticity.h"

/*
TODO
Il faut generaliser ce code :
    - Ajouter l'axisymetrique (obligatoire)
    - Ajouter les conditions de Neumann (obligatoire)
    - Ajouter les conditions en normal et tangentiel (fortement conseille)
    - Remplacer le solveur plein par un solveur plus subtil (obligatoire)
        => Solveur bande + renumerotation / solveur frontal / solveur multigrid / ...
*/

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
    double **A, *B;

    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;
    double rho = theProblem->rho;
    double gx  = theProblem->gx;
    double gy  = theProblem->gy;

    if (theSolver->type == FEM_FULL)
    {
        femFullSystem *theFullSystem = (femFullSystem *) theSolver->solver;
        A = theFullSystem->A;
        B = theFullSystem->B;
    }
    else if (theSolver->type == FEM_BAND)
    {
        femBandSystem *theBandSystem = (femBandSystem *) theSolver->solver;
        A = theBandSystem->A;
        B = theBandSystem->B;
        size = theBandSystem->size;
    }
    else if (theSolver->type == FEM_ITER)
    {
        femFullSystem *theFullSystem = theSolver->local;
        A = theFullSystem->A;
        B = theFullSystem->B;
        size = theFullSystem->size;
    }
    else { Error("Unexpected solver type !"); }

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
                xLoc += x[i] * phi[i];
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

            if (theProblem->planarStrainStress == PLANAR_STRAIN)
            {
                for (i = 0; i < theSpace->n; i++)
                {
                    for (j = 0; j < theSpace->n; j++)
                    {
                        if (theSolver->type == FEM_FULL || theSolver->type == FEM_ITER)
                        {
                            A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * weightedJac;
                            A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * weightedJac;
                            A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * weightedJac;
                            A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * weightedJac;
                        }
                        else if (theSolver->type == FEM_BAND)
                        {
                            A[mapX[i]][mapX[j]] += (mapX[j] >= mapX[i]) ? (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * weightedJac : 0.0;
                            A[mapX[i]][mapY[j]] += (mapY[j] >= mapX[i]) ? (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * weightedJac : 0.0;
                            A[mapY[i]][mapX[j]] += (mapX[j] >= mapY[i]) ? (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * weightedJac : 0.0;
                            A[mapY[i]][mapY[j]] += (mapY[j] >= mapY[i]) ? (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * weightedJac : 0.0;
                        }
                        else { Error("Unexpected solver type !"); }

                    B[mapX[i]] += phi[i] * gx * rho * jac * weight;
                    B[mapY[i]] += phi[i] * gy * rho * jac * weight;
                    }
                }
            }
            else if (theProblem->planarStrainStress == PLANAR_STRESS)
            {
                // TODO : A completer (pareil que PLANAR_STRAIN mais avec des termes supplementaires ?)
            }
            else if (theProblem->planarStrainStress == AXISYM)
            {
                for (i = 0; i < theSpace->n; i++)
                {
                    for (j = 0; j < theSpace->n; j++)
                    {
                        if (theSolver->type == FEM_FULL || theSolver->type == FEM_ITER)
                        {
                            A[mapX[i]][mapX[j]] += (dphidx[i] * a * xLoc * dphidx[j] + dphidy[i] * c * xLoc * dphidy[j] + phi[i] * ((b * dphidx[j]) + (a * phi[j] / xLoc)) + dphidx[i] * b * phi[j]) * jac * weight;
                            A[mapX[i]][mapY[j]] += (dphidx[i] * b * xLoc * dphidy[j] + dphidy[i] * c * xLoc * dphidx[j] + phi[i] * b * dphidy[j]) * jac * weight;
                            A[mapY[i]][mapX[j]] += (dphidy[i] * b * xLoc * dphidx[j] + dphidx[i] * c * xLoc * dphidy[j] + dphidy[i] * b * phi[j]) * jac * weight;
                            A[mapY[i]][mapY[j]] += (dphidy[i] * a * xLoc * dphidy[j] + dphidx[i] * c * xLoc * dphidx[j]) * jac * weight;
                        }
                        else if (theSolver->type == FEM_BAND)
                        {
                            A[mapX[i]][mapX[j]] += (mapX[j] >= mapX[i]) ? (dphidx[i] * a * xLoc * dphidx[j] + dphidy[i] * c * xLoc * dphidy[j] + phi[i] * ((b * dphidx[j]) + (a * phi[j] / xLoc)) + dphidx[i] * b * phi[j]) * jac * weight : 0.0;
                            A[mapX[i]][mapY[j]] += (mapY[j] >= mapX[i]) ? (dphidx[i] * b * xLoc * dphidy[j] + dphidy[i] * c * xLoc * dphidx[j] + phi[i] * b * dphidy[j]) * jac * weight : 0.0;
                            A[mapY[i]][mapX[j]] += (mapX[j] >= mapY[i]) ? (dphidy[i] * b * xLoc * dphidx[j] + dphidx[i] * c * xLoc * dphidy[j] + dphidy[i] * b * phi[j]) * jac * weight : 0.0;
                            A[mapY[i]][mapY[j]] += (mapY[j] >= mapY[i]) ? (dphidy[i] * a * xLoc * dphidy[j] + dphidx[i] * c * xLoc * dphidx[j]) * jac * weight : 0.0;
                        }
                        else { Error("Unexpected solver type !"); }
                    }
                    B[mapX[i]] -= phi[i] * xLoc * gx * rho * jac * weight;
                    B[mapY[i]] -= phi[i] * xLoc * gy * rho * jac * weight;
                }
            }
            else { Error("Unexpected planarStrainStress value !"); }
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
    double *B  = theSolver->local->B;

    for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++)
    {
        // Attention, pour le normal tangent on calcule la normale (sortante) au SEGMENT, surtout PAS celle de constrainedNodes
        
        // Une petite aide pour le calcul de la normale :
        // double tx = theNodes->X[node1] - theNodes->X[node0];
        // double ty = theNodes->Y[node1] - theNodes->Y[node0];
        // double nx = ty;
        // double ny = -tx;

        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        femDomain *domain = theCondition->domain;
        double value1 = theCondition->value1;
        double value2 = theCondition->value2;

        if (type == NEUMANN_N || type == NEUMANN_T) 
        {
            // TODO : Implement the Neumann boundary conditions (normal and tangent)
            continue;
        }
        else if (type == NEUMANN_X || type == NEUMANN_Y)
        {
            for (iEdge = 0; iEdge < domain->nElem; iEdge++)
            {
                iElem = domain->elem[iEdge];

                for (j = 0; j < nLocal; j++)
                {
                    map[j] = theEdges->elem[iElem * nLocal + j];
                    mapU[j] = 2 * map[j] + 1;
                    x[j] = theNodes->X[map[j]];
                    y[j] = theNodes->Y[map[j]];
                }

                double dx = x[1] - x[0];
                double dy = y[1] - y[0];
                double length = sqrt(dx * dx + dy * dy);
                double jac = length / 2;
-
                for (iInteg = 0; iInteg < theRule->n; iInteg++)
                {
                    double xsi    = theRule->xsi[iInteg];
                    double weight = theRule->weight[iInteg];

                    femDiscretePhi(theSpace, xsi, phi);

                    if (theProblem->planarStrainStress == AXISYM)
                    {
                        // TODO : A completer
                    } else if (theProblem->planarStrainStress == PLANAR_STRAIN)
                    {
                        // TODO : A completer
                        if (theSolver->type == FEM_FULL || theSolver->type == FEM_ITER)
                        {
                            // TODO : A completer
                            for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value1 * jac * weight; }
                        }
                        else if (theSolver->type == FEM_BAND)
                        {
                            // TODO : A completer
                            for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value1 * jac * weight; }
                        }
                        else { Error("Unexpected solver type !"); }
                    }
                    else if (theProblem->planarStrainStress == PLANAR_STRESS)
                    {
                        // TODO : A completer
                    }
                    else { Error("Unexpected planarStrainStress value !"); }
                }
            }
        }
        else { continue; }
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

            // A completer (COMPRENDS RIEN)
            // double tx = theNodes->X[node + 1] - theNodes->X[node];
            // double ty = theNodes->Y[node + 1] - theNodes->Y[node];
            // double nx = ty;
            // double ny = -tx;
        }
        else if (type == DIRICHLET_T)
        {
            double value = theConstrainedNode->value1;
            double nx = theConstrainedNode->nx;
            double ny = theConstrainedNode->ny;
            double tx = ny;
            double ty = -nx;
            // A completer
        }
        else if (type == DIRICHLET_NT)
        {
            double value_n = theConstrainedNode->value1;
            double value_t = theConstrainedNode->value2;
            double nx = theConstrainedNode->nx;
            double ny = theConstrainedNode->ny;
            double tx = ny;
            double ty = -nx;
            // A completer
        }
    }
}

double *femElasticitySolve(femProblem *theProblem)
{
    femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem);
    femElasticityApplyDirichlet(theProblem);

    double *soluce = femSolverEliminate(theProblem->solver);
    memcpy(theProblem->soluce, soluce, theProblem->solver->local->size * sizeof(double));
    return theProblem->soluce;
}