#include "fem.h"

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
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule  = theProblem->rule;
    femDiscrete *theSpace    = theProblem->space;
    femGeo *theGeometry      = theProblem->geometry;
    femNodes *theNodes       = theGeometry->theNodes;
    femMesh *theMesh         = theGeometry->theElements;

    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];

    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;
    double rho = theProblem->rho;
    double gx  = theProblem->gx;
    double gy  = theProblem->gy;
    double **A = theSystem->A;
    double *B  = theSystem->B;

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
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            double dxdxsi = 0.0; double dydxsi = 0.0;
            double dxdeta = 0.0; double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++)
            {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }

            double jac = dxdxsi * dydeta - dxdeta * dydxsi;
            if (jac < 0.0) { printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n"); }
            jac = fabs(jac); // Useless if the mesh is well oriented

            for (i = 0; i < theSpace->n; i++)
            {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }
            for (i = 0; i < theSpace->n; i++)
            {
                for (j = 0; j < theSpace->n; j++)
                {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight;
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight;
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight;
                }
                B[mapX[i]] += phi[i] * gx * rho * jac * weight;
                B[mapY[i]] += phi[i] * gy * rho * jac * weight;
            }
        }
    }
}

void femElasticityAssembleNeumann(femProblem *theProblem)
{
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule  = theProblem->ruleEdge;
    femDiscrete *theSpace    = theProblem->spaceEdge;
    femGeo *theGeometry      = theProblem->geometry;
    femNodes *theNodes       = theGeometry->theNodes;
    femMesh *theEdges        = theGeometry->theEdges;

    double x[2], y[2], phi[2];
    int iBnd, iElem, iInteg, iEdge, i, j, d, map[2], mapU[2];

    int nLocal = 2;
    double *B  = theSystem->B;

    for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++)
    {
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        double value = theCondition->value1;

        //
        // A completer
        // Attention, pour le normal tangent on calcule la normale (sortante) au SEGMENT, surtout PAS celle de constrainedNodes
        
        // Une petite aide pour le calcul de la normale :
        // double tx = theNodes->X[node1] - theNodes->X[node0];
        // double ty = theNodes->Y[node1] - theNodes->Y[node0];
        // double nx = ty;
        // double ny = -tx;
    }
}

void femElasticityApplyDirichlet(femProblem *theProblem)
{
    femFullSystem *theSystem = theProblem->system;
    femGeo *theGeometry      = theProblem->geometry;
    femNodes *theNodes       = theGeometry->theNodes;

    for (int node = 0; node < theNodes->nNodes; node++)
    {
        femConstrainedNode *theConstrainedNode = &theProblem->constrainedNodes[node];
        if (theConstrainedNode->type == UNDEFINED) { continue; }
        femBoundaryType type = theConstrainedNode->type;

        if (type == DIRICHLET_X)
        {
            double value = theConstrainedNode->value1;
            femFullSystemConstrain(theSystem, 2 * node + 0, value);
        }
        else if (type == DIRICHLET_Y)
        {
            double value = theConstrainedNode->value1;
            femFullSystemConstrain(theSystem, 2 * node + 1, value);
        }
        else if (type == DIRICHLET_XY)
        {
            double value_x = theConstrainedNode->value1;
            double value_y = theConstrainedNode->value2;
            femFullSystemConstrain(theSystem, 2 * node + 0, value_x);
            femFullSystemConstrain(theSystem, 2 * node + 1, value_y);
        }

        else if (type == DIRICHLET_N)
        {
            double value = theConstrainedNode->value1;
            double nx = theConstrainedNode->nx;
            double ny = theConstrainedNode->ny;
            // A completer
        }
        else if (type == DIRICHLET_T)
        {
            double value = theConstrainedNode->value1;
            double nx = theConstrainedNode->nx;
            double ny = theConstrainedNode->ny;
            // A completer
        }
        else if (type == DIRICHLET_NT)
        {
            double value_n = theConstrainedNode->value1;
            double value_t = theConstrainedNode->value2;
            double nx = theConstrainedNode->nx;
            double ny = theConstrainedNode->ny;
            // A completer
        }
    }
}

double *femElasticitySolve(femProblem *theProblem)
{
  femElasticityAssembleElements(theProblem);
  femElasticityAssembleNeumann(theProblem);
  femElasticityApplyDirichlet(theProblem);

  double *soluce = femFullSystemEliminate(theProblem->system);
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}