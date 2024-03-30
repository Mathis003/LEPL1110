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
    /*
    femMesh *theMesh = theProblem->geo->theElements;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
    femSolver *theSolver = theProblem->solver;
    int *number = theMesh->nodes->number;
    double source = theProblem->sourceValue;
    double dirichlet = theProblem->dirichletValue;

    if (theSpace->n > 4) Error("Unexpected discrete space size !");    
    double Xloc[4],Yloc[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    double Uloc[4];
    int iEdge,iElem,iInteg,i,j,map[4],ctr[4];
    double **A = theSolver->local->A;
    double *Aloc = theSolver->local->A[0];
    double *Bloc = theSolver->local->B;
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (i = 0; i < theSpace->n; i++)      Bloc[i] = 0;
        for (i = 0; i < (theSpace->n)*(theSpace->n); i++) Aloc[i] = 0;
        femDiffusionMeshLocal(theProblem,iElem,map,ctr,Xloc,Yloc,Uloc);  
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0; 
            double dydeta = 0;
            for (i = 0; i < theSpace->n; i++) {    
                dxdxsi += Xloc[i]*dphidxsi[i];       
                dxdeta += Xloc[i]*dphideta[i];   
                dydxsi += Yloc[i]*dphidxsi[i];   
                dydeta += Yloc[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[i][j] += (dphidx[i] * dphidx[j] 
                              + dphidy[i] * dphidy[j]) * jac * weight; }}                                                                                            
            for (i = 0; i < theSpace->n; i++) {
                Bloc[i] += phi[i] * jac * source * weight; }}
        for (i = 0; i < theSpace->n; i++) 
            if (ctr[i] == 1) femFullSystemConstrain(theSolver->local,i,dirichlet);
        femSolverAssemble(theSolver,Aloc,Bloc,Uloc,map,theSpace->n); } 
 
    double *soluce = femSolverEliminate(theSolver);
    for (i = 0; i < theProblem->size; i++)
        theProblem->soluce[i] += soluce[number[i]];
    */
    femSolver *theSolver     = theProblem->solver;
    femIntegration *theRule  = theProblem->rule;
    femDiscrete *theSpace    = theProblem->space;
    femGeometry *theGeometry = theProblem->geometry;
    femNodes *theNodes       = theGeometry->theNodes;
    femMesh *theMesh         = theGeometry->theElements;

    // femSolverAssemble(theSolver, ...);

    if (theSpace->n > 4) Error("Unexpected discrete space size !");

    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];

    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;
    double rho = theProblem->rho;
    double gx  = theProblem->gx;
    double gy  = theProblem->gy;
    double **A = theSolver->local->A;
    double *B  = theSolver->local->B;

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
    femSolver *theSolver = theProblem->solver;
    femIntegration *theRule  = theProblem->ruleEdge;
    femDiscrete *theSpace    = theProblem->spaceEdge;
    femGeometry *theGeometry      = theProblem->geometry;
    femNodes *theNodes       = theGeometry->theNodes;
    femMesh *theEdges        = theGeometry->theEdges;

    double x[2], y[2], phi[2];
    int iBnd, iElem, iInteg, iEdge, i, j, d, map[2], mapU[2];

    int nLocal = 2;
    double *B  = theSolver->local->B;

    for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++)
    {
        // Strip : BEGIN
        // TODO : Changer le code ci-dessous

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

        // // Skip Dirichlet boundary conditions
        // if (type == DIRICHLET_X || type == DIRICHLET_Y || type == DIRICHLET_XY ||
        //     type == DIRICHLET_N || type == DIRICHLET_T || type == DIRICHLET_NT) { continue; }
        
        // // SKIP FOR NOW (NOY YET IMPLMENTED)
        // if (type == NEUMANN_N || type == NEUMANN_T) { continue; }

        // // Iterate over the elements of the domain
        // for (iEdge = 0; iEdge < domain->nElem; iEdge++)
        // {
        //     // Get the element index (mapping)
        //     iElem = domain->elem[iEdge];

        //     // Mapping local nodes to global nodes
        //     for (j = 0; j < nLocal; j++)
        //     {
        //         map[j] = theEdges->elem[iElem * nLocal + j];
        //         mapU[j] = 2 * map[j] + 1;
        //         x[j] = theNodes->X[map[j]];
        //         y[j] = theNodes->Y[map[j]];
        //     }

        //     // Compute the constant Jacobian
        //     double dx = x[1] - x[0];
        //     double dy = y[1] - y[0];
        //     double jac = sqrt(dx * dx + dy * dy) / 2;

        //     // Iterate over the integration points
        //     for (iInteg = 0; iInteg < theRule->n; iInteg++)
        //     {
        //         // Get the integration point coordinates and weight
        //         double xsi    = theRule->xsi[iInteg];
        //         double weight = theRule->weight[iInteg];

        //         // Compute the shape functions
        //         femDiscretePhi(theSpace, xsi, phi);

        //         // Compute the forces and add them to the load vector
        //         for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value1 * jac * weight; }
        //     }
        // }
        // Strip : END
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
            femFullSystemConstrain(theSolver->local, 2 * node + 0, value);
        }
        else if (type == DIRICHLET_Y)
        {
            double value = theConstrainedNode->value1;
            femFullSystemConstrain(theSolver->local, 2 * node + 1, value);
        }
        else if (type == DIRICHLET_XY)
        {
            double value_x = theConstrainedNode->value1;
            double value_y = theConstrainedNode->value2;
            
            femFullSystemConstrain(theSolver->local, 2 * node + 0, value_x);
            femFullSystemConstrain(theSolver->local, 2 * node + 1, value_y);
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

    double *soluce = femSolverEliminate(theProblem->solver);
    memcpy(theProblem->soluce, soluce, theProblem->solver->local->size * sizeof(double));
    return theProblem->soluce;
}