/*
*  main.c
*  Projet 2023-2024
*  Elasticite lineaire plane
*
*  Code de calcul
*
*  Copyright (C) 2023 UCL-IMMC : Vincent Legat, Miguel De Le Court
*  All rights reserved.
*
*/


// geoMeshRead("../data/mesh.txt"); // DO NOT CHANGE THIS PATHNAME
// femProblem *theProblem = femElasticityRead(theGeometry, "../data/problem.txt"); // DO NOT CHANGE THIS PATHNAME
// femSolutionWrite(nNodes, 2, theSoluce, "../data/UV.txt"); // DO NOT CHANGE THIS PATHNAME

#include "../../../fem_library/include/fem_elasticity.h"
#include "../../../fem_library/include/fem_geometry.h"

int main(void)
{
    femSolverType typeSolver = FEM_FULL; // Chosing the solver type

    femGeometry *theGeometry = geoGetGeometry();
    geoMeshRead("../../../data/mesh.txt");
    femProblem *theProblem = femElasticityRead(theGeometry, typeSolver, "../../../data/problem.txt");
    femElasticityPrint(theProblem);
    double *theSoluce = femElasticitySolve(theProblem);
    int nNodes = theGeometry->theNodes->nNodes;
    femSolutionWrite(nNodes, 2, theSoluce, "../../../data/UV.txt");
    femElasticityFree(theProblem);
    geoFree();
    return EXIT_SUCCESS;
}