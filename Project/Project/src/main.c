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

#include "fem.h"

int main(void)
{
    femGeo *theGeometry = geoGetGeometry();
    geoMeshRead("../data/mesh.txt"); // DO NOT CHANGE THIS PATHNAME
    femProblem *theProblem = femElasticityRead(theGeometry, "../data/problem.txt"); // DO NOT CHANGE THIS PATHNAME
    femElasticityPrint(theProblem);
    double *theSoluce = femElasticitySolve(theProblem);
    int nNodes = theGeometry->theNodes->nNodes;

    // Write the solution to a file in this directory and in the ProjectPostProcess directory
    femSolutionWrite(nNodes, 2, theSoluce, "../data/UV.txt"); // DO NOT CHANGE THIS PATHNAME
    femSolutionWrite(nNodes, 2, theSoluce, "../../ProjectPostProcess/data/UV.txt"); // TO REMOVE IN THE FINAL VERSION (FOR SUBMISSION)

    femElasticityFree(theProblem);
    geoFree();
    return EXIT_SUCCESS;
}