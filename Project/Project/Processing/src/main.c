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

/*
* BE CAREFULL TO THESE PATHNAMES (DO NOT REMOVE IT!)
*   geoMeshRead("../data/mesh.txt");                                                // DO NOT CHANGE THIS PATHNAME
*   femProblem *theProblem = femElasticityRead(theGeometry, "../data/problem.txt"); // DO NOT CHANGE THIS PATHNAME
*   femSolutionWrite(nNodes, 2, theSoluce, "../data/UV.txt");                       // DO NOT CHANGE THIS PATHNAME
*/

#include "../../../fem_library/include/fem.h"

#include <getopt.h>

int main(int argc, char *argv[])
{
    // Deal with the options arguments
    int opt;
    int exampleUsage   = FALSE;
    while ((opt = getopt(argc, argv, "eh")) != -1)
    {
        switch (opt)
        {;
            case 'e':
                exampleUsage = TRUE;
                break;
            case 'h':
                printf("Usage: %s [-e] [-h]\n", argv[0]);
                printf("Options:\n");
                printf("  -e : Start the program with the example mesh\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-e] [-h]\n", argv[0]);
                return EXIT_FAILURE;
        }
    }

    femSolverType typeSolver = FEM_FULL; // Chosing the solver type

    femGeometry *theGeometry = geoGetGeometry();
    femProblem *theProblem;
    if (exampleUsage == TRUE)
    {
        geoMeshRead("../data/mesh_example.txt");
        theProblem = femElasticityRead(theGeometry, typeSolver, "../data/problem_example.txt");
    }
    else
    {
        geoMeshRead("../data/mesh.txt");
        theProblem = femElasticityRead(theGeometry, typeSolver, "../data/problem.txt");
    }
    
    femElasticityPrint(theProblem);

    double *theSoluce = femElasticitySolve(theProblem);
    double *theForces = femElasticityForces(theProblem); // TODO : Use theForces

    int nNodes = theGeometry->theNodes->nNodes;
    if (exampleUsage == TRUE) { femSolutionWrite(nNodes, 2, theSoluce, "../data/UV_example.txt"); }
    else                      { femSolutionWrite(nNodes, 2, theSoluce, "../data/UV.txt"); }
    femElasticityFree(theProblem);
    geoFree();
    return EXIT_SUCCESS;
}