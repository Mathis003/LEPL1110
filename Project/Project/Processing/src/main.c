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
#include "../../var.h"

#include <getopt.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
    // Deal with the options arguments
    int opt;
    int exampleUsage = FALSE;
    int animaion     = FALSE;
    while ((opt = getopt(argc, argv, "eah")) != -1)
    {
        switch (opt)
        {
            case 'e':
                exampleUsage = TRUE;
                break;
            case 'a':
                animaion = TRUE;
                break;
            case 'h':
                printf("Usage: %s [-e] [-h]\n", argv[0]);
                printf("Options:\n");
                printf("  -a : Launch the program 50 times to create an animation in Post-Processor\n");
                printf("  -e : Start the program with the example mesh\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-e] [-h]\n", argv[0]);
                return EXIT_FAILURE;
        }
    }

    femGeometry *theGeometry = geoGetGeometry();
    
    femProblem *theProblem;
    if (exampleUsage == TRUE)
    {
        geoMeshRead("../data/mesh_example.txt", discretType);
        theProblem = femElasticityRead(theGeometry, typeSolver, "../data/problem_example.txt", renumType, discretType, TRUE);
    }
    else
    {
        geoMeshRead("../data/mesh.txt", discretType);
        theProblem = femElasticityRead(theGeometry, typeSolver, "../data/problem.txt", renumType, discretType, TRUE);
    }
    
    femElasticityPrint(theProblem);

    if (!animaion)
    {
        double *theSoluce = femElasticitySolve(theProblem, renumType, 1.0);

        int nNodes = theGeometry->theNodes->nNodes;
        if (exampleUsage == TRUE) { femSolutionWrite(nNodes, 2, theSoluce, "../data/UV_example.txt"); }
        else                      { femSolutionWrite(nNodes, 2, theSoluce, "../data/UV.txt"); }
        femElasticityFree(theProblem);
        geoFree();
        return EXIT_SUCCESS;
    }
    else
    {
        const int NB_IMAGES_ANIMATION = 50;
        double FACTOR = 1.0 / NB_IMAGES_ANIMATION;
        double *theSoluce;
        int nNodes;
        char filename[100];
        if (access("../data/animations", F_OK) == -1) { system("mkdir -p ../data/animations"); }

        for (int i = 1; i <= NB_IMAGES_ANIMATION; i++)
        {
            theSoluce = femElasticitySolve(theProblem, renumType, FACTOR * i);
            nNodes = theGeometry->theNodes->nNodes;
            if (exampleUsage == TRUE) { sprintf(filename, "../data/animations/UV_example_%d.txt", i); }
            else                      { sprintf(filename, "../data/animations/UV_%d.txt", i); }         
            femSolutionWrite(nNodes, 2, theSoluce, filename);
        }

        femElasticityFree(theProblem);
        geoFree();
        return EXIT_SUCCESS;
    }
}