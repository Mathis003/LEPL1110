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

#include <getopt.h>
#include <unistd.h>
#include <time.h>

femElementType elementType  = FEM_TRIANGLE; // FEM_QUAD or FEM_TRIANGLE
femDiscreteType discretType = FEM_DISCRETE_TYPE_LINEAR; // FEM_DISCRETE_TYPE_LINEAR or FEM_DISCRETE_TYPE_QUADRATIC
femSolverType typeSolver    = FEM_BAND;  // FEM_FULL or FEM_BAND
femRenumType renumType      = FEM_RCMK;  // FEM_NO or FEM_XNUM or FEM_YNUM or FEM_RCMK


/* FINAL PATH */

/* "../data/mesh.txt" */
/* "../data/problem.txt" */
/* ""../data/UV.txt" */

int main(int argc, char *argv[])
{
    // Deal with the options arguments
    int opt;
    int bridgeSimplified = FALSE;
    int exampleUForm_Usage = FALSE;
    int exampleBeam_Usage = FALSE;
    int animation     = FALSE;
    int showRunTime = FALSE;
    while ((opt = getopt(argc, argv, "subtah")) != -1)
    {
        switch (opt)
        {
            case 's':
                bridgeSimplified = TRUE;
                break;
            case 'u':
                exampleUForm_Usage = TRUE;
                break;
            case 'b':
                exampleBeam_Usage = TRUE;
                break;
            case 't':
                showRunTime = TRUE;
                break;
            case 'a':
                animation = TRUE;
                break;
            case 'h':
                printf("Usage: %s [-s] [-u] [-b] [-t] [-a] [-h]\n", argv[0]);
                printf("Options:\n");
                printf("  -s : Start the program with the bridge without stay cables and pylon\n");
                printf("  -u : Start the program with the U mesh\n");
                printf("  -b : Start the program with the beam mesh\n");
                printf("  -t : Time the program execution\n");
                printf("  -a : Launch the program 50 times to create an animation in ProjectPostProcessor\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-s] [-u] [-b] [-t] [-a] [-h]\n", argv[0]);
                return EXIT_FAILURE;
        }
    }

    femGeometry *theGeometry = geoGetGeometry();
    
    femProblem *theProblem;
    if (exampleUForm_Usage == TRUE)
    {
        geoMeshRead("../../Rapport/data/mesh_example.txt", discretType);
        theProblem = femElasticityRead(theGeometry, typeSolver, "../../Rapport/data/problem_example.txt", renumType, discretType, TRUE);
    }
    else if (exampleBeam_Usage == TRUE)
    {
        geoMeshRead("../../Rapport/data/mesh_beam.txt", discretType);
        theProblem = femElasticityRead(theGeometry, typeSolver, "../../Rapport/data/problem_beam.txt", renumType, discretType, TRUE);
    }
    else if (bridgeSimplified == TRUE)
    {
        geoMeshRead("../../Rapport/data/mesh_simplified.txt", discretType);
        theProblem = femElasticityRead(theGeometry, typeSolver, "../../Rapport/data/problem_simplified.txt", renumType, discretType, TRUE);
    }
    else
    {
        geoMeshRead("../data/mesh.txt", discretType);
        theProblem = femElasticityRead(theGeometry, typeSolver, "../data/problem.txt", renumType, discretType, TRUE);
    }
    
    femElasticityPrint(theProblem);

    if (!animation)
    {
        clock_t start;
        if (showRunTime == TRUE)
        { 
            printf("Solving the system...\n"); 
            start = clock();    
        }
        double *theSoluce = femElasticitySolve(theProblem, renumType, 1.0);
        int nNodes = theGeometry->theNodes->nNodes;

        if( showRunTime == TRUE) { printf("Run time: %f seconds\n", (double)(clock() - start) / CLOCKS_PER_SEC); }

        if (exampleUForm_Usage == TRUE)     { femSolutionWrite(nNodes, 2, theSoluce, "../../Rapport/data/UV_example.txt"); }
        else if (exampleBeam_Usage == TRUE) { femSolutionWrite(nNodes, 2, theSoluce, "../../Rapport/data/UV_beam.txt"); }
        else if (bridgeSimplified == TRUE)  { femSolutionWrite(nNodes, 2, theSoluce, "../../Rapport/data/UV_simplified.txt"); }
        else                                { femSolutionWrite(nNodes, 2, theSoluce, "../data/UV.txt"); }
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
        if (access("../data/animations", F_OK) == -1) { system("mkdir -p ../../Rapport/data/animations"); }

        for (int i = 1; i <= NB_IMAGES_ANIMATION; i++)
        {
            theSoluce = femElasticitySolve(theProblem, renumType, FACTOR * i);
            nNodes = theGeometry->theNodes->nNodes;
            if (exampleUForm_Usage == TRUE)     { sprintf(filename, "../../Rapport/data/animations/UV_example_%d.txt", i); }
            else if (exampleBeam_Usage == TRUE) { sprintf(filename, "../../Rapport/data/animations/UV_beam_%d.txt", i); }
            else if (bridgeSimplified == TRUE)  { sprintf(filename, "../../Rapport/data/animations/UV_simplified_%d.txt", i); }
            else                                { sprintf(filename, "../../Rapport/data/animations/UV_%d.txt", i); }         
            femSolutionWrite(nNodes, 2, theSoluce, filename);
        }

        femElasticityFree(theProblem);
        geoFree();
        return EXIT_SUCCESS;
    }
}