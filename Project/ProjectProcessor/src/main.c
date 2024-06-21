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
femDiscreteType discretType = FEM_DISCRETE_TYPE_LINEAR; // FEM_DISCRETE_TYPE_LINEAR or FEM_DISCRETE_TYPE_QUADRATIC (FEM_DISCRETE_TYPE_QUADRATIC not implemented yet)
femSolverType typeSolver    = FEM_BAND;  // FEM_FULL or FEM_BAND
femRenumType renumType      = FEM_RCMK;  // FEM_NO or FEM_XNUM or FEM_YNUM or FEM_RCMK

int main(int argc, char *argv[])
{
    int opt;
    int bridgeSimplified_Usage = FALSE;
    int exampleUForm_Usage     = FALSE;
    int exampleBeam_Usage      = FALSE;
    int animation              = FALSE;
    int animationPosition      = FALSE;
    int showRunTime            = FALSE;
    while ((opt = getopt(argc, argv, "subtaxh")) != -1)
    {
        switch (opt)
        {
            case 's':
                bridgeSimplified_Usage = TRUE;
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
            case 'x':
                animationPosition = TRUE;
                animation         = TRUE;
                break;
            case 'h':
                printf("Usage: %s [-s] [-u] [-b] [-t] [-a] [-x] [-h]\n", argv[0]);
                printf("Options:\n");
                printf("  -s : Start the program with the bridge without stay cables and pylon\n");
                printf("  -u : Start the program with the U example mesh\n");
                printf("  -b : Start the program with the beam mesh\n");
                printf("  -t : Time the program execution\n");
                printf("  -a : Launch the program 50 times to create an animation of bridge's deformation (in ProjectPostProcessor)\n");
                printf("  -x : Launch the program 50 times to create an animation of car's movement (in ProjectPostProcessor)\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-s] [-u] [-b] [-t] [-a] [-x] [-h]\n", argv[0]);
                return EXIT_FAILURE;
        }
    }

    femGeometry *theGeometry = geoGetGeometry();
    femProblem *theProblem;

    char *meshPath;
    char *problemPath;
    char *solutionPath;

    if (animationPosition) { bridgeSimplified_Usage = TRUE; }

    if      (bridgeSimplified_Usage) { meshPath = "../../Rapport/data/mesh_simplified.txt"; problemPath = "../../Rapport/data/problem_simplified.txt";  solutionPath = "../../Rapport/data/UV_simplified.txt"; }
    else if (exampleBeam_Usage)      { meshPath = "../../Rapport/data/mesh_beam.txt";       problemPath = "../../Rapport/data/problem_beam.txt";        solutionPath = "../../Rapport/data/UV_beam.txt"; }
    else if (exampleUForm_Usage)     { meshPath = "../../Rapport/data/mesh_example.txt";    problemPath = "../../Rapport/data/problem_example.txt";     solutionPath = "../../Rapport/data/UV_example.txt"; }
    else                             { meshPath = "../../Rapport/data/mesh.txt"; problemPath = "../../Rapport/data/problem.txt"; solutionPath = "../../Rapport/data/UV.txt"; }
    
    geoMeshRead(meshPath, discretType);
    femMeshRenumber(theGeometry->theElements, renumType); // Renumbering the mesh nodes
    theProblem = femElasticityRead(theGeometry, typeSolver, problemPath, renumType, discretType);
    
    // femElasticityPrint(theProblem);

    int nNodes = theGeometry->theNodes->nNodes;
    if (!animation)
    {
        clock_t start;
        if (showRunTime) { printf("Solving the system...\n"); start = clock(); }

        double *theSoluce = femElasticitySolve(theProblem, renumType, 1.0, -1);

        if( showRunTime) { printf("Run time: %f seconds\n", (double)(clock() - start) / CLOCKS_PER_SEC); }

        femSolutionWrite(nNodes, 2, theSoluce, solutionPath);
    }
    else
    {
        int NB_IMAGES_ANIMATION = 50;
        double FACTOR = 1.0 / NB_IMAGES_ANIMATION; // Factor to increase the load
        
        double *theSoluce;
        char filename[100];

        // Create the directory if it doesn't exist
        if (access("../data/animations", F_OK) == -1) { system("mkdir -p ../../Rapport/data/animations"); }

        for (int i = 1; i <= NB_IMAGES_ANIMATION; i++)
        {
            femSolverInit(theProblem->solver);
            
            if (animationPosition) { theSoluce = femElasticitySolve(theProblem, renumType, 1.0, i); }
            else                   { theSoluce = femElasticitySolve(theProblem, renumType, FACTOR * i, -1); }

            if (animationPosition) { sprintf(filename, "../../Rapport/data/animations/UV_animation_position_%d.txt", i); }
            else
            {
                if      (exampleUForm_Usage)     { sprintf(filename, "../../Rapport/data/animations/UV_example_%d.txt", i); }
                else if (exampleBeam_Usage)      { sprintf(filename, "../../Rapport/data/animations/UV_beam_%d.txt", i); }
                else if (bridgeSimplified_Usage) { sprintf(filename, "../../Rapport/data/animations/UV_simplified_%d.txt", i); }
                else                             { sprintf(filename, "../../Rapport/data/animations/UV_%d.txt", i); }
            }

            theSoluce = femElasticitySolve(theProblem, renumType, 1.0, i);  
            femSolutionWrite(nNodes, 2, theSoluce, filename);
            printf("Solution %d created\n", i);
        }
    }
    femElasticityFree(theProblem); theProblem = NULL;
    geoFree();
    return EXIT_SUCCESS;
}