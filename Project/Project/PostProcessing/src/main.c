/*
*  main.c
*  Projet 2023-2024
*  Elasticite lineaire plane
*
*  Postprocesseur
*
*  Copyright (C) 2023 UCL-IMMC : Vincent Legat, Miguel De Le Court
*  All rights reserved.
*
*/

#include <getopt.h>

#include "../../../fem_library/include/fem.h"
#include "../../../glfem_library/glfem.h"
#include "../../var.h"
#include "../../Processing/src/homework.c"

double constFunct(double x, double y) { return 1.0; }

// void femElasticitySigma(femProblem *theProblem, double *sigmaXX, double *sigmaYY, double *sigmaXY)
// {
//     femIntegration *theRule  = theProblem->rule;
//     femDiscrete *theSpace    = theProblem->space;
//     femGeometry *theGeometry = theProblem->geometry;
//     femNodes *theNodes       = theGeometry->theNodes;
//     femMesh *theMesh         = theGeometry->theElements;

//     double *theSoluce = theProblem->soluce;
    
//     int nLocal = theSpace->n;
//     // int *number = theMesh->nodes->number;    

//     double xLoc, x[nLocal], y[nLocal], phi[nLocal], dphidxsi[nLocal], dphideta[nLocal], dphidx[nLocal], dphidy[nLocal], u[nLocal], v[nLocal];
//     int iElem, iInteg, iEdge, i, map[nLocal], mapX[nLocal], mapY[nLocal];
    
//     // Premier systeme pour epsilon_x_x
//     femSolver *theSolver = femSolverFullCreate(theNodes->nNodes);
//     double **A = femSolverGetA(theSolver);
//     double *B  = femSolverGetB(theSolver);

//     for (iElem = 0; iElem < theMesh->nElem; iElem++)
//     {
//         for (i = 0; i < theSpace->n; i++)
//         {
//             map[i] = theMesh->elem[iElem * nLocal + i];
//             x[i] = theNodes->X[map[i]];
//             y[i] = theNodes->Y[map[i]];
//             u[i] = theSoluce[2 * map[i]];
//             // map[i] = number[map[i]];
//         }

//         for (iInteg = 0; iInteg < theRule->n; iInteg++)
//         {
//             double xsi    = theRule->xsi[iInteg];
//             double eta    = theRule->eta[iInteg];
//             double weight = theRule->weight[iInteg];

//             femDiscretePhi2(theSpace, xsi, eta, phi);
//             femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

//             double dxdxsi = 0.0; double dydxsi = 0.0;
//             double dxdeta = 0.0; double dydeta = 0.0;
//             for (i = 0; i < theSpace->n; i++)
//             {
//                 dxdxsi += x[i] * dphidxsi[i];
//                 dxdeta += x[i] * dphideta[i];
//                 dydxsi += y[i] * dphidxsi[i];
//                 dydeta += y[i] * dphideta[i];
//             }

//             double jac = dxdxsi * dydeta - dxdeta * dydxsi;
//             if (jac < 0.0) { printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n"); }
//             jac = fabs(jac);

//             for (i = 0; i < theSpace->n; i++)
//             {
//                 dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
//                 dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
//             }

//             double weightedJac = jac * weight;

//             for (int i = 0; i < theSpace->n; i++)
//             {
//                 for (int j = 0; j < theSpace->n; j++)
//                 {
//                     A[map[i]][map[i]] += phi[i] * phi[j] * weightedJac;
//                 }
//                 B[map[i]] += phi[i] * u[i] * dphidx[i] * weightedJac;
//             }
//         }
//     }

//     femSolverEliminate(theSolver);
//     double *epsilon_x_x = (double *) malloc(theSolver->size * sizeof(double));
//     if (epsilon_x_x == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }
//     memcpy(epsilon_x_x, femSolverGetB(theSolver), theSolver->size * sizeof(double));


//     // Deuxieme systeme pour epsilon_y_y

//     theSolver = femSolverFullCreate(theNodes->nNodes);

//     for (iElem = 0; iElem < theMesh->nElem; iElem++)
//     {
//         for (i = 0; i < theSpace->n; i++)
//         {
//             map[i] = theMesh->elem[iElem * nLocal + i];
//             x[i] = theNodes->X[map[i]];
//             y[i] = theNodes->Y[map[i]];
//             v[i] = theProblem->soluce[2 * map[i] + 1];
//             // map[i] = number[map[i]];
//         }

//         for (iInteg = 0; iInteg < theRule->n; iInteg++)
//         {
//             double xsi    = theRule->xsi[iInteg];
//             double eta    = theRule->eta[iInteg];
//             double weight = theRule->weight[iInteg];

//             femDiscretePhi2(theSpace, xsi, eta, phi);
//             femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

//             double dxdxsi = 0.0; double dydxsi = 0.0;
//             double dxdeta = 0.0; double dydeta = 0.0;

//             for (i = 0; i < theSpace->n; i++)
//             {
//                 dxdxsi += x[i] * dphidxsi[i];
//                 dxdeta += x[i] * dphideta[i];
//                 dydxsi += y[i] * dphidxsi[i];
//                 dydeta += y[i] * dphideta[i];
//             }

//             double jac = dxdxsi * dydeta - dxdeta * dydxsi;
//             if (jac < 0.0) { printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n"); }
//             jac = fabs(jac);

//             for (i = 0; i < theSpace->n; i++)
//             {
//                 dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
//                 dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
//             }

//             double weightedJac = jac * weight;
//             double **A = femSolverGetA(theSolver);
//             double *B = femSolverGetB(theSolver);

//             for (int i = 0; i < theSpace->n; i++)
//             {
//                 for (int j = 0; j < theSpace->n; j++)
//                 {
//                     A[map[i]][map[i]] += phi[i] * phi[j] * weightedJac;
//                 }
//                 B[map[i]] += phi[i] * dphidy[i] * v[i] * weightedJac;
//             }
//         }
//     }

//     femSolverEliminate(theSolver);
//     double *epsilon_y_y = (double *) malloc(theSolver->size * sizeof(double));
//     if (epsilon_y_y == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }
//     memcpy(epsilon_y_y, femSolverGetB(theSolver), theSolver->size * sizeof(double));

//     for (int i = 0; i < theNodes->nNodes; i++)
//     {
//         if (epsilon_y_y[i] != 0.0) { printf("%f\n", epsilon_y_y[i]); }    
//     }

//     // Troisieme systeme pour epsilon_x_y

//     theSolver = femSolverFullCreate(theNodes->nNodes);

//     for (iElem = 0; iElem < theMesh->nElem; iElem++)
//     {
//         for (i = 0; i < theSpace->n; i++)
//         {
//             map[i] = theMesh->elem[iElem * nLocal + i];
//             x[i] = theNodes->X[map[i]];
//             y[i] = theNodes->Y[map[i]];
//             u[i] = theProblem->soluce[2 * map[i]];
//             v[i] = theProblem->soluce[2 * map[i] + 1];
//             // map[i] = number[map[i]];
//         }

//         for (iInteg = 0; iInteg < theRule->n; iInteg++)
//         {
//             double xsi    = theRule->xsi[iInteg];
//             double eta    = theRule->eta[iInteg];
//             double weight = theRule->weight[iInteg];

//             femDiscretePhi2(theSpace, xsi, eta, phi);
//             femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

//             double dxdxsi = 0.0; double dydxsi = 0.0;
//             double dxdeta = 0.0; double dydeta = 0.0;

//             for (i = 0; i < theSpace->n; i++)
//             {
//                 dxdxsi += x[i] * dphidxsi[i];
//                 dxdeta += x[i] * dphideta[i];
//                 dydxsi += y[i] * dphidxsi[i];
//                 dydeta += y[i] * dphideta[i];
//             }

//             double jac = dxdxsi * dydeta - dxdeta * dydxsi;
//             if (jac < 0.0) { printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n"); }
//             jac = fabs(jac);

//             for (i = 0; i < theSpace->n; i++)
//             {
//                 dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
//                 dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
//             }

//             double weightedJac = jac * weight;
//             double **A = femSolverGetA(theSolver);
//             double *B = femSolverGetB(theSolver);

//             for (int i = 0; i < theSpace->n; i++)
//             {
//                 for (int j = 0; j < theSpace->n; j++)
//                 {
//                     A[map[i]][map[i]] += phi[i] * phi[j] * weightedJac;
//                 }
//                 B[map[i]] += phi[i] * (dphidy[i] * u[i] + dphidx[i] * v[i]) * weightedJac / 2.0;
//             }
//         }
//     }

//     femSolverEliminate(theSolver);
    
//     double *epsilon_x_y = (double *) malloc(theSolver->size * sizeof(double));
//     if (epsilon_x_y == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }
//     memcpy(epsilon_x_y, femSolverGetB(theSolver), theSolver->size * sizeof(double));

//     for (int i = 0; i < theNodes->nNodes; i++)
//     {
//         if (epsilon_x_y[i] != 0.0) { printf("%f\n", epsilon_x_y[i]); }    
//     }

//     double a = theProblem->A;
//     double b = theProblem->B;
//     double c = theProblem->C;

//     for (int i = 0; i < theSolver->size; i++)
//     {
//         sigmaXX[i] = a * epsilon_x_x[i] + b * epsilon_y_y[i];
//         sigmaXY[i] = 2 * c * epsilon_x_y[i];
//         sigmaYY[i] = b * epsilon_x_x[i] + a * epsilon_y_y[i];
//     }
// }

double *femElasticityVonMises(femProblem *theProblem, double *sigmaXX, double *sigmaYY, double *sigmaXY, int nNodes)
{
    double *vonMises = (double *) malloc(theProblem->geometry->theNodes->nNodes * sizeof(double));
    if (vonMises == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return NULL; }

    for (int i = 0; i < nNodes; i++)
    {
        vonMises[i] = sqrt(sigmaXX[i] * sigmaXX[i] - sigmaXX[i] * sigmaYY[i] + sigmaYY[i] * sigmaYY[i] + 3 * sigmaXY[i] * sigmaXY[i]);
    }
    return vonMises;
}

void femElasticitySigma(femProblem *theProblem, double *sigmaXX, double *sigmaYY, double *sigmaXY)
{
    femIntegration *theRule  = theProblem->rule;
    femDiscrete *theSpace    = theProblem->space;
    femGeometry *theGeometry = theProblem->geometry;
    femNodes *theNodes       = theGeometry->theNodes;
    femMesh *theMesh         = theGeometry->theElements;

    int nLocal = theSpace->n;

    femSolver *theSolver;
    double **M1, **M2, **M3, *B1, *B2, *B3, *epsilonXX, *epsilonYY, *epsilonXY, *theSoluce;
    double x[nLocal], y[nLocal], phi[nLocal], dphidxsi[nLocal], dphideta[nLocal], dphidx[nLocal], dphidy[nLocal];
    double u[nLocal], v[nLocal], xsi, eta, weight, jac, weightedJac, dxdxsi, dxdeta, dydxsi, dydeta, a, b, c;
    int iElem, iInteg, iEdge, i, j, map[nLocal], nNodes;

    nNodes    = theNodes->nNodes;
    theSoluce = theProblem->soluce;    
    theSolver = femSolverFullCreate(nNodes);

    M1 = femSolverGetA(theSolver);
    B1 = femSolverGetB(theSolver);

    B2 = (double *) malloc(nNodes * sizeof(double));
    B3 = (double *) malloc(nNodes * sizeof(double));

    if (B2 == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }
    if (B3 == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }

    for (i = 0; i < nNodes; i++) { B2[i] = 0.0; B3[i] = 0.0; }
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++)
    {
        for (i = 0; i < theSpace->n; i++)
        {
            map[i] = theMesh->elem[iElem * nLocal + i];
            x[i] = theNodes->X[map[i]];
            y[i] = theNodes->Y[map[i]];
            u[i] = theSoluce[2 * map[i]];
            v[i] = theSoluce[2 * map[i] + 1];
        }

        for (iInteg = 0; iInteg < theRule->n; iInteg++)
        {
            xsi    = theRule->xsi[iInteg];
            eta    = theRule->eta[iInteg];
            weight = theRule->weight[iInteg];

            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            dxdxsi = 0.0; dydxsi = 0.0;
            dxdeta = 0.0; dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++)
            {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }

            jac = dxdxsi * dydeta - dxdeta * dydxsi;
            if (jac < 0.0) { printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n"); }
            jac = fabs(jac);

            for (i = 0; i < theSpace->n; i++)
            {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            weightedJac = jac * weight;

            for (i = 0; i < theSpace->n; i++)
            {
                for (j = 0; j < theSpace->n; j++) { M1[map[i]][map[j]] += phi[i] * phi[j] * weightedJac; }
                B1[map[i]] += phi[i] * dphidx[i] * u[i] * weightedJac;
                B2[map[i]] += phi[i] * dphidy[i] * v[i] * weightedJac;
                B3[map[i]] += phi[i] * (dphidy[i] * u[i] + dphidx[i] * v[i]) * weightedJac / 2.0;
            }
        }
    }

    M2 = (double **) malloc(sizeof(double *) * nNodes);
    M3 = (double **) malloc(sizeof(double *) * nNodes);

    if (M2 == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    if (M3 == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

    for (i = 0; i < nNodes; i++)
    {
        M2[i] = (double *) malloc(nNodes * sizeof(double));
        M3[i] = (double *) malloc(nNodes * sizeof(double));

        if (M2[i] == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        if (M3[i] == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

        memcpy(M2[i], M1[i], nNodes * sizeof(double));
        memcpy(M3[i], M1[i], nNodes * sizeof(double));
    }

    epsilonXX = (double *) malloc(nNodes * sizeof(double));
    epsilonYY = (double *) malloc(nNodes * sizeof(double));
    epsilonXY = (double *) malloc(nNodes * sizeof(double));

    if (epsilonXX == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }
    if (epsilonYY == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }
    if (epsilonXY == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }

    femSolverEliminate(theSolver);
    memcpy(epsilonXX, femSolverGetB(theSolver), nNodes * sizeof(double));

    femSolverSet(theSolver, M2, B2);
    femSolverEliminate(theSolver);
    memcpy(epsilonYY, femSolverGetB(theSolver), nNodes * sizeof(double));

    femSolverSet(theSolver, M3, B3);
    femSolverEliminate(theSolver);
    memcpy(epsilonXY, femSolverGetB(theSolver), nNodes * sizeof(double));

    a = theProblem->A;
    b = theProblem->B;
    c = theProblem->C;

    for (i = 0; i < nNodes; i++)
    {
        sigmaXX[i] = a * epsilonXX[i] + b * epsilonYY[i];
        sigmaXY[i] = 2 * c * epsilonXY[i];
        sigmaYY[i] = b * epsilonXX[i] + a * epsilonYY[i];
    }

    free(epsilonXX); epsilonXX = NULL;
    free(epsilonYY); epsilonYY = NULL;
    free(epsilonXY); epsilonXY = NULL;
    for (i = 0; i < nNodes; i++) { free(M2[i]); M2[i] = NULL; free(M3[i]); M3[i] = NULL; }
    free(M2); M2 = NULL;
    free(M3); M3 = NULL;
    free(B2); B2 = NULL;
    free(B3); B3 = NULL;
}

double *femElasticityForces(femProblem *theProblem)
{
    femNodes *theNodes   = theProblem->geometry->theNodes;
    double *theSoluce    = theProblem->soluce;
    femSolver *theSolver = theProblem->solver;
    int size = theSolver->size;

    double *soluce       = (double *) malloc(size * sizeof(double));
    double *residuals    = (double *) malloc(size * sizeof(double));
    int *inverted_number = (int *) malloc(theNodes->nNodes * sizeof(int));
    if (soluce == NULL)          { Error("Allocation Error\n"); exit(EXIT_FAILURE); return NULL; }
    if (residuals == NULL)       { Error("Allocation Error\n"); exit(EXIT_FAILURE); return NULL; }
    if (inverted_number == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return NULL; }

    for (int i = 0; i < size; i++) { residuals[i] = 0.0; }

    // Numerote the nodes of the solution
    for (int i = 0; i < theNodes->nNodes; i++)
    {
        soluce[2 * theNodes->number[i]]     = theSoluce[2 * i];
        soluce[2 * theNodes->number[i] + 1] = theSoluce[2 * i + 1];
    }

    femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem, 1.0);

    double **A = femSolverGetA(theSolver);
    double *B = femSolverGetB(theSolver);

    femSolverGetResidual(theSolver, residuals, soluce);

    // Remove the numerotation of the nodes of the residuals
    for (int i = 0; i < theNodes->nNodes; i++) { inverted_number[theNodes->number[i]] = i; }

    for (int i = 0; i < theNodes->nNodes; i++)
    {
        theProblem->residuals[2 * inverted_number[i]]     = residuals[2 * i];
        theProblem->residuals[2 * inverted_number[i] + 1] = residuals[2 * i + 1];
    }

    // Apply the Dirichlet boundary conditions to the system to plot the final system
    femElasticityApplyDirichlet(theProblem, 1.0);

    free(residuals); residuals = NULL;
    free(soluce); soluce = NULL;
    free(inverted_number); inverted_number = NULL;
    return theProblem->residuals;
}

int main(int argc, char *argv[])
{
    // Deal with the options arguments
    int opt;
    int resultVisualizer = TRUE;
    int exampleUsage     = FALSE;
    int animation        = FALSE;
    int plot             = TRUE;
    while ((opt = getopt(argc, argv, "areph")) != -1)
    {
        switch (opt)
        {
            case 'a':
                animation = TRUE;
            case 'r':
                resultVisualizer = FALSE;
                break;
            case 'e':
                exampleUsage = TRUE;
                break;
            case 'p':
                plot = FALSE;
                break;
            case 'h':
                printf("Usage: %s [-r] [-e] [-h]\n", argv[0]);
                printf("Options:\n");
                printf("  -r : Disable the result visualizer\n");
                printf("  -e : Start the program with the example mesh\n");
                printf("  -p : Disable the plot\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-r] [-e] [-h]\n", argv[0]);
                return EXIT_FAILURE;
        }
    }
    
    printf("\n\n");
    printf("    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    N : Next domain highlighted\n");
    printf("    S : Display the matrix \n");
    printf("    X : Display the horizontal forces \n");
    printf("    Y : Display the vertical forces \n");
    printf("\n\n");

    /***************************/
    /* 1 : Lecture des donnees */ 
    /***************************/

    femGeometry *theGeometry = geoGetGeometry();

    int n;
    double *theSoluce;
    femProblem *theProblem;

    if (exampleUsage == TRUE)
    {
        geoMeshRead("../../Processing/data/mesh_example.txt", discretType);
        theProblem = femElasticityRead(theGeometry, typeSolver, "../../Processing/data/problem_example.txt", renumType, discretType, TRUE);
        theSoluce = theProblem->soluce;
        n = theGeometry->theNodes->nNodes;
        femSolutiondRead(2 * n, theSoluce, "../../Processing/data/UV_example.txt");
        femElasticityPrint(theProblem);
    }
    else
    {
        geoMeshRead("../../Processing/data/mesh.txt", discretType);
        theProblem = femElasticityRead(theGeometry, typeSolver, "../../Processing/data/problem.txt", renumType, discretType, TRUE);
        theSoluce = theProblem->soluce;
        n = theGeometry->theNodes->nNodes;
        femSolutiondRead(2 * n, theSoluce, "../../Processing/data/UV.txt");
        femElasticityPrint(theProblem);
    }

    double rho = theProblem->rho;
    double gy = theProblem->gy;

    /*********************************/
    /* 2 : Calcul des forces         */
    /*     Création du système final */
    /*********************************/
    
    double *theForces = femElasticityForces(theProblem);
    double area       = femElasticityIntegrate(theProblem, constFunct);

    femNodes *theNodes = theGeometry->theNodes;
    int nNodes = theNodes->nNodes;

    double *sigmaXX = (double *) malloc(theGeometry->theNodes->nNodes * sizeof(double));
    if (sigmaXX == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    double *sigmaYY = (double *) malloc(theGeometry->theNodes->nNodes * sizeof(double));
    if (sigmaYY == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    double *sigmaXY = (double *) malloc(theGeometry->theNodes->nNodes * sizeof(double));
    if (sigmaXY == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    femElasticitySigma(theProblem, sigmaXX, sigmaYY, sigmaXY);

    double *vonMises = femElasticityVonMises(theProblem, sigmaXX, sigmaYY, sigmaXY, nNodes);
    
    FILE *file = fopen("../../Processing/data/vonMises.txt", "w");
    if (file == NULL) { Error("File opening error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    for (int i = 0; i < nNodes; i++) { fprintf(file, "%f\n", vonMises[i]); }
    fclose(file);
    
    /****************************************************/
    /* 3 : Deformation du maillage pour le plot final   */ 
    /*     Creation du champ de la norme du deplacement */
    /****************************************************/

    double deformationFactor;
    if (exampleUsage == TRUE) { deformationFactor = 1e5; } // To change the deformation factor
    else                      { deformationFactor = 1e4; } // To change the deformation factor

    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    if (normDisplacement == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    double *forcesX = malloc(theNodes->nNodes * sizeof(double));
    if (forcesX == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    double *forcesY = malloc(theNodes->nNodes * sizeof(double));
    if (forcesY == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }

    for (int i = 0; i < n; i++)
    {
        theNodes->X[i] += theSoluce[2 * i + 0] * deformationFactor;
        theNodes->Y[i] += theSoluce[2 * i + 1] * deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2 * i + 0] * theSoluce[2 * i + 0] + theSoluce[2 * i + 1] * theSoluce[2 * i + 1]);
        forcesX[i] = theForces[2 * i + 0];
        forcesY[i] = theForces[2 * i + 1];
    }

    double hMin = femMin(normDisplacement, n);
    double hMax = femMax(normDisplacement, n);

    printf(" ==== Minimum displacement          : %14.7e \n", hMin);
    printf(" ==== Maximum displacement          : %14.7e \n", hMax);

    /*********************************************/
    /* 4 : Calcul de la force globale resultante */ 
    /*********************************************/

    double theGlobalForce[2] = {0, 0};
    for (int i = 0; i < theProblem->geometry->theNodes->nNodes; i++)
    {
        theGlobalForce[0] += theForces[2 * i];
        theGlobalForce[1] += theForces[2 * i + 1];
    }

    printf(" ==== Global horizontal force       : %14.7e [N] \n",theGlobalForce[0]);
    printf(" ==== Global vertical force         : %14.7e [N] \n",theGlobalForce[1]);
    printf(" ==== Weight                        : %14.7e [N] \n", area * rho * gy);

    /******************************/
    /* 5 : Lance le script python */ 
    /******************************/

    if (plot == TRUE)
    {
        char command[100] = "python3 ../src/plot.py";
        if (exampleUsage == TRUE) { strcat(command, " -e"); }
        if (animation == TRUE)    { strcat(command, " -a"); }
        printf(" ==== Command : %s \n", command);
        system(command);
    }

    /*********************/
    /* 6 : Visualisation */ 
    /*********************/

    if (!resultVisualizer) { return EXIT_SUCCESS; }

    int mode = 1;
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME + 5];

    GLFWwindow *window = glfemInit("EPL1110 : Project 2022-23 ");
    glfwMakeContextCurrent(window);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    do {
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glfemReshapeWindows(window, theGeometry->theNodes, w, h);

        t = glfwGetTime();
        if (glfwGetKey(window, 'D') == GLFW_PRESS) { mode = 0; }
        if (glfwGetKey(window, 'V') == GLFW_PRESS) { mode = 1; }
        if (glfwGetKey(window, 'S') == GLFW_PRESS) { mode = 2; }
        if (glfwGetKey(window, 'X') == GLFW_PRESS) { mode = 3; }
        if (glfwGetKey(window, 'Y') == GLFW_PRESS) { mode = 4; }
        if (glfwGetKey(window, 'A') == GLFW_PRESS) { mode = 5; }
        if (glfwGetKey(window, 'B') == GLFW_PRESS) { mode = 6; }
        if (glfwGetKey(window, 'C') == GLFW_PRESS) { mode = 7; }
        if (glfwGetKey(window, 'U') == GLFW_PRESS) { mode = 8; }
        if (glfwGetKey(window, 'N') == GLFW_PRESS && freezingButton == FALSE)
        {
            domain++;
            freezingButton = TRUE;
            told = t;
        }

        if (t - told > 0.5) { freezingButton = FALSE; }
        if (mode == 1)
        {
            glfemPlotField(theGeometry->theElements, normDisplacement);
            glfemPlotMesh(theGeometry->theElements);
            sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
            glColor3f(1.0, 0.0, 0.0);
            glfemMessage(theMessage);
        }
        else if (mode == 0)
        {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain(theGeometry->theDomains[domain]);
            sprintf(theMessage, "%s : %d ", theGeometry->theDomains[domain]->name, domain);
            glColor3f(1.0, 0.0, 0.0);
            glfemMessage(theMessage);
        } 
        else if (mode == 2)
        {
            glColor3f(1.0, 0.0, 0.0);
            glfemPlotSolver(theProblem->solver, theProblem->solver->size, w, h);
        }
        else if (mode == 3)
        {
            glfemPlotField(theGeometry->theElements, forcesX);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if (mode == 4)
        {
            glfemPlotField(theGeometry->theElements, forcesY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if (mode == 5)
        {
            glfemPlotField(theGeometry->theElements, sigmaXX);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if (mode == 6)
        {
            glfemPlotField(theGeometry->theElements, sigmaYY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if (mode == 7)
        {
            glfemPlotField(theGeometry->theElements, sigmaXY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if (mode == 8)
        {
            printf("%.12f\n", femMax(vonMises, nNodes));
            glfemPlotField(theGeometry->theElements, vonMises);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else { printf("Mode %d not implemented\n", mode); }  

        glfwSwapBuffers(window);
        glfwPollEvents();
    } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) != 1);
    // Check if the ESC key was pressed or the window was closed

    free(normDisplacement);
    femElasticityFree(theProblem);
    geoFree();
    glfwTerminate();

    exit(EXIT_SUCCESS);
    return EXIT_SUCCESS;
}