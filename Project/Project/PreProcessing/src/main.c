/*
*  main.c
*  Projet 2023-2024
*  Elasticite lineaire plane
*
*  Pre-processeur
*
*  Copyright (C) 2023 UCL-IMMC : Vincent Legat, Miguel De Le Court
*  All rights reserved.
*
*/

#include "../../../glfem_library/glfem.h"
#include "../../../fem_library/include/fem_gmsh.h"

#include <getopt.h>

int main(int argc, char *argv[])
{
    // Deal with the options arguments
    int opt;
    int meshVisualizer = TRUE;
    while ((opt = getopt(argc, argv, "mh")) != -1)
    {
        switch (opt)
        {
            case 'm':
                meshVisualizer = FALSE;
                break;
            case 'h':
                printf("Usage: %s [-m] [-h]\n", argv[0]);
                printf("Options:\n");
                printf("  -m : Disable the mesh visualizer\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-m] [-h]\n", argv[0]);
                return EXIT_FAILURE;
        }
    }

    /************************************/
    /* 1 : Construction de la geometrie */ 
    /************************************/

    geoInitialize();
    femGeometry *theGeometry = geoGetGeometry();

    geoMeshGenerate();
    geoMeshImport();

    // Define the domain's name
    geoSetDomainName(5, "PILAR R 1");
    geoSetDomainName(6, "PILAR D 1");
    geoSetDomainName(7, "PILAR L 1");
    geoSetDomainName(1, "PILAR R 2");
    geoSetDomainName(2, "PILAR D 2");
    geoSetDomainName(3, "PILAR L 2");

    geoSetDomainName(10, "SUB ROADWAY U 1");
    geoSetDomainName(33, "SUB ROADWAY U 2");
    geoSetDomainName(38, "SUB ROADWAY U 3");
    geoSetDomainName(42, "SUB ROADWAY U 4");
    geoSetDomainName(51, "SUB ROADWAY U 5");
    geoSetDomainName(49, "SUB ROADWAY U 6");
    geoSetDomainName(56, "SUB ROADWAY U 7");
    geoSetDomainName(27, "SUB ROADWAY U 8");
    geoSetDomainName(8, "SUB ROADWAY D 1");
    geoSetDomainName(4, "SUB ROADWAY D 2");
    geoSetDomainName(0, "SUB ROADWAY D 3");
    geoSetDomainName(9, "SUB ROADWAY L");
    geoSetDomainName(28, "SUB ROADWAY R");

    geoSetDomainName(13, "ROADWAY L");
    geoSetDomainName(24, "ROADWAY R");
    // geoSetDomainName(13, "ROADWAY U 1"); // NUMERO
    // geoSetDomainName(24, "ROADWAY U 2"); // NUMERO

    geoSetDomainName(29, "WINDOW D 1");
    geoSetDomainName(30, "WINDOW L 1");
    geoSetDomainName(31, "WINDOW R 1");
    geoSetDomainName(32, "WINDOW U 1");
    geoSetDomainName(43, "WINDOW D 2");
    geoSetDomainName(44, "WINDOW L 2");
    geoSetDomainName(45, "WINDOW R 2");
    geoSetDomainName(46, "WINDOW U 2");

    geoSetDomainName(11, "PILE L 1");
    geoSetDomainName(35, "PILE R 1");
    geoSetDomainName(37, "PILE L 2");
    geoSetDomainName(41, "PILE R 2");
    geoSetDomainName(40, "PILE L 3");
    geoSetDomainName(52, "PILE R 3");
    geoSetDomainName(48, "PILE L 4");
    geoSetDomainName(55, "PILE R 4");
    geoSetDomainName(54, "PILE L 5");
    geoSetDomainName(26, "PILE R 5");

    geoSetDomainName(12, "ARC 1");
    geoSetDomainName(34, "ARC 2");
    geoSetDomainName(36, "ARC 3");
    geoSetDomainName(39, "ARC 4");
    geoSetDomainName(50, "ARC 5");
    geoSetDomainName(47, "ARC 6");
    geoSetDomainName(53, "ARC 7");
    geoSetDomainName(25, "ARC 8");

    // geoSetDomainName(0, "PYLON 1 L"); // NUMERO
    // geoSetDomainName(0, "PYLON 1 R"); // NUMERO
    // geoSetDomainName(0, "PYLON 1 UL"); // NUMERO
    // geoSetDomainName(0, "PYLON 1 UR"); // NUMERO
    // geoSetDomainName(0, "PYLON 2 L"); // NUMERO
    // geoSetDomainName(0, "PYLON 2 R"); // NUMERO
    // geoSetDomainName(0, "PYLON 2 UL"); // NUMERO
    // geoSetDomainName(0, "PYLON 2 UR"); // NUMERO
    // geoSetDomainName(0, "PYLON 3 L"); // NUMERO
    // geoSetDomainName(0, "PYLON 3 R"); // NUMERO

    // geoSetDomainName(0, "TOP BALL 1"); // NUMERO
    // geoSetDomainName(0, "TOP BALL 2"); // NUMERO


    geoMeshWrite("../../../data/mesh.txt");

    /******************************/
    /* 2 : Definition du probleme */ 
    /******************************/

    double E_steel = 200.e9;
    double E_reinforced_concrete = 35e9;

    double nu_steel = 0.3;
    double nu_reinforced_concrete = 0.2;

    double rho_steel = 7.85e3;
    double rho_reinforced_concrete = 2.5e3;
    
    double gx = 0;
    double gy = -9.81;

    // TODO : Rendre possible le multi-materiaux
    femProblem *theProblem = femElasticityCreate(theGeometry, E_steel, nu_steel, rho_steel, gx, gy, PLANAR_STRAIN);

    // Boundary conditions

    femElasticityAddBoundaryCondition(theProblem, "PILAR R 1", DIRICHLET_XY, 0.0, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "PILAR L 1", DIRICHLET_XY, 0.0, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "PILAR D 1", DIRICHLET_XY, 0.0, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "PILAR R 2", DIRICHLET_XY, 0.0, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "PILAR L 2", DIRICHLET_XY, 0.0, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "PILAR D 2", DIRICHLET_XY, 0.0, 0.0);

    femElasticityAddBoundaryCondition(theProblem, "SUB ROADWAY L", DIRICHLET_XY, 0.0, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "SUB ROADWAY R", DIRICHLET_XY, 0.0, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "SUB ROADWAY D 1", NEUMANN_Y, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "SUB ROADWAY D 2", NEUMANN_Y, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "SUB ROADWAY D 3", NEUMANN_Y, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "SUB ROADWAY U 1", NEUMANN_Y, 800.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "SUB ROADWAY U 2", NEUMANN_Y, 1000.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "SUB ROADWAY U 3", NEUMANN_Y, 1200.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "SUB ROADWAY U 4", NEUMANN_Y, 700.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "SUB ROADWAY U 5", NEUMANN_Y, 20.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "SUB ROADWAY U 6", NEUMANN_Y, 1000.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "SUB ROADWAY U 7", NEUMANN_Y, 1300.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "SUB ROADWAY U 8", NEUMANN_Y, 800.0, NAN);

    femElasticityAddBoundaryCondition(theProblem, "ROADWAY L", DIRICHLET_XY, 0.0, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "ROADWAY R", DIRICHLET_XY, 0.0, 0.0);

    femElasticityAddBoundaryCondition(theProblem, "WINDOW D 1", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "WINDOW L 1", NEUMANN_Y, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "WINDOW R 1", NEUMANN_Y, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "WINDOW U 1", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "WINDOW D 2", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "WINDOW L 2", NEUMANN_Y, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "WINDOW R 2", NEUMANN_Y, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "WINDOW U 2", NEUMANN_X, 0.0, NAN);

    femElasticityAddBoundaryCondition(theProblem, "PILE L 1", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "PILE R 1", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "PILE L 2", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "PILE R 2", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "PILE L 3", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "PILE R 3", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "PILE L 4", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "PILE R 4", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "PILE L 5", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "PILE R 5", NEUMANN_X, 0.0, NAN);

    femElasticityAddBoundaryCondition(theProblem, "ARC 1", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "ARC 2", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "ARC 3", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "ARC 4", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "ARC 5", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "ARC 6", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "ARC 7", NEUMANN_X, 0.0, NAN);
    femElasticityAddBoundaryCondition(theProblem, "ARC 8", NEUMANN_X, 0.0, NAN);


    femElasticityPrint(theProblem);

    femElasticityWrite(theProblem, "../../../data/problem.txt");

    /***************************************************/
    /* 3 : Champ de la taille de référence du maillage */
    /*     (uniquement pour la visualisation)          */
    /***************************************************/

    double *meshSizeField = malloc(theGeometry->theNodes->nNodes * sizeof(double));
    if (meshSizeField == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }

    femNodes *theNodes = theGeometry->theNodes;
    for (int i = 0; i < theNodes->nNodes; ++i) { meshSizeField[i] = theGeometry->geoSize(theNodes->X[i], theNodes->Y[i]); }

    double hMin = femMin(meshSizeField, theNodes->nNodes);
    double hMax = femMax(meshSizeField, theNodes->nNodes);

    printf(" ==== Global requested h : %14.7e \n", theGeometry->h);
    printf(" ==== Minimum h          : %14.7e \n", hMin);
    printf(" ==== Maximum h          : %14.7e \n", hMax);

    /*********************/
    /* 4 : Visualisation */ 
    /*********************/

    if (!meshVisualizer) { return EXIT_SUCCESS; }

    int mode = 1;
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];

    GLFWwindow *window = glfemInit("EPL1110 : Project 2023-24 ");
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
        if (glfwGetKey(window, 'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t; }

        if (t - told > 0.5) { freezingButton = FALSE; }
        if (mode == 1)
        {
            glfemPlotField(theGeometry->theElements, meshSizeField);
            glfemPlotMesh(theGeometry->theElements);
            sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
            glColor3f(1.0, 0.0, 0.0);
            glfemMessage(theMessage);
        }
        else if (mode == 0)
        {
            domain %= theGeometry->nDomains;
            glfemPlotDomain(theGeometry->theDomains[domain]);
            sprintf(theMessage, "%s : %d ", theGeometry->theDomains[domain]->name, domain);
            glColor3f(1.0, 0.0, 0.0);
            glfemMessage(theMessage);
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) != 1);
    // Check if the ESC key was pressed or the window was closed

    free(meshSizeField);
    meshSizeField = NULL;
    femElasticityFree(theProblem);
    geoFree();
    glfwTerminate();

    exit(EXIT_SUCCESS);
    return EXIT_SUCCESS;
}