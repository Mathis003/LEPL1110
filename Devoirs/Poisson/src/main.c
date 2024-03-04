/*
 *  main.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2022 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

#include "glfem.h"



int main(void)
{   
    // char *path_filename = "../data/mesh2520.txt";
    // char *path_filename = "../data/tiny.txt";
    // char *path_filename = "../data/example.txt";
    // char *path_filename = "../data/mesh254.txt";
    char *path_filename = "../data/quads418.txt";

 
    femPoissonProblem *theProblem = femPoissonCreate(path_filename);
    femPoissonFindBoundaryNodes(theProblem);
    // geoMeshPrint(theProblem->geo);
    
    // Pour Windows, remplacer l'argument :
    // ("../data/triangles_166.txt") 
    // par :
    // ("..\\data\\triangles_166.txt") 
    //
    // Sorry for the inconvenience :-)
    // On reflechit pour rendre cela plus transparent dans les homeworks suivants :-)
    // Be patient !
    
    
    femMesh *theMesh = theProblem->geo->theElements;
    femDomain *theBoundary = geoGetDomain(theProblem->geo,"Boundary");
    printf("Number of elements       : %4d\n", theMesh->nElem);
    printf("Number of local nodes    : %4d\n", theMesh->nLocalNode);
    printf("Number of boundary nodes : %4d\n", theBoundary->nElem);
    printf("Number of unknowns       : %4d\n", theProblem->system->size);


    femPoissonSolve(theProblem);   
 
    printf("Maximum value : %.4f\n", femMax(theProblem->system->B,theProblem->system->size));
    fflush(stdout);
    
    char theMessage[256];
    sprintf(theMessage, "Max : %.4f", femMax(theProblem->system->B,theProblem->system->size));
  
    GLFWwindow* window = glfemInit("EPL1110 : Poisson");
    
    
    glfwMakeContextCurrent(window);
    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theMesh->nodes,w,h);
        glfemPlotField(theMesh,theProblem->system->B);            
        glColor3f(1.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);              
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
           
    // Check if the ESC key was pressed or the window was closed
            
    femPoissonFree(theProblem);
    exit(EXIT_SUCCESS);
}