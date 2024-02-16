/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Utilisation de l'API de GMSH pour créer un maillage
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "glfem.h"


int main(void)
{  
    printf("\n\n    V : Mesh and size mesh field \n");
    printf("    D : Domains \n");
    printf("    N : Next domain highlighted\n");



 
    double Lx = 1.0;
    double Ly = 2.0;
      
    int ierr;
    
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
    
    theGeometry->LxPlate =  Lx;
    theGeometry->LyPlate =  Ly;
    theGeometry->xPlate  =  0.0;
    theGeometry->yPlate  =  0.0;
    theGeometry->xHole   =  Lx / 4.0;
    theGeometry->yHole   =  Ly / 4.0;
    theGeometry->rHole   =  Lx / 8.0;
    theGeometry->xNotch  = -Lx / 2.0;
    theGeometry->yNotch  = -Ly / 2.0;
    theGeometry->rNotch  =  Lx / 2.0;
    
    theGeometry->h       =  Lx * 0.1;    // 0.1 c'est plus joli :-)
    theGeometry->hHole   =  theGeometry->h * 0.2;
    theGeometry->hNotch  =  theGeometry->h * 0.05;
    theGeometry->dHole   =  theGeometry->h * 1.0;
    theGeometry->dNotch  =  theGeometry->h * 4.0;

   
    geoMeshGenerate();
    geoMeshImport();
    
    geoSetDomainName(0,"Outer Disk");
    geoSetDomainName(1,"Bottom");
    geoSetDomainName(2,"Left");
    geoSetDomainName(3,"Right");
    geoSetDomainName(4,"Top");
    geoSetDomainName(5,"Inner Disk");
    

//
//  -2- Creation du fichier du maillage
//
    
    char filename[] = "../data/mesh.txt";
    geoMeshWrite(filename);

//
//  -3- Champ de la taille de référence du maillage
//

    double *meshSizeField = malloc(theGeometry->theNodes->nNodes*sizeof(double));
    femNodes *theNodes = theGeometry->theNodes;
    for(int i=0; i < theNodes->nNodes; ++i)
        meshSizeField[i] = geoSize(theNodes->X[i], theNodes->Y[i]);
    double hMin = femMin(meshSizeField,theNodes->nNodes);  
    double hMax = femMax(meshSizeField,theNodes->nNodes);  
    printf(" ==== Global requested h : %14.7e \n",theGeometry->h);
    printf(" ==== Minimum h          : %14.7e \n",hMin);
    printf(" ==== Maximum h          : %14.7e \n",hMax);
 
//
//  -4- Visualisation du maillage
//  
    
    int mode = 1; // Change mode by pressing "j", "k", "l"
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[256];
    double pos[2] = {20,460};
 
 
    GLFWwindow* window = glfemInit("EPL1110 : Mesh generation ");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
    
    
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
    //    glfemChangeState(&mode, theMeshes->nMesh);
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}

        
        if (t-told > 0.5) {freezingButton = FALSE; }
            
        
        
         
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements, meshSizeField);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
 
            
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
 
            
            
            }
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            
            
            
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
 
            
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
            }
            
         glfwSwapBuffers(window);
         glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    free(meshSizeField);  
    geoFinalize();
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}

 
