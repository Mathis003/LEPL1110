
/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2024 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

static const double _gaussQuad4Xsi[4]    = {-0.577350269189626,-0.577350269189626, 0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Eta[4]    = { 0.577350269189626,-0.577350269189626,-0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Weight[4] = { 1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000};
static const double _gaussTri3Xsi[3]     = { 0.166666666666667, 0.666666666666667, 0.166666666666667};
static const double _gaussTri3Eta[3]     = { 0.166666666666667, 0.166666666666667, 0.666666666666667};
static const double _gaussTri3Weight[3]  = { 0.166666666666667, 0.166666666666667, 0.166666666666667};


femGeo* geoMeshCreate(const char *filename) 
{
   FILE* file = fopen(filename,"r");
   
   int trash, *elem;
   
   femGeo *theGeometry = malloc(sizeof(femGeo));
   femNodes *theNodes = malloc(sizeof(femNodes));
   theGeometry->theNodes = theNodes;
   ErrorScan(fscanf(file, "Number of nodes %d \n", &theNodes->nNodes));
   theNodes->X = malloc(sizeof(double)*(theNodes->nNodes));
   theNodes->Y = malloc(sizeof(double)*(theNodes->nNodes));
   for (int i = 0; i < theNodes->nNodes; i++) {
       ErrorScan(fscanf(file,"%d : %le %le \n",&trash,&theNodes->X[i],&theNodes->Y[i]));} 

   femMesh *theEdges = malloc(sizeof(femMesh));
   theGeometry->theEdges = theEdges;
   theEdges->nLocalNode = 2;
   theEdges->nodes = theNodes;
   ErrorScan(fscanf(file, "Number of edges %d \n", &theEdges->nElem));
   theEdges->elem = malloc(sizeof(int)*theEdges->nLocalNode*theEdges->nElem);
   for(int i=0; i < theEdges->nElem; ++i) {
        elem = theEdges->elem;
        ErrorScan(fscanf(file, "%6d : %6d %6d \n", &trash,&elem[2*i],&elem[2*i+1])); }
  
   femMesh *theElements = malloc(sizeof(femMesh));
   theGeometry->theElements = theElements;
   theElements->nLocalNode = 0;
   theElements->nodes = theNodes;
   char elementType[MAXNAME];  
   ErrorScan(fscanf(file, "Number of %s %d \n",elementType,&theElements->nElem));  
   if (strncasecmp(elementType,"triangles",MAXNAME) == 0) {
      theElements->nLocalNode = 3;
      theElements->elem = malloc(sizeof(int)*theElements->nLocalNode*theElements->nElem);
      for(int i=0; i < theElements->nElem; ++i) {
          elem = theElements->elem;
          ErrorScan(fscanf(file, "%6d : %6d %6d %6d \n", 
                    &trash,&elem[3*i],&elem[3*i+1],&elem[3*i+2])); }}
   if (strncasecmp(elementType,"quads",MAXNAME) == 0) {
      theElements->nLocalNode = 4;
      theElements->elem = malloc(sizeof(int)*theElements->nLocalNode*theElements->nElem);
      for(int i=0; i < theElements->nElem; ++i) {
          elem = theElements->elem;
          ErrorScan(fscanf(file, "%6d : %6d %6d %6d %6d \n", 
                    &trash,&elem[4*i],&elem[4*i+1],&elem[4*i+2],&elem[4*i+3])); }}
           
   ErrorScan(fscanf(file, "Number of domains %d\n", &theGeometry->nDomains));
   int nDomains = theGeometry->nDomains;
   theGeometry->theDomains = malloc(sizeof(femDomain*)*nDomains);
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = malloc(sizeof(femDomain)); 
      theGeometry->theDomains[iDomain] = theDomain;
      theDomain->mesh = theEdges; 
      ErrorScan(fscanf(file,"  Domain : %6d \n", &trash));
      ErrorScan(fscanf(file,"  Name : %[^\n]s \n", (char*)&theDomain->name));
      ErrorScan(fscanf(file,"  Number of elements : %6d\n", &theDomain->nElem));
      theDomain->elem = malloc(sizeof(int)*2*theDomain->nElem); 
      for (int i=0; i < theDomain->nElem; i++){
          ErrorScan(fscanf(file,"%6d",&theDomain->elem[i]));
          if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) ErrorScan(fscanf(file,"\n")); }}
    
   fclose(file);
   return theGeometry;
}

void geoMeshFree(femGeo* theGeometry) 
{
   femNodes *theNodes = theGeometry->theNodes;
   free(theNodes->X);
   free(theNodes->Y);
   free(theNodes);
   femMesh *theEdges = theGeometry->theEdges;
   free(theEdges->elem);
   free(theEdges);
   femMesh *theElements = theGeometry->theElements;
   free(theElements->elem);
   free(theElements);
   int nDomains = theGeometry->nDomains;
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = theGeometry->theDomains[iDomain];
      free(theDomain->elem);
      free(theDomain); }
   free(theGeometry->theDomains);
   free(theGeometry);
}

void geoMeshPrint(femGeo* theGeometry) 
{
   femNodes *theNodes = theGeometry->theNodes;
   if (theNodes != NULL) {
      printf("Number of nodes %d \n", theNodes->nNodes);
      for (int i = 0; i < theNodes->nNodes; i++) {
        printf("%6d : %14.7e %14.7e \n",i,theNodes->X[i],theNodes->Y[i]); }}
   femMesh *theEdges = theGeometry->theEdges;
   if (theEdges != NULL) {
     printf("Number of edges %d \n", theEdges->nElem);
     int *elem = theEdges->elem;
     for (int i = 0; i < theEdges->nElem; i++) {
        printf("%6d : %6d %6d \n",i,elem[2*i],elem[2*i+1]); }}
   femMesh *theElements = theGeometry->theElements;
   if (theElements != NULL) {
     if (theElements->nLocalNode == 3) {
        printf("Number of triangles %d \n", theElements->nElem);
        int *elem = theElements->elem;
        for (int i = 0; i < theElements->nElem; i++) {
            printf("%6d : %6d %6d %6d\n",i,elem[3*i],elem[3*i+1],elem[3*i+2]); }}
     if (theElements->nLocalNode == 4) {
        printf("Number of quads %d \n", theElements->nElem);
        int *elem = theElements->elem;
        for (int i = 0; i < theElements->nElem; i++) {
            printf("%6d : %6d %6d %6d %6d\n",i,elem[4*i],elem[4*i+1],elem[4*i+2],elem[4*i+3]); }}}
   int nDomains = theGeometry->nDomains;
   printf("Number of domains %d\n", nDomains);
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = theGeometry->theDomains[iDomain];
      printf("  Domain : %6d \n", iDomain);
      printf("  Name : %s\n", theDomain->name);
      printf("  Number of elements : %6d\n", theDomain->nElem);
      for (int i=0; i < theDomain->nElem; i++){
 //         if (i != theDomain->nElem  && (i % 10) != 0)  printf(" - ");
          printf("%6d",theDomain->elem[i]);
          if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) printf("\n"); }
      printf("\n"); }
  
  
}

femDomain* geoGetDomain(femGeo* theGeometry, char *name)
{
    int theIndex = -1;
    int nDomains = theGeometry->nDomains;
    for (int iDomain = 0; iDomain < nDomains; iDomain++) {
        femDomain *theDomain = theGeometry->theDomains[iDomain];
        if (strncasecmp(name,theDomain->name,MAXNAME) == 0)
            theIndex = iDomain;  }
    return theGeometry->theDomains[theIndex];
            
}

femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (type == FEM_QUAD && n == 4) {
        theRule->n      = 4;
        theRule->xsi    = _gaussQuad4Xsi;
        theRule->eta    = _gaussQuad4Eta;
        theRule->weight = _gaussQuad4Weight; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theRule->n      = 3;
        theRule->xsi    = _gaussTri3Xsi;
        theRule->eta    = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight; }
    else Error("Cannot create such an integration rule !");
    return theRule; 
}

void femIntegrationFree(femIntegration *theRule)
{
    free(theRule);
}

void _q1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  1.0;  eta[0] =  1.0;
    xsi[1] = -1.0;  eta[1] =  1.0;
    xsi[2] = -1.0;  eta[2] = -1.0;
    xsi[3] =  1.0;  eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;  
    phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
    phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
    phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] =   (1.0 + eta) / 4.0;  
    dphidxsi[1] = - (1.0 + eta) / 4.0;
    dphidxsi[2] = - (1.0 - eta) / 4.0;
    dphidxsi[3] =   (1.0 - eta) / 4.0;
    dphideta[0] =   (1.0 + xsi) / 4.0;  
    dphideta[1] =   (1.0 - xsi) / 4.0;
    dphideta[2] = - (1.0 - xsi) / 4.0;
    dphideta[3] = - (1.0 + xsi) / 4.0;

}

void _p1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  0.0;  eta[0] =  0.0;
    xsi[1] =  1.0;  eta[1] =  0.0;
    xsi[2] =  0.0;  eta[2] =  1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;  
    phi[1] = xsi;
    phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = -1.0;  
    dphidxsi[1] =  1.0;
    dphidxsi[2] =  0.0;
    dphideta[0] = -1.0;  
    dphideta[1] =  0.0;
    dphideta[2] =  1.0;

}


femDiscrete *femDiscreteCreate(int n, femElementType type)
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete));
    if (type == FEM_QUAD && n == 4) {
        theSpace->n       = 4;
        theSpace->x2      = _q1c0_x;
        theSpace->phi2    = _q1c0_phi;
        theSpace->dphi2dx = _q1c0_dphidx; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theSpace->n       = 3;
        theSpace->x2      = _p1c0_x;
        theSpace->phi2    = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx; }
    else Error("Cannot create such a discrete space !");
    return theSpace; 
}

void femDiscreteFree(femDiscrete *theSpace)
{
    free(theSpace);
}

void femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta)
{
    mySpace->x2(xsi,eta);
}

void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi)
{
    mySpace->phi2(xsi,eta,phi);
}

void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta)
{
    mySpace->dphi2dx(xsi,eta,dphidxsi,dphideta);
}

void femDiscretePrint(femDiscrete *mySpace)
{
    int i,j;
    int n = mySpace->n;
    double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];
    
    femDiscreteXsi2(mySpace,xsi,eta);
    for (i=0; i < n; i++) {
        
        femDiscretePhi2(mySpace,xsi[i],eta[i],phi);
        femDiscreteDphi2(mySpace,xsi[i],eta[i],dphidxsi,dphideta);

        for (j=0; j < n; j++)  {
            printf("(xsi=%+.1f,eta=%+.1f) : ",xsi[i],eta[i]);
            printf(" phi(%d)=%+.1f",j,phi[j]);  
            printf("   dphidxsi(%d)=%+.1f",j,dphidxsi[j]);  
            printf("   dphideta(%d)=%+.1f \n",j,dphideta[j]);  }
        printf(" \n"); }
}




femFullSystem *femFullSystemCreate(int size)
{
    femFullSystem *theSystem = malloc(sizeof(femFullSystem));
    femFullSystemAlloc(theSystem, size);
    femFullSystemInit(theSystem);

    return theSystem; 
}

void femFullSystemFree(femFullSystem *theSystem)
{
    free(theSystem->A);
    free(theSystem->B);
    free(theSystem);
}

void femFullSystemAlloc(femFullSystem *mySystem, int size)
{
    int i;  
    double *elem = malloc(sizeof(double) * size * (size+1)); 
    mySystem->A = malloc(sizeof(double*) * size); 
    mySystem->B = elem;
    mySystem->A[0] = elem + size;  
    mySystem->size = size;
    for (i=1 ; i < size ; i++) 
        mySystem->A[i] = mySystem->A[i-1] + size;
}

void femFullSystemInit(femFullSystem *mySystem)
{
    int i,size = mySystem->size;
    for (i=0 ; i < size*(size+1) ; i++) 
        mySystem->B[i] = 0;}


void femFullSystemPrint(femFullSystem *mySystem)
{
    double  **A, *B;
    int     i, j, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) {
        for (j=0; j < size; j++)
            if (A[i][j] == 0)  printf("         ");   
            else               printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
}

double* femFullSystemEliminate(femFullSystem *mySystem)
{
    double  **A, *B, factor;
    int     i, j, k, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    /* Gauss elimination */
    
    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-16 ) {
            printf("Pivot index %d  ",k);
            printf("Pivot value %e  ",A[k][k]);
            Error("Cannot eliminate with such a pivot"); }
        for (i = k+1 ; i <  size; i++) {
            factor = A[i][k] / A[k][k];
            for (j = k+1 ; j < size; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
    
    /* Back-substitution */
    
    for (i = size-1; i >= 0 ; i--) {
        factor = 0;
        for (j = i+1 ; j < size; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }
    
    return(mySystem->B);    
}

void  femFullSystemConstrain(femFullSystem *mySystem, 
                             int myNode, double myValue) 
{
    double  **A, *B;
    int     i, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0; }
    
    for (i=0; i < size; i++) 
        A[myNode][i] = 0; 
    
    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}


double femMin(double *x, int n) 
{
    double myMin = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMin = fmin(myMin,x[i]);
    return myMin;
}

double femMax(double *x, int n)
{
    double myMax = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMax = fmax(myMax,x[i]);
    return myMax;
}

void femError(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);                                                 
}

void femErrorScan(int test, int line, char *file)                                  
{ 
    if (test >= 0)  return;
    
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");   
    exit(69);                                       
}

void femWarning(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");                                              
}

