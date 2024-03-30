#ifndef _FEM_SOLVER_H_
#define _FEM_SOLVER_H_

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "fem_other.h"

typedef enum {FEM_FULL,FEM_BAND,FEM_ITER} femSolverType;
typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM} femRenumType;

typedef struct {
    double *B;
    double **A;
    int size;
} femFullSystem;

typedef struct {
    double *B;
    double **A;        
    int size;
    int band;        
} femBandSystem;

typedef struct {
    double *R;
    double *D;
    double *S;
    double *X; 
    double error;      
    int size;
    int iter;        
} femIterativeSolver;

typedef struct {
    femSolverType type;
    femFullSystem *local;
    void *solver;
} femSolver;

femFullSystem *femFullSystemCreate(int size);
void femFullSystemFree(femFullSystem *mySystem);
void femFullSystemAlloc(femFullSystem *mySystem, int size);
void femFullSystemInit(femFullSystem *mySystem);
void femFullSystemAssemble(femFullSystem *mySystem, double *Aloc, double *Bloc, int *map, int nLoc);
void femFullSystemConstrain(femFullSystem *mySystem, int myNode, double value);
double femFullSystemGet(femFullSystem* myFullSystem, int myRow, int myCol);
double *femFullSystemEliminate(femFullSystem *mySystem);
void femFullSystemPrint(femFullSystem *mySystem);
void femFullSystemPrintInfos(femFullSystem *mySystem);

femBandSystem *femBandSystemCreate(int size, int band);
void femBandSystemFree(femBandSystem *myBandSystem);
void femBandSystemInit(femBandSystem *myBandSystem);
int comparPosNode(const void *a, const void *b);
void femMeshRenumber(femMesh *theMesh, femRenumType renumType);
int femMeshComputeBand(femMesh *theMesh);
void femBandSystemAssemble(femBandSystem *myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc);
double femBandSystemGet(femBandSystem* myBandSystem, int myRow, int myCol);
double *femBandSystemEliminate(femBandSystem *myBand);
void femBandSystemPrint(femBandSystem *myBand);
void femBandSystemPrintInfos(femBandSystem *myBand);

femIterativeSolver *femIterativeSolverCreate(int size);
void femIterativeSolverFree(femIterativeSolver *mySolver);
void femIterativeSolverInit(femIterativeSolver *mySolver);
void femIterativeSolverPrint(femIterativeSolver *mySolver);
void femIterativeSolverPrintInfos(femIterativeSolver *mySolver);
int femIterativeSolverConverged(femIterativeSolver *mySolver);
void femIterativeSolverAssemble(femIterativeSolver *mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc);
double *femIterativeSolverEliminate(femIterativeSolver *mySolver);

femSolver *femSolverCreate(int sizeLoc);
femSolver *femSolverFullCreate(int size, int sizeLoc);
femSolver *femSolverBandCreate(int size, int sizeLoc, int band);
femSolver *femSolverIterativeCreate(int size, int sizeLoc);
void femSolverFree(femSolver *mySolver);
void femSolverInit(femSolver *mySolver);
double femSolverGet(femSolver *mySolver, int i, int j);
void femSolverPrint(femSolver *mySolver);
void femSolverPrintInfos(femSolver *mySolver);
void femSolverAssemble(femSolver *mySolver, double *Aloc, double *Bloc, double *Uloc,int *map, int nLoc);
double *femSolverEliminate(femSolver *mySolver);
int femSolverConverged(femSolver *mySolver);

#endif // _FEM_SOLVER_H_