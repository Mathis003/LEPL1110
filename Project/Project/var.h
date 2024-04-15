#ifndef _VAR_H_
#define _VAR_H_

#include "../fem_library/include/fem.h"

femElementType elementType  = FEM_TRIANGLE; // FEM_QUAD or FEM_TRIANGLE
femDiscreteType discretType = FEM_DISCRETE_TYPE_LINEAR; // FEM_DISCRETE_TYPE_LINEAR or FEM_DISCRETE_TYPE_QUADRATIC
femElasticCase theCase      = PLANAR_STRESS; // PLANAR_STRESS or PLANAR_STRAIN or AXISYM (PLANAR_STRESS for our bridge problem)
femSolverType typeSolver    = FEM_BAND;  // FEM_FULL or FEM_BAND
femRenumType renumType      = FEM_RCMK;  // FEM_NO or FEM_XNUM or FEM_YNUM or FEM_RCMK

#endif // _VAR_H_