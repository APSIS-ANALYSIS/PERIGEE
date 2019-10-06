#ifndef KREFINEMENT_HPP
#define KREFINEMENT_HPP
// ========================================================
// kRefinement.hpp
// Object:
// Perfrom k-refinement for given geometry.
//
// Output:
// refined knot vectors and control points.
//
// Date: Sept. 19th 2013
// ========================================================
#include "Vec_Tools.hpp"
#include "NURBS_Tools.hpp"

// generate insert knot. evenly discribute the num_inserted into 
// each knot span.
void hRefine_newKnot_Generator(const std::vector<double> &existKnot,
    std::vector<double> &insertKnot, const int num_inserted,
    double &h_max, double &h_min ); 

// generate insert knot. insert num_inserted[ii] knots into the iith
// knot span.
void hRefine_newKnot_Generator(const std::vector<double> &existKnot,
    std::vector<double> &insertKnot, const std::vector<int> &num_inserted,
    double &h_max, double &h_min ); 

// check the format of existing knot vector and return the number
// of knot spans with nonzero measure if the format is correct.
int knotVec_check(const std::vector<double> &knotVec,
    const int degree );

// Use insertKnot and addDegree to do k refinement for the given mesh.
// Output: new knot vectors, and new control points.
void kRefinement( const int &addSDegree, 
    const int &addTDegree, const int &addUDegree,
    const std::vector<double> &insertSKnot,
    const std::vector<double> &insertTKnot,
    const std::vector<double> &insertUKnot,
    std::vector<double> &sknots, 
    std::vector<double> &tknots, 
    std::vector<double> &uknots,
    std::vector<double> &ctrlPts, 
    const int &dim, 
    int &sdegree, int &tdegree, int &udegree );


void kRefinement( const int &addSDegree, const int &addTDegree,
    const std::vector<double> &insertSKnot,
    const std::vector<double> &insertTKnot,
    std::vector<double> &sknots, std::vector<double> &tknots,
    std::vector<double> &ctrlPts, 
    const int &dim, 
    int &sdegree, int &tdegree );


// Use addDegree to do p-refinement as a special case of
// k-refinement.
// The knots in the given knot vector will be added `addXDegree'
// times to maintain the continuity across knots.
// Then degree elevation will be called.
void pRefinement( const int &dim, const int &addSDegree,
    const int &addTDegree, const int &addUDegree,
    std::vector<double> &sKnots,
    std::vector<double> &tKnots,
    std::vector<double> &uKnots,
    std::vector<double> &ctrlPts,
    int &sdegree, int &tdegree, int &udegree );

#endif
