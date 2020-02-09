#ifndef AINT_WEIGHT_HPP
#define AINT_WEIGHT_HPP
// ==================================================================
// AInt_Weight.hpp
// Analysis Tools: Integration Weights.
// Stores the Integration weights.
//
// Note: This class is used if ParentElement is not used. This class
//       together with BernsteinBasis class perofrm the same job as
//       ParentElement class.
//
// Date: Nov. 23 2013.
// ==================================================================
#include "IQuadPts.hpp"

class AInt_Weight
{
	public:
		AInt_Weight( const IQuadPts * const &qua_s,
				const IQuadPts * const &qua_t,
				const IQuadPts * const &qua_u );
		
    AInt_Weight( const IQuadPts * const &qua_s,
				const IQuadPts * const &qua_t );
		
    // Obtain the weights from a single quadrature rule.
    // This can be used for initializing a weight object for 1D problem
    // or can be used for initializing a weight object for a multi-dim
    // problem with all quadrature rule written in one quad file, like
    // the Quad_Gauss_Tet class.
    AInt_Weight( const IQuadPts * const &qua );
		
    ~AInt_Weight();

    double get_weight(const int &index) const {return Weight[index];}

    int get_num() const {return num;}

    void print_info() const;

	private:
		int num; // number of the quadrature points/weights
		double * Weight; // dynamic array storing the weights
};

#endif
