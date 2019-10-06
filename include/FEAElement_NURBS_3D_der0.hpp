#ifndef FEAELEMENT_NURBS_3D_DER0_HPP
#define FEAELEMENT_NURBS_3D_DER0_HPP
// ==================================================================
// FEAElement_NURBS_3D_der0.hpp
// It is a finite element implementation of the physical 3D NURBS 
// element, only the function values are calculated and cached.
// 
// Note: This is especially useful when doing postprocessing.
//
// Date:
// Nov. 19th 2013
// ==================================================================
#include <vector>
#include <iostream>

#include "Sys_Tools.hpp"
#include "Vec_Tools.hpp"
#include "ALocal_IEN.hpp"
#include "IALocal_meshSize.hpp"
#include "FEANode.hpp"
#include "IAExtractor.hpp"
#include "ParentElement.hpp"
#include "FEAElement.hpp"

using namespace std;

class FEAElement_NURBS_3D_der0 : public FEAElement
{
  public:
    FEAElement_NURBS_3D_der0( const int &in_eIndex,
        const IALocal_meshSize * const &mSize,
        const ParentElement * const &pElem,
        const FEANode * const &feaNode,
        const IAExtractor * const &extractor,
        const ALocal_IEN * const &locIEN );

    virtual ~FEAElement_NURBS_3D_der0();

    virtual int get_elemIndex() const {return elem_index;}
    
    // This is a three-dimensional element class
    virtual int get_elemDim() const {return 3;}

    // Type = 0: 3D NURBS H1 conforming basis with up to 1st order
    //           derivatives.
    virtual int get_Type() const {return 0;}

    // numType = 1: This is a basis for scalar variable.
    virtual int get_numType() const {return 1;}

    virtual int get_nLocBas(const int &elemDir)
      const {return nLocBas;}

    virtual bool is_sizeNonzero() const {return is_sNonzero;}

    virtual void clearBasisCache();
    virtual void buildBasis( const IALocal_meshSize * const &mSize, 
        const ParentElement * const &pElem, const FEANode * const &feaNode,
        const IAExtractor * const &extractor, 
        const ALocal_IEN * const &locIEN );

    virtual FEAElement * get_BCElement(const int &sideID) const;

    virtual int get_numQuapts() const {return numQuapts;}

    virtual void get_R(const int &elemDir,
        const int &quaindex, vector<double> &basis ) const;

    virtual double get_detJac(const int &quaindex) const;

    virtual void print() const;
    
    virtual double get_memory_usage() const;

  private:
    int elem_index;
    int numQuapts;
    int nLocBas;
    bool is_sNonzero;

    vector<double> R;
    vector<double> detJac;

    // This is a private funtion that evaluates shape functions at
    // given quadrature point index.
    void BuildShape_atQua( int quaindex, const ParentElement * const &pElem,
        double hx, double hy, double hz,
        const vector<double> &ext_x,
        const vector<double> &ext_y,
        const vector<double> &ext_z,
        vector<double> ctrl_x, vector<double> ctrl_y,
        vector<double> ctrl_z, vector<double> ctrl_w );
}; 
#endif
