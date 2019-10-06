#ifndef FEAELEMENT_NURBS_3D_DER1_V2_HPP
#define FEAELEMENT_NURBS_3D_DER1_V2_HPP
// ==================================================================
// FEAElement_NURBS_3D_der1_v2.hpp
// It is a finite element implementation of the physical 3D NURBS 
// element, with derivative up to first order.
// 
// _v2: version 2. Explore fast element implementation.
//
// Date:
// Nov. 7th 2013
// ==================================================================
#include <vector>
#include <iostream>

#include "Sys_Tools.hpp"
#include "Vec_Tools.hpp"
#include "Matrix_double_3by3_Array.hpp"
#include "ALocal_IEN.hpp"
#include "IALocal_meshSize.hpp"
#include "FEANode.hpp"
#include "IAExtractor.hpp"
#include "ParentElement.hpp"
#include "FEAElement.hpp"

using namespace std;

class FEAElement_NURBS_3D_der1_v2 : public FEAElement
{
  public:
    FEAElement_NURBS_3D_der1_v2( const int &in_eIndex,
        const IALocal_meshSize * const &mSize,
        const ParentElement * const &pElem,
        const FEANode * const &feaNode,
        const IAExtractor * const &extractor,
        const ALocal_IEN * const &locIEN );

    virtual ~FEAElement_NURBS_3D_der1_v2();

    virtual int get_elemIndex() const {return elem_index;}
    virtual int get_elemDim() const {return 3;}

    // Type = 1: 3D NURBS H1 conforming basis with up to 1st order
    //           derivatives.
    virtual int get_Type() const {return 1;}

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

    virtual void get_dR_dx( const int &elemDir, 
        const int &quaindex, vector<double> &basis_x ) const;

    virtual void get_dR_dy( const int &elemDir, 
        const int &quaindex, vector<double> &basis_y ) const;

    virtual void get_dR_dz( const int &elemDir, 
        const int &quaindex, vector<double> &basis_z ) const;
    
    virtual void get_gradR( const int &elemDir, const int &quaindex, 
        vector<double> &basis_x, vector<double> &basis_y, 
        vector<double> &basis_z ) const;

    virtual double get_detJac(const int &quaindex) const;

    virtual void print() const;
    
    virtual double get_memory_usage() const;

  private:
    int elem_index;
    int numQuapts;
    int nLocBas;
    bool is_sNonzero;

    double * R;
    double * dR_dx;
    double * dR_dy;
    double * dR_dz;
    double * detJac;

    // This is a private funtion that evaluates shape functions at
    // given quadrature point index.
    void BuildShape_atQua( int quaindex, const ParentElement * const &pElem,
        double hx, double hy, double hz,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z,
        const double * const &ctrl_x, 
        const double * const &ctrl_y, 
        const double * const &ctrl_z, 
        const double * const &ctrl_w ); 
}; 
#endif
