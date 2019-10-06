#ifndef FEAELEMENT_NURBS_3D_DER1_V3_HPP
#define FEAELEMENT_NURBS_3D_DER1_V3_HPP
// ==================================================================
// FEAElement_NURBS_3D_der1_v3.hpp
// It is a finite element implementation of the physical 3D NURBS 
// element, with derivative up to first order.
// 
// _v3: version 3. Explore fast element implementation: This routine
//      is based on univariate bernstein basis instead of 
//      parentElement.
//
// Date:
// Nov. 21th 2013
// ==================================================================
#include "Sys_Tools.hpp"
#include "ALocal_IEN.hpp"
#include "IALocal_meshSize.hpp"
#include "FEANode.hpp"
#include "IAExtractor.hpp"
#include "BernsteinBasis_Array.hpp"
#include "FEAElement.hpp"

class FEAElement_NURBS_3D_der1_v3 : public FEAElement
{
  public:
    FEAElement_NURBS_3D_der1_v3( const int &in_eIndex,
        const IALocal_meshSize * const &mSize,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const FEANode * const &feaNode,
        const IAExtractor * const &extractor,
        const ALocal_IEN * const &locIEN );

    virtual ~FEAElement_NURBS_3D_der1_v3();

    virtual int get_elemIndex() const {return elem_index;}
    virtual int get_elemDim() const {return 3;}

    // Type = 304: 3D NURBS H1 conforming basis with up to 1st order
    //           derivatives.
    virtual int get_Type() const {return 304;}

    // numType = 1: This is a basis for scalar variable.
    virtual int get_numType() const {return 1;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual bool is_sizeNonzero() const {return is_sNonzero;}

    virtual void clearBasisCache();
    virtual void buildBasis( 
        const int &ssdp1, const int &ttdp1, const int &uudp1,
        const int &num_qs, const int &num_qt, const int &num_qu,
        const IALocal_meshSize * const &mSize, 
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const FEANode * const &feaNode,
        const IAExtractor * const &extractor, 
        const ALocal_IEN * const &locIEN );

    virtual int get_numQuapts() const {return numQuapts;}

    virtual void get_R(const int &quaindex, 
        double * const &basis ) const;

    virtual void get_gradR( const int &quaindex, 
        double * const &basis_x, double * const &basis_y, 
        double * const &basis_z ) const;

    virtual void get_R_gradR( const int &quaindex, 
        double * const &basis, double * const &basis_x, 
        double * const &basis_y, double * const &basis_z ) const;

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
    void BuildShape_atQua( 
        const int &ssdp1, const int &ttdp1, const int &uudp1,
        const int &quaindex,
        const int &quaindex_s, const int &quaindex_t, const int &quaindex_u,
        const double &hx, const double &hy, const double &hz,
        const double * const &ctrl_x, 
        const double * const &ctrl_y, 
        const double * const &ctrl_z, 
        const double * const &ctrl_w,
        double * &N, double * &dN_ds, double * &dN_dt, double * &dN_du,
        double * &dR_ds, double * &dR_dt, double * &dR_du,
        const double * const &Ns,
        const double * const &Nt,
        const double * const &Nu,
        const double * const &dNs_ds,
        const double * const &dNt_dt, 
        const double * const &dNu_du ); 
}; 
#endif
