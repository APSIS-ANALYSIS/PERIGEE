#ifndef FEAELEMENT_NURBS_2D_DER0_HPP
#define FEAELEMENT_NURBS_2D_DER0_HPP
// ==================================================================
// FEAElement_NURBS_2D_der0.hpp
// It is a finite element implementation of 2D NURBS element. Only the
// function values are computed, which is useful for fast
// postprocessing.
//
// Date: April 23 2014
// ==================================================================
#include <vector>
#include "Sys_Tools.hpp"
#include "ALocal_IEN.hpp"
#include "IALocal_meshSize.hpp"
#include "FEANode.hpp"
#include "IAExtractor.hpp"
#include "BernsteinBasis_Array.hpp"
#include "FEAElement.hpp"

class FEAElement_NURBS_2D_der0 : public FEAElement
{
  public:
    FEAElement_NURBS_2D_der0( const int &in_eIndex,
        const IALocal_meshSize * const &mSize,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const FEANode * const &feaNode,
        const IAExtractor * const &extractor,
        const ALocal_IEN * const &locIEN );

    virtual ~FEAElement_NURBS_2D_der0();

    virtual int get_elemIndex() const {return elem_index;}
    virtual int get_elemDim() const {return 2;}

    // Type = 201: 2D NURBS H1 conforming basis with upto 1st order 
    //             derivatives
    virtual int get_Type() const {return 203;}

    // scalar type element
    virtual int get_numType() const {return 1;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual bool is_sizeNonzero() const {return is_sNonzero;}

    virtual void clearBasisCache();

    virtual void buildBasis( const int &ssdp1, const int &ttdp1,
        const int &num_qs, const int &num_qt,
        const IALocal_meshSize * const &mSize,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const FEANode * const &feaNode,
        const IAExtractor * const &extractor,
        const ALocal_IEN * const &locIEN );

    virtual int get_numQuapts() const {return numQuapts;}

    // get_R: return basis value at qua pt quaindex
    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual double get_detJac(const int &quaindex) const;

    virtual void print() const;

    virtual double get_memory_usage() const;

  private:
    int elem_index;
    int numQuapts;
    int nLocBas;
    bool is_sNonzero;

    double * R;
    double * detJac;
    
    // BuildShape_atQua: evaluate shape functions at given quadrature pt.
    void BuildShape_atQua(
        const int &ssdp1, const int &ttdp1,
        const int &quaindex,
        const int &quaindex_s, const int &quaindex_t,
        const double &hx, const double &hy,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z,
        const double * const &ctrl_w,
        double * &N, double * &dN_ds, double * &dN_dt,
        double * &dR_ds, double * &dR_dt,
        const double * const &Ns,
        const double * const &Nt,
        const double * const &dNs_ds,
        const double * const &dNt_dt );

};
#endif
