#ifndef FEAELEMENT_NURBS_2D_DERLAP_HPP
#define FEAELEMENT_NURBS_2D_DERLAP_HPP
// ==================================================================
// FEAElement_NURBS_2D_derLap.hpp
// It is a finite element implementation of 2D NURBS element, with
// 1st order derivatives and laplacian operator.
//
// Date: April 25 2014
// ==================================================================
#include <vector>
#include "Sys_Tools.hpp"
#include "ALocal_IEN.hpp"
#include "IALocal_meshSize.hpp"
#include "FEANode.hpp"
#include "IAExtractor.hpp"
#include "BernsteinBasis_Array.hpp"
#include "FEAElement.hpp"

class FEAElement_NURBS_2D_derLap : public FEAElement
{
  public:
    FEAElement_NURBS_2D_derLap( const int &in_eIndex,
        const IALocal_meshSize * const &mSize,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const FEANode * const &feaNode,
        const IAExtractor * const &extractor,
        const ALocal_IEN * const &locIEN );

    virtual ~FEAElement_NURBS_2D_derLap();

    virtual int get_elemIndex() const {return elem_index;}
    virtual int get_elemDim() const {return 2;}

    // Type = 204: 2D NURBS H1 conforming basis with 1st order 
    //             derivatives and laplacians.
    virtual int get_Type() const {return 204;}

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

    virtual void get_R_gradR( const int &quaindex,
        double * const &basis, double * const &basis_x,
        double * const &basis_y ) const;

    virtual void get_2D_R_gradR_LaplacianR( const int &quaindex,
        double * const &basis, double * const &basis_x, double * const &basis_y,
        double * const &basis_xx, double * const &basis_yy ) const;

    virtual double get_detJac(const int &quaindex) const
    {return detJac[quaindex];}

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
    double * d2R_dxx;
    double * d2R_dyy;
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
        double * &d2N_dss, double * &d2N_dtt, double * &d2N_dst,
        double * &dR_ds, double * &dR_dt,
        double * &d2R_dss, double * &d2R_dtt, double * &d2R_dst,
        const double * const &Ns,
        const double * const &Nt,
        const double * const &dNs_ds,
        const double * const &dNt_dt,
        const double * const &d2Ns_dss,
        const double * const &d2Nt_dtt );

};
#endif
