#ifndef MESHPART_NURBS_MULTIPATCH_3D_HPP
#define MESHPART_NURBS_MULTIPATCH_3D_HPP
// ============================================================================
// MeshPart_NURBS_multiPatch_3D.hpp
// This is a derived class from IMeshPart, which mainly reads in the global
// partition infomation and decides the local information for each subdomain
// (cpu).
//
// Date: Sept. 7 2015
// ============================================================================
#include "IMeshPart.hpp"

class MeshPart_NURBS_multiPatch_3D : public IMeshPart
{
  public:
    MeshPart_NURBS_multiPatch_3D( const IMesh * const &mesh,
        const IGlobal_Part * const &gpart,
        const int &in_cpu_size, const int &in_dofNum, const int &in_elemType );

    virtual ~MeshPart_NURBS_multiPatch_3D();

    virtual int get_cpu_rank() const {return cpu_rank;}

    virtual int get_nlocalnode() const {return nlocalnode;}

    virtual int get_nghostnode() const {return nghostnode;}

    virtual int get_ntotalnode() const {return ntotalnode;}

    virtual int get_nbadnode() const {return nbadnode;}

    virtual int get_nlocghonode() const {return nlocghonode;}
    
    virtual int get_local_to_global(const int &pos) const {return local_to_global[pos];}

    virtual int get_elemLocIndex(const int &pos) const;

    virtual bool isElemInPart(const int &eindex) const;

    virtual bool isNodeInPart(const int &nindex) const;


    virtual void set_part_info(const int &in_cpu_rank,
        const IMesh * const &mesh,
        const IGlobal_Part * const &gpart,
        const Map_Node_Index * const &mnindex,
        const IIEN * const &IEN,
        const std::vector<double> &ctrlPts,
        const bool &isPrintInfo  );

    virtual void write( const std::string &inputFileName,
        const IMesh * const &mesh,
        const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_x,
        const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_y,
        const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_z ) const;

    virtual int get_elem_loc(const int &pos) const {return elem_loc[pos];}

  private:
    const int numPat;

    std::vector<int> elem_loc;
    int nlocalele;

    std::vector<int> node_loc, node_loc_original, node_ghost, local_to_global;
    int nlocalnode, nghostnode, ntotalnode, nbadnode, nlocghonode;

    int cpu_rank;
    const int cpu_size, dual_edge_ncommon;

    const int nElem, nFunc, sDegree, tDegree, uDegree, nLocBas;
    const double hx_max, hy_max, hz_max, hx_min, hy_min, hz_min;

    std::vector<double> hx, hy, hz;

    const int probDim, dofNum;

    std::vector<int> LIEN;

    std::vector<double> ctrlPts_x_loc, ctrlPts_y_loc, ctrlPts_z_loc, ctrlPts_w_loc;

    const int elemType;

    // ------------------------------------------------------------------------
    // Private functions
    // ------------------------------------------------------------------------
    std::vector< std::vector<int> > nonzero_x, nonzero_y, nonzero_z;

    virtual void Get_bezier_ext( const int &elem, 
        const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_x,
        const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_y,
        const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_z,
        const IMesh * const &mesh, std::vector<double> &ext_x,
        std::vector<double> &ext_y, std::vector<double> &ext_z ) const;

};

#endif
