#ifndef PART_NURBS_1PATCH_3D_METIS_HPP
#define PART_NURBS_1PATCH_3D_METIS_HPP
// ==================================================================
// Part_NURBS_1Patch_3D_METIS.hpp
// Object:
// Partition the 3D NURBS single patch mesh into subdomains, using 
// METIS library.
//
// Date: Oct. 1 2013
// ==================================================================
#include <cassert>
#include "NURBS_Bezier.hpp"
#include "IMesh.hpp"
#include "IPart.hpp"
#include "HDF5_Writer.hpp"

class Part_NURBS_1Patch_3D_METIS : public IPart
{
  public:
    Part_NURBS_1Patch_3D_METIS( const class IMesh * const &mesh,
        const class IGlobal_Part * const &gpart,
        const class Map_Node_Index * const &mnindex,
        const class IIEN * const &IEN,
        const std::vector<double> &ctrlPts,
        const std::vector<NURBS_T::BezierElem*> &in_seg_x,
        const std::vector<NURBS_T::BezierElem*> &in_seg_y,
        const std::vector<NURBS_T::BezierElem*> &in_seg_z,
        const int &in_cpu_rank, const int &in_cpu_size,
        const int &in_dofNum, const int &in_elemType,
        const bool isPrintInfo );

    virtual ~Part_NURBS_1Patch_3D_METIS();

    // 1. function access element partition.
    virtual s_int get_elem_loc(int pos) const {return elem_loc[pos];}
    virtual int get_nlocalele() const {return nlocalele;}

    // 2. function access node partition
    virtual s_int get_node_loc(int pos) const {return node_loc[pos];}
    virtual s_int get_node_loc_original(int pos) const {return node_loc_original[pos];}
    virtual s_int get_node_ghost(int pos) const {return node_ghost[pos];}
    virtual s_int get_local_to_global(int pos) const {return local_to_global[pos];}
    
    virtual int get_nlocalnode() const {return nlocalnode;}
    virtual int get_nghostnode() const {return nghostnode;}
    virtual int get_ntotalnode() const {return ntotalnode;}
    virtual int get_nbadnode() const {return nbadnode;}
    virtual int get_nlocghonode() const {return nlocghonode;}

    virtual bool isElemInPart(s_int gloindex) const;
    virtual bool isNodeInPart(s_int gloindex) const;
    
    virtual int get_elemLocIndex(const s_int &gloindex) const;

    // 3. function that return the partition method parameters
    virtual int get_cpu_rank() const {return cpu_rank;}
    virtual int get_cpu_size() const {return cpu_size;}
    virtual bool get_part_isDual() const {return part_isdual;}
    virtual bool get_isMETIS() const {return isMETIS;}
    virtual int get_dual_edge_ncommon() const {return dual_edge_ncommon;}


    // 3. Global mesh information
    virtual s_int get_nElem() const {return nElem;}
    virtual s_int get_nElem_x() const {return nElem_x;}
    virtual s_int get_nElem_y() const {return nElem_y;}
    virtual s_int get_nElem_z() const {return nElem_z;}

    virtual s_int get_nFunc() const {return nFunc;}
    virtual s_int get_nFunc_x() const {return nFunc_x;}
    virtual s_int get_nFunc_y() const {return nFunc_y;}
    virtual s_int get_nFunc_z() const {return nFunc_z;}

    virtual l_int get_sDegree() const {return sDegree;}
    virtual l_int get_tDegree() const {return tDegree;}
    virtual l_int get_uDegree() const {return uDegree;}

    virtual l_int get_nLocBas() const {return nLocBas;}

    virtual double get_hx_max() const {return hx_max;}
    virtual double get_hy_max() const {return hy_max;}
    virtual double get_hz_max() const {return hz_max;}

    virtual double get_hx_min() const {return hx_min;}
    virtual double get_hy_min() const {return hy_min;}
    virtual double get_hz_min() const {return hz_min;}

    virtual double get_hx(int ee) const {return hx[ee];}
    virtual double get_hy(int ee) const {return hy[ee];}
    virtual double get_hz(int ee) const {return hz[ee];}

    // 4. LIEN array
    virtual int get_LIEN(int ee, int ii) const {return LIEN[ee][ii];}

    // 5. Control points in this processor
    virtual double get_ctrlPts_x_loc(int pos) const {return ctrlPts_x_loc[pos];}
    virtual double get_ctrlPts_y_loc(int pos) const {return ctrlPts_y_loc[pos];}
    virtual double get_ctrlPts_z_loc(int pos) const {return ctrlPts_z_loc[pos];}
    virtual double get_ctrlPts_w_loc(int pos) const {return ctrlPts_w_loc[pos];}


    // 6. print on screen the partition results 
    virtual void print_info() const;
    virtual void print_part_ele() const;
    virtual void print_part_node() const;
    virtual void print_part_ghost_node() const;
    virtual void print_part_local_to_global() const;
    virtual void print_part_LIEN() const;
    virtual void print_part_loadbalance_edgecut() const;

    virtual void print_part_bezier_ext_x(int elem) const;
    virtual void print_part_bezier_ext_y(int elem) const;
    virtual void print_part_bezier_ext_z(int elem) const;

    // 7. write info
    virtual void write( const char * inputFileName ) const; 

  private:
    // ----------------------- Private Data -------------------------------
    // 1. local element
    std::vector<s_int> elem_loc;
    int nlocalele;

    // 2. local node 
    std::vector<s_int> node_loc;
    std::vector<s_int> node_loc_original;
    std::vector<s_int> node_ghost;
    std::vector<s_int> local_to_global;
    
    int nlocalnode;
    int nghostnode;
    int ntotalnode;
    int nbadnode;
    int nlocghonode; 
    
    // 3. cpu info and partition key parameters
    int cpu_rank, cpu_size;
    bool isMETIS, part_isdual;
    int dual_edge_ncommon;

    // 4. global mesh info
    s_int nElem, nElem_x, nElem_y, nElem_z;
    s_int nFunc, nFunc_x, nFunc_y, nFunc_z;
    int sDegree, tDegree, uDegree;
    int nLocBas;
    
    double hx_max, hy_max, hz_max, hx_min, hy_min, hz_min;

    std::vector<double> hx, hy, hz;
    
    int probDim, dofNum;

    // 5. LIEN
    int ** LIEN;

    // 6. local control points
    std::vector<double> ctrlPts_x_loc, ctrlPts_y_loc, ctrlPts_z_loc, ctrlPts_w_loc;
    
    // 7. pointer to extraction operator objects
    std::vector<NURBS_T::BezierElem*> pseg_x, pseg_y, pseg_z;
    std::vector<int> nonzero_x, nonzero_y, nonzero_z;
    int elemType;

    // ----------------------- Private Data -------------------------------

    // Private Function:
    // 1. Get global mesh information
    void Get_global_meshinfo( const class IMesh * const &mesh );

    // 2. Generate partition by calling METIS 
    void Generate_Partition( const class IMesh * const &mesh,
       const class IGlobal_Part * const &gpart,
       const class Map_Node_Index * const &mnindex,
       const class IIEN * const &IEN,
       const std::vector<double> &ctrlPts,
       const bool isPrintinfo );

    // 3. Private writers for all private data.
    // 3.1 writer for local element
    void write_local_elem( hid_t file_id ) const;

    // 3.2 writer for local node index
    void write_local_node( hid_t file_id ) const;
    
    // 3.3 writer for cpu info and partition parameters
    void write_part_info( hid_t file_id ) const;
    
    // 3.4 writer for global mesh info
    void write_global_mesh_info( hid_t file_id ) const;

    // 3.5 writer for LIEN
    void write_LIEN( hid_t file_id ) const;

    // 3.6 writer for local control points
    void write_ctrlPts_loc( hid_t file_id ) const;

    // 3.7 write extraction operators
    void write_extractor( hid_t file_id ) const;

    // 3.8 write mesh size hx hy hz
    void write_hxyz( hid_t file_id ) const;

    // 4. Private get function to obtain the extraction operator 
    //    for each local element.
    void Get_bezier_ext_x(const int &elem_x, std::vector<double> &ext_x) const;
    void Get_bezier_ext_y(const int &elem_y, std::vector<double> &ext_y) const;
    void Get_bezier_ext_z(const int &elem_z, std::vector<double> &ext_z) const;
};

#endif
