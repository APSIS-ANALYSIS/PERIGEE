#ifndef PART_FEM_HPP
#define PART_FEM_HPP
// ============================================================================
// Part_FEM.hpp
//
// Object: Partition 3D finite element mesh into subdomains, and record each 
//         subdomain in hdf5 file format.
//
// Date: Jan. 12 2017
// ============================================================================
#include "Vec_Tools.hpp"
#include "IMesh.hpp"
#include "IPart.hpp"
#include "Map_Node_Index.hpp"
#include "IIEN.hpp"
#include "Field_Property.hpp"

class Part_FEM : public IPart
{
  public:
    Part_FEM( const IMesh * const &mesh,
        const IGlobal_Part * const &gpart,
        const Map_Node_Index * const &mnindex,
        const IIEN * const &IEN,
        const std::vector<double> &ctrlPts,
        const int &in_cpu_rank, const int &in_cpu_size,
        const int &in_elemType, 
        const Field_Property &in_fp );

    // Constructor that load the partition info from h5 file on disk
    Part_FEM( const std::string &fileName, const int &in_cpu_rank );

    virtual ~Part_FEM();

    virtual void write( const std::string &inputFileName ) const;
    
    virtual bool isElemInPart(const int &gloindex) const
    {return VEC_T::is_invec(elem_loc, gloindex);}
    
    virtual bool isNodeInPart(const int &gloindex) const
    {return VEC_T::is_invec(node_loc, gloindex);}
   
    // Determine the position of a given index in the elem_loc array 
    virtual int get_elemLocIndex(const int &gloindex) const
    {return VEC_T::get_pos(elem_loc, gloindex);}

    // Determine the position of a given index in the local_to_global array
    virtual int get_nodeLocGhoIndex(const int &gloindex) const
    {return VEC_T::get_pos(local_to_global, gloindex);}

    virtual void print_part_ele() const;

    virtual void print_part_node() const;

    virtual void print_part_ghost_node() const;

    virtual void print_part_local_to_global() const;

    virtual void print_part_LIEN() const;

    virtual void print_part_loadbalance_edgecut() const;

    virtual int get_elem_loc(const int &pos) const {return elem_loc[pos];}
    virtual int get_nlocalele() const {return nlocalele;}
    virtual int get_node_loc(const int &pos) const {return node_loc[pos];}
    virtual int get_node_loc_original(const int &pos) const {return node_loc_original[pos];}
    virtual int get_node_ghost(const int &pos) const {return node_ghost[pos];}
    virtual int get_local_to_global(const int &pos) const {return local_to_global[pos];}
    virtual int get_nlocalnode() const {return nlocalnode;}
    virtual int get_nghostnode() const {return nghostnode;}
    virtual int get_ntotalnode() const {return ntotalnode;}
    virtual int get_nbadnode() const {return nbadnode;}
    virtual int get_nlocghonode() const {return nlocghonode;}
    virtual int get_cpu_rank() const {return cpu_rank;}
    virtual int get_cpu_size() const {return cpu_size;}

    virtual int get_nElem() const {return nElem;}
    virtual int get_nFunc() const {return nFunc;}
    virtual int get_sDegree() const {return sDegree;}
    virtual int get_tDegree() const {return tDegree;}
    virtual int get_uDegree() const {return uDegree;}
    virtual int get_nLocBas() const {return nLocBas;}
    virtual int get_LIEN(const int &ee, const int &ii) const {return LIEN[ee][ii];}

    virtual double get_ctrlPts_x_loc(const int &pos) const {return ctrlPts_x_loc[pos];}
    virtual double get_ctrlPts_y_loc(const int &pos) const {return ctrlPts_y_loc[pos];}
    virtual double get_ctrlPts_z_loc(const int &pos) const {return ctrlPts_z_loc[pos];}

  protected:
    // ------------------------------------------------------------------------
    // Data
    // 1. local element
    std::vector<int> elem_loc {};
    int nlocalele;

    // 2. local node
    std::vector<int> node_loc {};
    std::vector<int> node_loc_original {};
    std::vector<int> node_ghost {};
    std::vector<int> local_to_global {};

    int nlocalnode, nghostnode, ntotalnode, nbadnode, nlocghonode;

    // 3. CPU info and partition parameters
    int cpu_rank, cpu_size;

    // 4. global mesh info
    int nElem, nFunc, sDegree, tDegree, uDegree, nLocBas;
    int probDim, dofNum, elemType;

    // 5. LIEN
    int ** LIEN;

    // 6. local point coordinates (i.e. control point geometry)
    const bool is_geo_field; // tag for geometry-related field
    std::vector<double> ctrlPts_x_loc {};
    std::vector<double> ctrlPts_y_loc {};
    std::vector<double> ctrlPts_z_loc {};

    // ------------------------------------------------------------------------
    // Function
    void Generate_Partition( const IMesh * const &mesh,
        const IGlobal_Part * const &gpart,
        const Map_Node_Index * const &mnindex,
        const IIEN * const &IEN,
        const int &field );
    
    Part_FEM() = delete;
};

#endif
