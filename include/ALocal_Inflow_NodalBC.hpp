#ifndef ALOCAL_INFLOW_NODALBC_HPP
#define ALOCAL_INFLOW_NODALBC_HPP
// ============================================================================
// ALocal_Inflow_NodalBC.hpp
//
// Analysis-use data that is related to the inflow surface nodes.
//
// Date Created: July 8 2017
// Author: Ju Liu
// ============================================================================
#include "HDF5_Reader.hpp"
#include "Math_Tools.hpp"

class ALocal_Inflow_NodalBC
{
  public:
    ALocal_Inflow_NodalBC( const std::string &fileBaseName, const int &cpu_rank );

    virtual ~ALocal_Inflow_NodalBC();

    // Get global index of a Dirichlet node in the local partition
    // \para node ranges [ 0 , Num_LD ).
    virtual int get_LDN( const int &node ) const {return LDN[node];}

    // Get the number of Dirichlet nodes in the local partition
    virtual int get_Num_LD() const {return Num_LD;}

    // get the outward normal vector components.
    // ii=0 : x-component; ii=1 : y-component; ii=2 : z-component
    virtual double get_outvec( const int &ii ) const {return outnormal( ii );}

    // get the active area of the surface
    virtual double get_actarea() const {return act_area;}
    
    // get the full area of the surface
    virtual double get_fularea() const {return ful_area;}

    // determine whether a given index belongs to the LDN vector
    virtual bool is_inLDN( const int &ii) const
    {return VEC_T::is_invec(LDN, ii); }

    // ------------------------------------------------------------------------
    // get_radius: return the given point's (estimated) scaled radius for
    //             generating the inflow velocity profile.
    //             Algorithm: Find rc := point dist from centroid;
    //                        Find rb := point dist to nearest bc pt;
    //                        return rc / (rc + rb);
    //             If this partition does not contain any inflow bc
    //             nodes, this function will throw an error.
    // ------------------------------------------------------------------------
    virtual double get_radius( const double &x, const double &y, 
        const double &z ) const;

    // get number of nodes beloging to the local partition
    virtual int get_num_local_node() const {return num_local_node;}

    // get number of cells belonging to the local partition
    virtual int get_num_local_cell() const {return num_local_cell;}

    // get the number of nodes per element (nLocBas)
    virtual int get_cell_nLocBas() const {return cell_nLocBas;}

    // access coordinates of a node in the local partition by indexing
    // the local_pt_xyz array
    // 0 <= ii < 3 x num_local_node
    virtual double get_local_pt_xyz(const int &ii) const {return local_pt_xyz[ii];}

    // access coordinates of a node in the local partition using its 
    // local index and dir
    // 0 <= ii < num_local_node, 0 <= dir < 3
    virtual double get_local_pt_xyz(const int &ii, const int &dir) const 
    {return local_pt_xyz[3*ii+dir];}

    // access an element's IEN array by indexing the local_tri_ien array
    // 0 <= ii < cell_nLocBas x num_local_cell
    virtual int get_local_tri_ien(const int &ii) const {return local_tri_ien[ii];}

    // access an element's IEN array by its local element index and surface element
    // node index
    // 0 <= ee < num_local_cell, 0 <= ii < cell_nLocBas
    virtual int get_local_tri_ien(const int &ee, const int &ii) const
    {return local_tri_ien[ee * cell_nLocBas + ii];}

    // ------------------------------------------------------------------------
    // get_ctrlPts_xyz: given surface element index eindex, return the
    // control point coordinates.
    // Users are responsible for allocating/deallocating memory for ctrl_(x/y/z)
    // surface element id: 0 <= eindex < num_local_cell;
    // ctrl_x/y/z : output coordinate arrays, each of length cell_nLocBas.
    // ------------------------------------------------------------------------
    virtual void get_ctrlPts_xyz( const int &eindex, double * const &ctrl_x,
        double * const &ctrl_y, double * const &ctrl_z ) const;

    // ------------------------------------------------------------------------
    // get_SIEN: given surface element index eindex, return the IEN.
    // Users are responsible for allocating/deallocating memory for sien
    // eindex : 0 <= eindex < num_local_cell
    // sien : length cell_nLocBas.
    // ------------------------------------------------------------------------
    virtual void get_SIEN( const int &eindex, int * const &sien ) const;

    // ------------------------------------------------------------------------
    // Generate a filename for inlet data
    // Inlet_data.txt
    // ------------------------------------------------------------------------
    virtual std::string gen_flowfile_name() const
    {
      std::ostringstream ss;
      ss<<"Inlet_data.txt";
      return ss.str();
    }

  private:
    // Number of local Dirichlet nodes
    int Num_LD;

    // Global indices of the local Dirichlet nodes (inflow BC)
    std::vector<int> LDN;

    // Outward normal vector
    Vector_3 outnormal;

    // Inflow surface active area
    double act_area;

    // Inflow surface full area
    double ful_area;

    // Number of outer boundary points
    int num_out_bc_pts;

    // Coordinates of the outer boundary points
    std::vector<double> outline_pts;

    // Centroid point's coordinates
    Vector_3 centroid;

    // Number of local nodes & cells, and cell's nLocBas
    int num_local_node, num_local_cell, cell_nLocBas;

    // Nodal coordinates of all local nodes
    std::vector<double> local_pt_xyz;

    // Surface IEN array of all local elements
    std::vector<int> local_tri_ien;

    // Indices of all local nodes in the local volumetric partition's
    // local_to_global array
    std::vector<int> local_node_pos;
};

#endif
