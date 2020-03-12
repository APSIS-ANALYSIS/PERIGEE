#ifndef ALOCAL_INFLOW_NODALBC_HPP
#define ALOCAL_INFLOW_NODALBC_HPP
// ==================================================================
// ALocal_Inflow_NodalBC.hpp
//
// Analysis-use, inflow nodal indices.
//
// Date Created: July 8 2017
// Author: Ju Liu
// ==================================================================
#include "HDF5_Reader.hpp"
#include "Math_Tools.hpp"

class ALocal_Inflow_NodalBC
{
  public:
    ALocal_Inflow_NodalBC( const std::string &fileBaseName,
        const int &cpu_rank );

    virtual ~ALocal_Inflow_NodalBC();

    // Get the Dirichlet node's index in this local partition
    // para node ranges [ 0 , Num_LD ).
    virtual int get_LDN( const int &node ) const {return LDN[node];}

    // Get the number of Dirichlet nodes in this local partition
    virtual int get_Num_LD() const {return Num_LD;}

    // get the outward normal vector components.
    // ii=0 : x-component; ii=1 : y-component; ii=2 : z-component
    virtual double get_outvec( const int &ii ) const {return outvec[ii];}

    // get the active area of the surface
    virtual double get_actarea() const {return act_area;}
    
    // get the full area of the surface
    virtual double get_fularea() const {return ful_area;}

    // determine if a given index belong to LDN array
    virtual bool is_inLDN( const int &ii) const
    {return VEC_T::is_invec(LDN, ii); }

    // --------------------------------------------------------------    
    // get_radius: return the given point's (estimated) scaled radius for
    //             generating the inflow parabolic profile.
    //             Algorithm: Find rc := point dist from centroid;
    //                        Find rb := point dist to nearest bc pt;
    //                        return rc / (rc + rb);
    //             If this partition does not contain any inflow bc
    //             nodes, this function will throw and error.
    // --------------------------------------------------------------    
    virtual double get_radius( const double &x, const double &y,
                const double &z ) const;

  private:
    // Number of Local Dirichlet nodes
    int Num_LD;

    // Local Dirichlet Nodes (storing the inflow BC nodal indices)
    std::vector<int> LDN;

    // Outward Normal vector
    std::vector<double> outvec;

    // Inflow surface active area
    double act_area;

    // Inflow surface full area
    double ful_area;

    // Number of outer boundary points
    int num_out_bc_pts;

    // Coordinates of the outer boundary points
    std::vector<double> outline_pts;

    // Centroid point
    std::vector<double> centroid;

    // number of local nodes, cells, and cell's nLocBas
    int num_local_node, num_local_cell, cell_nLocBas;

    // local nodes' coordinates
    std::vector<double> local_pt_xyz;

    // local cell's IEN
    std::vector<int> local_tri_ien;
};

#endif
