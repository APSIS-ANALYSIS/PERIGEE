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

    // get the number of nodes beloging to this local partition
    virtual int get_num_local_node() const {return num_local_node;}

    // get the number of cells belonging to this local partition
    virtual int get_num_local_cell() const {return num_local_cell;}

    // get the element's number of nodes
    virtual int get_cell_nLocBas() const {return cell_nLocBas;}

    // get the coordinates of the nodes
    // direct access the data array local_pt_xyz by the naive array index
    // 0 <= ii < 3 x num_local_node
    virtual double get_local_pt_xyz(const int &ii) const {return local_pt_xyz[ii];}

    // access the point coordinates by the node surface index and dir
    // 0 <= ii < num_local_node, 0 <= dir < 3
    virtual double get_local_pt_xyz(const int &ii, const int &dir) const 
    {return local_pt_xyz[3*ii+dir];}

    // get the element's IEN array
    // direct access the data array local_tri_ien by its naive array index
    // 0 <= ii < cell_nLocBas x num_local_cell
    virtual int get_local_tri_ien(const int &ii) const {return local_tri_ien[ii];}

    // access the ien array by the local element index and surface element local
    // node index
    // 0 <= ee < num_local_cell, 0 <= ii < cell_nLocBas
    virtual int get_local_tri_ien(const int &ee, const int &ii) const
    {return local_tri_ien[ee * cell_nLocBas + ii];}

    // --------------------------------------------------------------
    // get_ctrlPts_xyz: given element index eindex,
    // return the control point's geometry.
    // The users are responsible for allocating the deleting the ctrl_xyz
    // array.
    // surface element id: 0 <= eindex < num_local_cell;
    // ctrl_x/y/z : output geometry array, length is cell_nLocBas.
    // --------------------------------------------------------------
    virtual void get_ctrlPts_xyz( const int &eindex, double * const &ctrl_x,
        double * const &ctrl_y, double * const &ctrl_z ) const;

    // --------------------------------------------------------------
    // get_SIEN: returns the surface element's IEN.
    // The users are responsible for allocating and deleting the sien
    // array.
    // eindex : 0 <= eindex < num_local_cell
    // sien : length cell_nLocBas.
    // --------------------------------------------------------------
    virtual void get_SIEN( const int &eindex, int * const &sien ) const;

    // --------------------------------------------------------------
    // Generate a file name for inlet face
    // Inlet_data.txt
    // --------------------------------------------------------------
    virtual std::string gen_flowfile_name() const
    {
      std::ostringstream ss;
      ss<<"Inlet_data.txt";

      return ss.str();
    }

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

    // local node's position in the volumetric local portion's local_to_global
    // array
    std::vector<int> local_node_pos;
};

#endif
