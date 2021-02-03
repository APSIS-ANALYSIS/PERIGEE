#ifndef ALOCAL_EBC_HPP
#define ALOCAL_EBC_HPP
// ==================================================================
// ALocal_EBC.hpp
//
// Analysis-use local subdomain's elemental boundary conditions.
//
// Author: Ju Liu
// Date: Jan. 16 2017
// ==================================================================
#include "HDF5_Reader.hpp"

class ALocal_EBC
{
  public:
    // Constructor. Read from part file, and the EBC info is
    // stored in the group /ebc by default. Give a group name if
    // the data is written to a different groupname.
    ALocal_EBC( const std::string &fileBaseName,
        const int &cpu_rank, const std::string &gname="/ebc" );

    virtual ~ALocal_EBC();

    virtual void print_info() const
    {
      SYS_T::commPrint("---- ALocal_EBC: \n");
      SYS_T::commPrint("     num_ebc = %d \n", num_ebc);
    }

    virtual int get_num_ebc() const {return num_ebc;}

    // The following are functions that access the geometrical data of
    // the ii-th surface that is prescribed with the elemental BC
    // 0 <= ii < num_ebc
    virtual int get_num_local_node(const int &ii) const 
    {return num_local_node[ii];}

    // 0 <= ii < num_ebc
    virtual int get_num_local_cell(const int &ii) const
    {return num_local_cell[ii];}

    // 0 <= ii < num_ebc
    virtual int get_cell_nLocBas(const int &ii) const
    {return cell_nLocBas[ii];}

    // 0 <= ii < num_ebc
    // 0 <= jj < 3 x num_local_node[ii]
    virtual double get_local_pt_xyz(const int &ii, const int &jj) const
    {return local_pt_xyz[ii][jj];}

    // 0 <= ii < num_ebc
    // 0 <= jj < cell_nLocBas[ii] x num_local_cell[ii]
    virtual int get_local_tri_ien(const int &ii, const int &jj) const
    {return local_tri_ien[ii][jj];}

    // 0 <= ii < num_ebc
    // 0 <= jj < num_local_node[ii]
    virtual int get_local_global_node(const int &ii, const int &jj) const
    {return local_global_node[ii][jj];}

    // 0 <= ii < num_ebc
    // 0 <= jj < num_local_node[ii] 
    virtual int get_local_node_pos(const int &ii, const int &jj) const
    {return local_node_pos[ii][jj];}

    // 0 <= ii < num_ebc
    // 0 <= jj < num_local_cell[ii]
    virtual int get_local_global_cell(const int &ii, const int &jj) const
    {return local_global_cell[ii][jj];}

    // --------------------------------------------------------------
    // get_ctrlPts_xyz: given the ebc_id ii, the element index eindex,
    // return the control points' geometry.
    // Users are responsible for allocating & deleting the ctrl_xyz
    // arrays.
    // ebc_id : 0 <= ii < num_ebc;
    // surface element id: 0 <= eindex < num_local_cell[ii];
    // ctrl_x/y/z : output geometry array, length is 
    //              cell_nLocBas[ii].
    // --------------------------------------------------------------
    virtual void get_ctrlPts_xyz(const int &ii, 
        const int &eindex, double * const &ctrl_x, 
        double * const &ctrl_y, double * const &ctrl_z ) const;

    // --------------------------------------------------------------
    // get_SIEN: returns the surface element's IEN.
    // The users are responsible for allocating and deleting the sien
    // array.
    // ebc_id : 0 <= ii < num_ebc;
    // eindex : 0 <= eindex < num_local_cell[ii]
    // sien : length cell_nLocBas[ii].
    // --------------------------------------------------------------
    virtual void get_SIEN( const int &ii,
        const int &eindex, int * const &sien ) const;

    // --------------------------------------------------------------
    // get_intPts_xyz: returns the surface element's interior point
    // coordinates.
    // ebc_id : 0 <= ii < num_ebc;
    // eindex : 0 <= eindex < num_local_cell[ii];
    // coor_x/y/z : output interior point coordinates
    // --------------------------------------------------------------
    virtual void get_intPts_xyz(const int &ii,
        const int &eindex, double &coor_x, double &coor_y,
        double &coor_z ) const
    {
      SYS_T::print_fatal("Error: ALocal_EBC::get_intPts_xyz is not implemented. \n");
    }

    // --------------------------------------------------------------
    // Generate a file name for outlet face ii
    // Outlet_xx_flowrate.txt
    // --------------------------------------------------------------
    virtual std::string gen_flowfile_name(const int &ii) const
    {
      std::ostringstream ss;
      ss<<"Outlet_";
      if( ii/10 == 0 ) ss<<"00";
      else if(ii/100 == 0) ss<<"0";

      ss<<ii<<"_data.txt";
      
      return ss.str();
    }

    // --------------------------------------------------------------
    // The following functions are in the derived _outflow classes
    // for the Dirichlet-to-Neumann type outflow boundary conditions.
    // --------------------------------------------------------------
    // get_num_face_nodes : return the number of nodes on the surface
    //                      if this partition owns any cell on this
    //                      surface; otherwise return 0.
    //                      ii : face id ranging from 0 <= ii < num_ebc.
    // --------------------------------------------------------------
    virtual int get_num_face_nodes(const int &ii) const
    {
      SYS_T::print_fatal("Error: ALocal_EBC::get_num_face_nodes is not implemented. \n");
      return -1;
    }

    // --------------------------------------------------------------
    // get_intNA : returns the integral of N_A basis of the face,
    //             if this partition owns any cell on this surface.
    //             ii : face_id ranging 0 <= ii < num_ebc.
    //             out : length is get_num_face_nodes(ii).
    // --------------------------------------------------------------
    virtual void get_intNA(const int &ii, std::vector<double> &out) const
    {
      SYS_T::print_fatal("Error: ALocal_EBC::get_intNA is not implemented. \n");
    }
    
    // --------------------------------------------------------------
    // get_LID : returns the LID for the nodes associated with intNA,
    //           if this partition owns any cell on this surface.
    //           ii : face_id ranging from 0 <= ii < num_ebc,
    //           out : length is 3 x get_num_face_nodes(ii),
    //           and is the LID of node 0 x, y, z, node1 x, y, z, ...
    // --------------------------------------------------------------
    virtual void get_LID(const int &ii, std::vector<int> &out ) const
    {
      SYS_T::print_fatal("Error: ALocal_EBC::get_LID is not implemented. \n");
    }

    // --------------------------------------------------------------
    // get_outvec : return the outward normal vector for face ii, if
    //              this partition owns any cell on this surface.
    //              ii : face_id ranging from 0 <= ii < num_ebc,
    // --------------------------------------------------------------
    virtual void get_outvec( const int &ii, double &nx, double &ny, double &nz ) const
    {
      SYS_T::print_fatal("Error: ALocal_EBC::get_outvec is not implemented. \n");
    }

    // --------------------------------------------------------------
    // get_thickness : return the wall thickness, if this partition
    //                 owns any cell on the wall.
    // Users are responsible for allocating & deleting the e_thickness array.
    // Only one surface per the assumption in wall ebc.
    // surface element id: 0 <= eindex < num_local_cell[0]
    // e_thickness : output thickness array, length is cell_nLocBas[0].
    // --------------------------------------------------------------
    virtual void get_thickness( const int &eindex, double * const &e_thickness ) const
    {
      SYS_T::print_fatal("Error: ALocal_EBC::get_thickness is not implemented. \n");
    }

    // --------------------------------------------------------------
    // get_youngsmod : return the wall young's modulus, if this partition
    //                 owns any cell on the wall.
    // Users are responsible for allocating & deleting the e_youngsmod array.
    // Only one surface per the assumption in wall ebc.
    // surface element id: 0 <= eindex < num_local_cell[0]
    // e_youngsmod : output youngsmod array, length is cell_nLocBas[0].
    // --------------------------------------------------------------
    virtual void get_youngsmod( const int &eindex, double * const &e_youngsmod ) const
    {
      SYS_T::print_fatal("Error: ALocal_EBC::get_youngsmod is not implemented. \n");
    }

    // --------------------------------------------------------------
    // get_fluid_density : return the fluid density used to compute young's modulus 
    // --------------------------------------------------------------
    virtual double get_fluid_density() const
    {
      SYS_T::print_fatal("Error: ALocal_EBC::get_fluid_density is not implemented. \n");
      return -1.0;
    }

  protected:
    // the number of different ebc domain on which one may prescribe different
    // elemental boundary conditions.
    int num_ebc;

    // num_local_node[ii] gives the ii-th ebc's local node number
    // num_local_cell[ii] gives the ii-th ebc's local cell number
    // cell_nLocBas[ii] gives the cell's number of node. e.g., 
    //                  triangle surface is 3,
    //                  quadralaterial surface is 4,
    //                  quadratic triangle surface is 6.
    std::vector<int> num_local_node, num_local_cell, cell_nLocBas;

    // local_pt_xyz[ii] gives a list of local node's coordinates
    // size: num_ebc x (3 x num_local_node[ii])
    std::vector< std::vector<double> > local_pt_xyz;

    // local_tri_ien[ii] gives the local cell's IEN array
    // size: num_ebc x (cell_nLocBas[ii] x num_local_cell[ii]) 
    std::vector< std::vector<int> > local_tri_ien;

    // local nodes' global indices
    // size: num_ebc x num_local_node[ii]
    std::vector< std::vector<int> > local_global_node;

    // local node's position in the volumetric local portion's local_to_global array.
    // size: num_ebc x num_local_node[ii]
    std::vector< std::vector<int> > local_node_pos;

    // local cell's corresponding volumetric element indices
    // size: num_ebc x num_local_cell[ii]
    std::vector< std::vector<int> > local_global_cell;
};

#endif
