#ifndef ALOCAL_EBC_HPP
#define ALOCAL_EBC_HPP
// ============================================================================
// ALocal_EBC.hpp
//
// Analysis-use local subdomain's elemental boundary conditions.
//
// Author: Ju Liu
// Date: Jan. 16 2017
// ============================================================================
#include "HDF5_Reader.hpp"
#include "Vector_3.hpp"

class ALocal_EBC
{
  public:
    // ------------------------------------------------------------------------
    // ! Constructor. 
    //   Read from part file, and the EBC info is stored in the group /ebc by 
    //   default. User may specify a group name for gname if the data is written
    //   under a different groupname.
    // ------------------------------------------------------------------------
    ALocal_EBC( const std::string &fileBaseName, const int &cpu_rank, 
        const std::string &gname="/ebc" );

    virtual ~ALocal_EBC() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("---- ALocal_EBC: \n");
      SYS_T::commPrint("     num_ebc = %d \n", num_ebc);
    }

    // ------------------------------------------------------------------------
    // ! Get the number of different elemental surfaces, that may potentially
    //   be associated with different material properties or boundary tractions.
    // ------------------------------------------------------------------------
    virtual int get_num_ebc() const {return num_ebc;}

    // ------------------------------------------------------------------------
    // The following are functions that access the geometrical data of
    // the ii-th surface that is prescribed with the elemental BC
    // 0 <= ii < num_ebc
    // ! get the number of all nodes associated with the surface cells within 
    //   this partition.
    // ------------------------------------------------------------------------
    virtual int get_num_local_cell_node(const int &ii) const 
    {return num_local_cell_node[ii];}

    // ------------------------------------------------------------------------
    // ! get the number of surface cells within this partition.
    //   \para 0 <= ii < num_ebc
    // ------------------------------------------------------------------------
    virtual int get_num_local_cell(const int &ii) const {return num_local_cell[ii];}

    // ------------------------------------------------------------------------
    // ! get the number of local basis functions of the surface cells. This is
    //   directly associated with the cell element type.
    //   \para 0 <= ii < num_ebc
    // ------------------------------------------------------------------------
    virtual int get_cell_nLocBas(const int &ii) const {return cell_nLocBas[ii];}

    // ------------------------------------------------------------------------
    // ! get the local cell node's spatial coordinates.
    //   \para 0 <= ii < num_ebc
    //   \para 0 <= jj < 3 x num_local_cell_node[ii]
    // ------------------------------------------------------------------------
    virtual double get_local_cell_node_xyz(const int &ii, const int &jj) const
    {return local_cell_node_xyz[ii][jj];}

    // ------------------------------------------------------------------------
    // ! get the local cell's IEN connectivity array, with ranges in the local
    //   node array of this class.
    //   \para 0 <= ii < num_ebc
    //   \para 0 <= jj < cell_nLocBas[ii] x num_local_cell[ii]
    // ------------------------------------------------------------------------
    virtual int get_local_cell_ien(const int &ii, const int &jj) const
    {return local_cell_ien[ii][jj];}

    // ------------------------------------------------------------------------
    // ! get the local cell node's volumetric mesh index
    //   \para 0 <= ii < num_ebc
    //   \para 0 <= jj < num_local_cell_node[ii]
    // ------------------------------------------------------------------------
    virtual int get_local_cell_node_vol_id(const int &ii, const int &jj) const
    {return local_cell_node_vol_id[ii][jj];}

    // ------------------------------------------------------------------------
    // ! get the local cell node's location in the local_to_global array.
    //   \para 0 <= ii < num_ebc
    //   \para 0 <= jj < num_local_cell_node[ii]
    // ------------------------------------------------------------------------
    virtual int get_local_cell_node_pos(const int &ii, const int &jj) const
    {return local_cell_node_pos[ii][jj];}

    // ------------------------------------------------------------------------
    // ! get the local cell's volumetric mesh index.
    //   \para 0 <= ii < num_ebc
    //   \para 0 <= jj < num_local_cell[ii]
    // ------------------------------------------------------------------------
    virtual int get_local_cell_vol_id(const int &ii, const int &jj) const
    {return local_cell_vol_id[ii][jj];}

    // ------------------------------------------------------------------------
    // get_ctrlPts_xyz: given the ebc_id ii, the element index eindex,
    // return the control points' geometry.
    // Users are responsible for allocating & deleting the ctrl_xyz arrays.
    // ebc_id : 0 <= ii < num_ebc;
    // surface element id: 0 <= eindex < num_local_cell[ii];
    // ctrl_x/y/z : output geometry array, length is 
    //              cell_nLocBas[ii].
    // ------------------------------------------------------------------------
    virtual void get_ctrlPts_xyz(const int &ii, const int &eindex, 
        double * const &ctrl_x, double * const &ctrl_y, double * const &ctrl_z ) const;

    // ------------------------------------------------------------------------
    // get_SIEN: returns the surface element's IEN.
    // The users are responsible for allocating and deleting the sien  array.
    // ebc_id : 0 <= ii < num_ebc;
    // eindex : 0 <= eindex < num_local_cell[ii]
    // sien : length cell_nLocBas[ii].
    // ------------------------------------------------------------------------
    virtual void get_SIEN( const int &ii, const int &eindex, int * const &sien ) const;

    virtual std::vector<int> get_SIEN( const int &ii, const int &eindex ) const;

    // ------------------------------------------------------------------------
    // get_intPts_xyz: returns the surface element's interior point
    // coordinates.
    // ebc_id : 0 <= ii < num_ebc;
    // eindex : 0 <= eindex < num_local_cell[ii];
    // coor_x/y/z : output interior point coordinates
    // ------------------------------------------------------------------------
    virtual void get_intPts_xyz(const int &ii, const int &eindex, 
        double &coor_x, double &coor_y, double &coor_z ) const
    {
      SYS_T::print_fatal("Error: ALocal_EBC::get_intPts_xyz is not implemented. \n");
    }

    // ------------------------------------------------------------------------
    // Generate a file name for outlet face ii as Outlet_xx_flowrate.txt
    // ------------------------------------------------------------------------
    virtual std::string gen_flowfile_name(const int &ii) const
    {
      std::ostringstream ss;
      ss<<"Outlet_";
      if( ii/10 == 0 ) ss<<"00";
      else if(ii/100 == 0) ss<<"0";

      ss<<ii<<"_data.txt";
      
      return ss.str();
    }

    // ------------------------------------------------------------------------
    // The following functions are in the derived _outflow classes
    // for the Dirichlet-to-Neumann type outflow boundary conditions.
    // ------------------------------------------------------------------------
    // get_num_face_nodes : return the number of nodes on the surface
    //                      if this partition owns any cell on this
    //                      surface; otherwise return 0.
    //                      ii : face id ranging from 0 <= ii < num_ebc.
    // ------------------------------------------------------------------------
    virtual int get_num_face_nodes(const int &ii) const
    {
      SYS_T::print_fatal("Error: ALocal_EBC::get_num_face_nodes is not implemented. \n");
      return -1;
    }

    // ------------------------------------------------------------------------
    // get_intNA : returns the integral of N_A basis of the face,
    //             if this partition owns any cell on this surface.
    //             ii : face_id ranging 0 <= ii < num_ebc.
    //             output vector's length is get_num_face_nodes(ii).
    // ------------------------------------------------------------------------
    virtual std::vector<double> get_intNA( const int &ii ) const
    {
      SYS_T::print_fatal("Error: ALocal_EBC::get_intNA is not implemented. \n");
      return {};
    }
    
    // ------------------------------------------------------------------------
    // get_LID : returns the LID for the nodes associated with intNA,
    //           if this partition owns any cell on this surface.
    //           ii : face_id ranging from 0 <= ii < num_ebc,
    //           output is a vector, whose length is 3 x get_num_face_nodes(ii),
    //           and is the LID of node 0 x, y, z, node1 x, y, z, ...
    // ------------------------------------------------------------------------
    virtual std::vector<int> get_LID( const int &ii ) const
    {
      SYS_T::print_fatal("Error: ALocal_EBC::get_LID is not implemented. \n");
      return {};
    }

    // ------------------------------------------------------------------------
    // get_outvec : return the outward normal vector for face ii, if
    //              this partition owns any cell on this surface.
    //              ii : face_id ranging from 0 <= ii < num_ebc,
    // ------------------------------------------------------------------------
    virtual Vector_3 get_outvec( const int &ii ) const
    {
      SYS_T::print_fatal("Error: ALocal_EBC::get_outvec is not implemented. \n");
      return Vector_3();
    }

  protected:
    // the number of different ebc domain on which one may prescribe different
    // elemental boundary conditions.
    int num_ebc;

    // num_local_cell_node[ii] gives the ii-th ebc's local cell node number. These are
    //                    nodes associated with local cells.
    // num_local_cell[ii] gives the ii-th ebc's local cell number.
    // cell_nLocBas[ii]   gives the cell's number of node. e.g., 
    //                    triangle surface is 3,
    //                    quadralaterial surface is 4,
    //                    quadratic triangle surface is 6.
    std::vector<int> num_local_cell_node, num_local_cell, cell_nLocBas;

    // local_cell_node_xyz[ii] gives a list of local cell node's coordinates
    // size: num_ebc x (3 x num_local_cell_node[ii])
    std::vector< std::vector<double> > local_cell_node_xyz;

    // local_cell_ien[ii] gives the local cell's IEN array
    // size: num_ebc x (cell_nLocBas[ii] x num_local_cell[ii]) 
    std::vector< std::vector<int> > local_cell_ien;

    // local cell nodes' global indices
    // size: num_ebc x num_local_cell_node[ii]
    std::vector< std::vector<int> > local_cell_node_vol_id;

    // local cell node's position in the volumetric local portion's local_to_global array.
    // size: num_ebc x num_local_cell_node[ii]
    std::vector< std::vector<int> > local_cell_node_pos;

    // local cell's corresponding volumetric element indices
    // size: num_ebc x num_local_cell[ii]
    std::vector< std::vector<int> > local_cell_vol_id;
    
    ALocal_EBC() = delete; 
};

#endif
