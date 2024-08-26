#ifndef ALOCAL_INTERFACE_HPP
#define ALOCAL_INTERFACE_HPP
// ============================================================================
// ALocal_Interface.hpp
//
// FEM-analysis-use Local subdomain's interface.
// Our model may contains more than one interface, therefore this class is 
// similar with ALocal_EBC, but it uses the volume elements' info.
//
// The volume elements attached to an interface would be from the fixed volume
// and the rotated volume. We distribute the fixed volume elements to each
// partition, but store the rotated volume elements' information to every
// partition h5 file which involves the interface. 

// When calculate the integral on the interface, we go through the fixed elements,
// use the quadrature points of the fixed face, and search for the opposite points
// on the rotated face.
//
// Author: Xuanming Huang
// Date Created: Jun. 24  2024
// ============================================================================
#include "PDNSolution.hpp"
#include "PETSc_Tools.hpp"
#include "Math_Tools.hpp"

class ALocal_Interface
{
  public:
    ALocal_Interface( const std::string &fileBaseName, const int &cpu_rank);

    virtual ~ALocal_Interface() = default;

    virtual void print_info() const;

    virtual int get_num_itf() const
    {return num_itf;}

    virtual int get_num_fixed_ele(const int &ii) const
    {return num_fixed_ele[ii];}

    virtual int get_num_tag(const int &ii) const
    {return num_tag[ii];}

    virtual int get_num_rotated_ele(const int &ii, const int &tag) const
    {return num_rotated_ele[ii][tag];}

    virtual int get_num_rotated_node(const int &ii) const
    {return num_rotated_node[ii];}

    virtual int get_fixed_face_id(const int &ii, const int &jj) const
    {return fixed_ele_face_id[ii][jj];}

    virtual int get_fixed_ele_tag(const int &ii, const int &jj) const
    {return fixed_ele_tag[ii][jj];}

    virtual int get_fixed_layer_ien(const int &ii, const int &jj) const
    {return fixed_layer_ien[ii][jj];}

    virtual double get_fixed_node_xyz(const int &ii, const int &jj) const
    {return fixed_node_xyz[ii][jj];}

    virtual int get_fixed_node_id(const int &ii, const int &jj) const
    {return fixed_node_id[ii][jj];}

    virtual int get_fixed_ID(const int &ii, const int &dof_index, const int &node) const
    {return fixed_ID[ii][dof_index * num_fixed_node[ii] + node];}

    // Return the local ien array and the local sol array of a fixed layer element
    virtual void get_fixed_local(const int &ii, const int &ee,
      int * const &local_ien, double * const &local_sol) const
    {
      for(int nn = 0; nn < nLocBas; ++nn)
      {
        const int node = get_fixed_layer_ien(ii, ee * nLocBas + nn);
        local_ien[nn] = node;

        for(int dd = 0; dd < dof_sol; ++dd)
          local_sol[dof_sol * nn + dd] = fixed_node_sol[ii][dof_sol * node + dd];
      }
    }

    virtual int get_rotated_layer_ien(const int &ii, const int &tag, const int &jj) const
    {return rotated_layer_ien[ii][tag][jj];}

    // Return the local ien array and the local sol array of a rotated layer element
    virtual void get_rotated_local(const int &ii, const int &tag, const int &ee,
      int * const &local_ien, double * const &local_sol) const
    {
      for(int nn = 0; nn < nLocBas; ++nn)
      {
        const int node = get_rotated_layer_ien(ii, tag, ee * nLocBas + nn);
        local_ien[nn] = node;

        for(int dd = 0; dd < dof_sol; ++dd)
          local_sol[dof_sol * nn + dd] = rotated_node_sol[ii][dof_sol * node + dd];
      }
    }

    // Return the local ien array and the local mvelo array of a rotated layer element
    virtual void get_rotated_mvelo(const int &ii, const int &tag, const int &ee,
      int * const &local_ien, double * const &local_mvleo) const
    {
      for(int nn = 0; nn < nLocBas; ++nn)
      {
        local_ien[nn] = get_rotated_layer_ien(ii, tag, ee * nLocBas + nn);

        for(int dd = 0; dd < 3; ++dd)
          {
            local_mvleo[3 * nn + dd] = rotated_node_mvelo[ii][3 * local_ien[nn] + dd];
          }
      }
    }

    // Return the local ien array and the local mdisp array of a rotated layer element
    virtual void get_rotated_disp(const int &ii, const int &tag, const int &ee,
      int * const &local_ien, double * const &local_disp) const
    {
      for(int nn = 0; nn < nLocBas; ++nn)
      {
        local_ien[nn] = get_rotated_layer_ien(ii, tag, ee * nLocBas + nn);

        for(int dd = 0; dd < 3; ++dd)
          {
            local_disp[3 * nn + dd] = rotated_node_disp[ii][3 * local_ien[nn] + dd];
          }
      }
    }

    virtual int get_rotated_face_id(const int &ii, const int &tag, const int &jj) const
    {return rotated_layer_face_id[ii][tag][jj];}

    virtual double get_init_rotated_node_xyz(const int &ii, const int &jj) const
    {return init_rotated_node_xyz[ii][jj];}

    virtual int get_rotated_node_id(const int &ii, const int &jj) const
    {return rotated_node_id[ii][jj];}

    virtual int get_rotated_ID(const int &ii, const int &dof_index, const int &node) const
    {return rotated_ID[ii][dof_index * num_rotated_node[ii] + node];}

    virtual void get_fixed_ele_ctrlPts(const int &ii, const int &ee,
      double * const volctrl_x,  double * const volctrl_y,  double * const volctrl_z) const;

    virtual void get_rotated_ele_ctrlPts(const int &ii, const int &tag,const int &ee, const double &tt,
      double * const volctrl_x,  double * const volctrl_y,  double * const volctrl_z) const;

    virtual void Zero_node_sol()
    {
      for(int ii = 0; ii < num_itf; ++ii)
      {
        for(int jj = 0; jj < dof_sol * num_fixed_node[ii]; ++jj)
          fixed_node_sol[ii][jj] = 0.0;
      }

      for(int ii = 0; ii < num_itf; ++ii)
      {
        for(int jj = 0; jj < dof_sol * num_rotated_node[ii]; ++jj)
          rotated_node_sol[ii][jj] = 0.0;
      }
    }

    virtual void Zero_node_mvelo()
    {
      for(int ii = 0; ii < num_itf; ++ii)
      {
        for(int jj = 0; jj < 3 * num_rotated_node[ii]; ++jj)
          rotated_node_mvelo[ii][jj] = 0.0;
      }
    }

    virtual void Zero_node_disp()
    {
      for(int ii = 0; ii < num_itf; ++ii)
      {
        for(int jj = 0; jj < 3 * num_rotated_node[ii]; ++jj)
          rotated_node_disp[ii][jj] = 0.0;
      }
    }

    virtual void restore_node_sol(const PDNSolution * const &sol);

    virtual void restore_node_mvelo(const PDNSolution * const &mvelo);

    virtual void restore_node_disp(const PDNSolution * const &disp);

    virtual void init_curr(const int &nqp_sur);

    virtual void set_curr(const double &itf_id, const int &fixed_ee, const int &qua,
      const int &ele_tag, const int &rotated_ee, const std::vector<double> &xi);

    virtual void get_curr(const double &itf_id, const int &fixed_ee, const int &qua,
      int &ele_tag, int &rotated_ee, std::vector<double> &xi) const;

  protected:
    // the number of interfaces
    int num_itf;

    // the number of local basis function
    int nLocBas;

    // the degree of freedom of the solution
    int dof_sol;

    // the number of local & ghost node in this part
    int nlgn;

    // cpu rank
    int cpu;

    int nqp_sur;

    // the number of fixed volume elements in this part
    // size: num_itf
    std::vector<int> num_fixed_ele;

    // the number of rotated volume elements of each interface
    // size: num_itf x num_tag[ii]
    std::vector<std::vector<int>> num_rotated_ele;

    // the number of the nodes from the fixed volume elements
    // size: num_itf
    std::vector<int> num_fixed_node;

    // the number of the nodes from the rotated volume elements
    // size: num_itf
    std::vector<int> num_rotated_node;

    // the number of interval tag of each pair of interfaces
    // size: num_itf
    std::vector<int> num_tag;

    // stores the face id of fixed volume element
    // size: num_itf x num_fixed_ele[ii]
    std::vector<std::vector<int>> fixed_ele_face_id;

    // stores the interval tag of fixed volume element
    // size: num_itf x num_fixed_ele[ii]
    std::vector<std::vector<int>> fixed_ele_tag;

    // stores the volume element's IEN array of the fixed "layer"
    // size: num_itf x (nlocbas x num_fixed_ele[ii])
    std::vector<std::vector<int>> fixed_layer_ien;

    // stores the initial coordinates of the nodes from the fixed volume elements
    // size: num_itf x (3 x num_fixed_node[ii])
    std::vector<std::vector<double>> fixed_node_xyz;

    // stores the pressure and velocity info of the nodes from the fixed volume elements
    // size: num_itf x (4 x num_fixed_node[ii])
    std::vector<std::vector<double>> fixed_node_sol;
    
    // stores the partition tag of the nodes from the rotated volume elements
    // size: num_itf x num_fixed_node[ii]
    std::vector<std::vector<int>> fixed_node_part_tag;

    // stores the local position in the partition of the nodes from the rotated volume elements
    // size: num_itf x num_fixed_node[ii]
    std::vector<std::vector<int>> fixed_node_loc_pos;

    // the (mapped) global node id corresponding to the fixed_node_xyz
    // size: num_itf x num_fixed_node[ii]
    std::vector<std::vector<int>> fixed_node_id;

    // ID array
    // size: num_itf x (num_fixed_node[ii] x dof)
    std::vector<std::vector<int>> fixed_ID;
    
    // stores the volume element's IEN array of the rotated "layer"
    // size: num_itf x (nlocbas x num_rotated_ele[ii])
    std::vector<std::vector<std::vector<int>>> rotated_layer_ien;

    // stores the face id of all the rotated volume element
    // size: num_itf x num_rotated_ele[ii]
    std::vector<std::vector<std::vector<int>>> rotated_layer_face_id;

    // stores the initial coordinates of the nodes from the rotated volume elements
    // size: num_itf x (3 x num_rotated_node[ii])
    std::vector<std::vector<double>> init_rotated_node_xyz;

    // stores the pressure and velocity info of the nodes from the rotated volume elements
    // size: num_itf x (4 x num_rotated_node[ii])
    std::vector<std::vector<double>> rotated_node_sol;

    // stores the mesh velocity info of the nodes from the rotated volume elements
    // size: num_itf x (3 x num_rotated_node[ii])    
    std::vector<std::vector<double>> rotated_node_mvelo;    

    // stores the mesh displacement info of the nodes from the rotated volume elements
    // size: num_itf x (3 x num_rotated_node[ii])    
    std::vector<std::vector<double>> rotated_node_disp;  

    // stores the partition tag of the nodes from the rotated volume elements
    // size: num_itf x num_rotated_node[ii]
    std::vector<std::vector<int>> rotated_node_part_tag;

    // stores the local position in the partition of the nodes from the rotated volume elements
    // size: num_itf x num_rotated_node[ii]
    std::vector<std::vector<int>> rotated_node_loc_pos;

    // the (mapped) global node id corresponding to the init_rotated_node_xyz
    // size: num_itf x num_rotated_node[ii]
    std::vector<std::vector<int>> rotated_node_id;

    // ID array
    // size: num_itf x (num_rotated_node[ii] x dof)
    std::vector<std::vector<int>> rotated_ID;

    // stores the current rotated element tag for each quadrature point
    // size: num_itf x num_fixed_ele[ii] x numQuadPts(surface)
    std::vector<std::vector<int>> curr_tag;

    // stores the current rotated element number for each quadrature point
    // size: num_itf x num_fixed_ele[ii] x numQuadPts(surface)
    std::vector<std::vector<int>> curr_ee;

    // stores the current rotated element xi = (r, s) for each quadrature point
    // size: num_itf x num_fixed_ele[ii] x numQuadPts(surface) x 2
    std::vector<std::vector<std::vector<double>>> curr_xi;

    ALocal_Interface() = delete;
};
#endif
