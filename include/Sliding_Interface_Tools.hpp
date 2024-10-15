#ifndef SLIDING_INTERFACE_TOOLS
#define SLIDING_INTERFACE_TOOLS
// ============================================================================
// Sliding_Interface_Tools.hpp
// This file defines help functions for the sliding-interface technique.
//
// Date Created: Aug. 16 2024
// ============================================================================
#include "ALocal_Interface.hpp"
#include "IPLocAssem.hpp"
#include "FE_Tools.hpp"
#include "PETSc_Tools.hpp"

namespace SI_T
{
  class SI_solution
  {
    public:
      SI_solution(const std::string &fileBaseName, const int &cpu_rank);

      ~SI_solution() = default;

      void get_fixed_local(const ALocal_Interface * const &itf,
        const int &ii, const int &ee,
        int * const &local_ien, double * const &local_sol) const
      {
        for(int nn = 0; nn < nLocBas; ++nn)
        {
          local_ien[nn] = itf->get_fixed_lien(ii, ee * nLocBas + nn);

          for(int dd = 0; dd < dof_sol; ++dd)
            local_sol[dof_sol * nn + dd] = fixed_node_sol[ii][dof_sol * local_ien[nn] + dd];
        }
      }

      void get_projected_rotated_local(const int &ii, const int * const &fixed_local_ien,
        double * const &local_sol, double * const &local_sol_x, double * const &local_sol_y,
        double * const &local_sol_z, double * const &local_mvleo) const
      {
        for(int nn=0; nn < nLocBas; ++nn)
        {
          for(int dd=0; dd < dof_sol; ++dd)
          {
            local_sol[dof_sol * nn + dd] = f_r_node_sol[ii][dof_sol * fixed_local_ien[nn] + dd];
          }

          for(int dd=0; dd < 3; ++dd)
          {
            local_sol_x[3 * nn + dd] = f_r_node_sol_x[ii][3 * fixed_local_ien[nn] + dd];
            local_sol_y[3 * nn + dd] = f_r_node_sol_y[ii][3 * fixed_local_ien[nn] + dd];
            local_sol_z[3 * nn + dd] = f_r_node_sol_z[ii][3 * fixed_local_ien[nn] + dd];
            local_mvleo[3 * nn + dd] = f_r_node_mvelo[ii][3 * fixed_local_ien[nn] + dd];
          }
        }
      }

      // Return the local ien array and the local mdisp array of a rotated layer element
      void get_rotated_mdisp(const ALocal_Interface * const &itf,
        const int &ii, const int &ee,
        int * const &local_ien, double * const &local_disp) const
      {
        for(int nn = 0; nn < nLocBas; ++nn)
        {
          local_ien[nn] = itf->get_rotated_lien(ii, ee * nLocBas + nn);

          for(int dd = 0; dd < 3; ++dd)
            local_disp[3 * nn + dd] = rotated_node_mdisp[ii][3 * local_ien[nn] + dd];
        }
      }

      void get_rotated_local(const ALocal_Interface * const &itf,
        const int &ii, const int &ee,
        int * const &local_ien, double * const &local_sol, double * const &local_mvleo) const
      {
        for(int nn = 0; nn < nLocBas; ++nn)
        {
          local_ien[nn] = itf->get_rotated_lien(ii, ee * nLocBas + nn);

          for(int dd = 0; dd < dof_sol; ++dd)
          {
            local_sol[dof_sol * nn + dd] = rotated_node_sol[ii][dof_sol * local_ien[nn] + dd];
          }
          for(int dd = 0; dd < 3; ++dd)
          {
            local_mvleo[3 * nn + dd] = rotated_node_mvelo[ii][3 * local_ien[nn] + dd];
          }
        }
      }

      void get_rotated_local(const int &ii, const int * const &local_ien,
        double * const &local_sol, double * const &local_mvleo) const
      {
        for(int nn = 0; nn < nLocBas; ++nn)
        {
          for(int dd = 0; dd < dof_sol; ++dd)
          {
            local_sol[dof_sol * nn + dd] = rotated_node_sol[ii][dof_sol * local_ien[nn] + dd];
          }
          for(int dd = 0; dd < 3; ++dd)
          {
            local_mvleo[3 * nn + dd] = rotated_node_mvelo[ii][3 * local_ien[nn] + dd];
          }
        }
      }

      void get_projected_fixed_local(const int &ii, const int * const &rotated_local_ien,
        double * const &local_sol, double * const &local_sol_x, double * const &local_sol_y,
        double * const &local_sol_z) const
      {
        for(int nn=0; nn < nLocBas; ++nn)
        {
          for(int dd=0; dd < dof_sol; ++dd)
          {
            local_sol[dof_sol * nn + dd] = r_f_node_sol[ii][dof_sol * rotated_local_ien[nn] + dd];
          }
          for(int dd=0; dd < 3; ++dd)
          {
            local_sol_x[3 * nn + dd] = r_f_node_sol_x[ii][3 * rotated_local_ien[nn] + dd];
            local_sol_y[3 * nn + dd] = r_f_node_sol_y[ii][3 * rotated_local_ien[nn] + dd];
            local_sol_z[3 * nn + dd] = r_f_node_sol_z[ii][3 * rotated_local_ien[nn] + dd];
          }
        }
      }

      void update_node_sol(const PDNSolution * const &sol);

      void update_node_mvelo(const PDNSolution * const &mvelo);

      void update_node_mdisp(const PDNSolution * const &mdisp);

      void zero_node_sol()
      {
        for(int ii = 0; ii < VEC_T::get_size(num_fixed_node); ++ii)
        {
          for(int jj = 0; jj < dof_sol * num_fixed_node[ii]; ++jj)
            fixed_node_sol[ii][jj] = 0.0;

          for(int jj = 0; jj < dof_sol * num_rotated_node[ii]; ++jj)
            rotated_node_sol[ii][jj] = 0.0;
        }
      }

      void Zero_node_mvelo()
      {
        for(int ii = 0; ii < VEC_T::get_size(num_fixed_node); ++ii)
        {
          for(int jj = 0; jj < 3 * num_rotated_node[ii]; ++jj)
            rotated_node_mvelo[ii][jj] = 0.0;
        }
      }

      void Zero_node_disp()
      {
        for(int ii = 0; ii < VEC_T::get_size(num_fixed_node); ++ii)
        {
          for(int jj = 0; jj < 3 * num_rotated_node[ii]; ++jj)
            rotated_node_mdisp[ii][jj] = 0.0;
        }
      }

      int get_num_all_fixed_node() const
      {
        return num_all_fixed_node;
      }

      int get_L2_proj_fixed_node_pos(const int &ii, const int &jj) const
      {
        return L2_proj_fixed_node_pos[ii][jj];
      }

      int get_num_all_fixed_inner_node() const
      {
        return num_all_fixed_inner_node;
      }

      int get_fixed_inner_node(const int &ii) const
      {
        return all_fixed_inner_node[ii];
      }

      int get_num_all_rotated_node() const
      {
        return num_all_rotated_node;
      }

      int get_L2_proj_rotated_node_pos(const int &ii, const int &jj) const
      {
        return L2_proj_rotated_node_pos[ii][jj];
      }

      int get_num_all_rotated_inner_node() const
      {
        return num_all_rotated_inner_node;
      }

      int get_rotated_inner_node(const int &ii) const
      {
        return all_rotated_inner_node[ii];
      }

      void update_L2_proj_result(const Vec &L2_proj_lhs, const int &which);
      // which: 0 -- f_r_node_sol
      //        1,2,3 -- f_r_node_sol_x/y/z
      //        4 -- f_r_node_mvelo
      //        5 -- r_f_node_sol
      //        6,7,8 -- r_f_node_sol_x/y/z

    private:
      const int cpu_rank;
      int nLocBas, dof_sol;

      // the number of the nodes from the fixed volume elements
      // size: num_itf
      std::vector<int> num_fixed_node;

      // stores the pressure and velocity info of the nodes from the fixed volume elements
      // size: num_itf x (4 x num_fixed_node[ii])
      std::vector<std::vector<double>> fixed_node_sol;

      std::vector<std::vector<double>> f_r_node_sol;
      std::vector<std::vector<double>> f_r_node_sol_x;
      std::vector<std::vector<double>> f_r_node_sol_y;
      std::vector<std::vector<double>> f_r_node_sol_z;
      
      std::vector<std::vector<double>> f_r_node_mvelo;

      std::vector<std::vector<int>> L2_proj_fixed_node_pos;

      int num_all_fixed_node;

      int num_all_fixed_inner_node;

      std::vector<int> all_fixed_inner_node;
    
      // stores the partition tag of the nodes from the rotated volume elements
      // size: num_itf x num_fixed_node[ii]
      std::vector<std::vector<int>> fixed_node_part_tag;

      // stores the local position in the partition of the nodes from the rotated volume elements
      // size: num_itf x num_fixed_node[ii]
      std::vector<std::vector<int>> fixed_node_loc_pos;

      // the number of the nodes from the rotated volume elements
      // size: num_itf
      std::vector<int> num_rotated_node;

      // stores the pressure and velocity info of the nodes from the rotated volume elements
      // size: num_itf x (4 x num_rotated_node[ii])
      std::vector<std::vector<double>> rotated_node_sol;

      std::vector<std::vector<double>> r_f_node_sol;
      std::vector<std::vector<double>> r_f_node_sol_x;
      std::vector<std::vector<double>> r_f_node_sol_y;
      std::vector<std::vector<double>> r_f_node_sol_z;

      // r_f_node_mvelo are zeros 

      std::vector<std::vector<int>> L2_proj_rotated_node_pos;

      int num_all_rotated_node;

      int num_all_rotated_inner_node;

      std::vector<int> all_rotated_inner_node;

      // stores the mesh velocity info of the nodes from the rotated volume elements
      // size: num_itf x (3 x num_rotated_node[ii])    
      std::vector<std::vector<double>> rotated_node_mvelo;    

      // stores the mesh displacement info of the nodes from the rotated volume elements
      // size: num_itf x (3 x num_rotated_node[ii])    
      std::vector<std::vector<double>> rotated_node_mdisp;      

      // stores the partition tag of the nodes from the rotated volume elements
      // size: num_itf x num_rotated_node[ii]
      std::vector<std::vector<int>> rotated_node_part_tag;

      // stores the local position in the partition of the nodes from the rotated volume elements
      // size: num_itf x num_rotated_node[ii]
      std::vector<std::vector<int>> rotated_node_loc_pos;

  };

  class SI_quad_point
  {
    public:
      SI_quad_point(const ALocal_Interface * const &itf, const int &nqp_sur_in);

      ~SI_quad_point() = default;

      void set_curr_rotated(const double &itf_id, const int &fixed_ee_index, const int &qua,
          const int &rotated_ee, const double &xi, const double &eta);

      void get_curr_rotated(const double &itf_id, const int &fixed_ee_index, const int &qua,
          int &rotated_ee, double &xi, double &eta) const;

      void set_curr_fixed(const double &itf_id, const int &rotated_ee_index, const int &qua,
          const int &fixed_ee, const double &xi, const double &eta);

      void get_curr_fixed(const double &itf_id, const int &rotated_ee_index, const int &qua,
          int &fixed_ee, double &xi, double &eta) const;

      void search_all_opposite_point(
          FEAElement * const &anchor_elementv,
          FEAElement * const &opposite_elementv,
          FEAElement * const &elements,
          const IQuadPts * const &quad_s,
          IQuadPts * const &free_quad,
          const ALocal_Interface * const &itf_part,
          const SI_solution * const &SI_sol );

      void search_opposite_rotated_point(
          const Vector_3 &fixed_pt,
          const ALocal_Interface * const &itf_part,
          const SI_solution * const &SI_sol,
          const int &itf_id,
          FEAElement * rotated_elementv,
          FEAElement * elements,
          int &tag,
          int &rotated_ee,
          IQuadPts * const &rotated_xi );

      void search_opposite_fixed_point(
          const Vector_3 &rotated_pt,
          const ALocal_Interface * const &itf_part,
          const SI_solution * const &SI_sol,
          const int &itf_id,
          FEAElement * fixed_elementv,
          FEAElement * elements,
          int &tag,
          int &fixed_ee,
          IQuadPts * const &fixed_xi );

    private:
      const int nqp_sur;

      // // stores the current rotated element tag for each quadrature point
      // // size: num_itf x num_fixed_ele[ii] x numQuadPts(surface)
      // std::vector<std::vector<int>> fixed_qp_curr_rotated_tag;

      // stores the current rotated element number for each quadrature point
      // size: num_itf x num_fixed_ele[ii] x numQuadPts(surface)
      std::vector<std::vector<int>> fixed_qp_curr_rotated_ee;

      // stores the current rotated element xi = (r, s) for each quadrature point
      // size: num_itf x num_fixed_ele[ii] x numQuadPts(surface) x 2
      std::vector<std::vector<double>> fixed_qp_curr_rotated_xi;
      std::vector<std::vector<double>> fixed_qp_curr_rotated_eta;

      // std::vector<std::vector<int>> rotated_qp_curr_fixed_tag;
      std::vector<std::vector<int>> rotated_qp_curr_fixed_ee;
      std::vector<std::vector<double>> rotated_qp_curr_fixed_xi;
      std::vector<std::vector<double>> rotated_qp_curr_fixed_eta;
  };

  // ancillary parameters for PGAssem_NS_FEM
  class SI_ancillary
  {
    public:
      SI_ancillary(
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &fixed_elementv,
        FEAElement * const &rotated_elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        IQuadPts * const &free_quad,
        const ALocal_Interface * const &itf_part,
        SI_T::SI_solution * const &SI_sol,
        SI_T::SI_quad_point * const &SI_qp ):
        A_quad_s{quad_s}, A_itf_part{itf_part}
      {
        A_dt = 0.0;
        A_lassemptr = lassem_ptr;
        A_fixed_elementv = fixed_elementv;
        A_rotated_elementv = rotated_elementv;
        A_elements = elements;
        A_quad_s = quad_s;
        A_free_quad = free_quad;
        A_itf_part = itf_part;
        A_SI_sol = SI_sol;
        A_SI_qp = SI_qp;
      }

      ~SI_ancillary()
      {
        A_dt = 0.0;
        A_lassemptr = nullptr;
        A_fixed_elementv = nullptr;
        A_rotated_elementv = nullptr;
        A_elements = nullptr;
        A_free_quad = nullptr;
        A_SI_sol = nullptr;
        A_SI_qp = nullptr;
      }

      double A_dt;

      IPLocAssem * A_lassemptr;

      FEAElement * A_fixed_elementv;

      FEAElement * A_rotated_elementv;

      FEAElement * A_elements;

      const IQuadPts * A_quad_s;

      IQuadPts * A_free_quad;

      const ALocal_Interface * A_itf_part;

      SI_T::SI_solution * A_SI_sol;

      SI_T::SI_quad_point * A_SI_qp;

  };

  void get_currPts( const double * const &ept_x,
    const double * const &ept_y,
    const double * const &ept_z,
    const double * const &disp,
    const int &len,
    double * const &currPt_x,
    double * const &currPt_y,
    double * const &currPt_z );
}

#endif
