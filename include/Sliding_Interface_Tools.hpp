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
#include "FEAElementFactory.hpp"
#include "QuadPtsFactory.hpp"

namespace SI_T
{
  class SI_solution
  {
    public:
      SI_solution(const std::string &fileBaseName, const int &cpu_rank);

      ~SI_solution() = default;

      // Return the local ien array and the local solution array of a fixed layer element
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

      // Return the local solution array and mesh velocity of a rotated layer element
      // Used after get_rotated_mdisp
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

    private:
      const int cpu_rank;
      int nLocBas, dof_sol;

      // the number of the nodes from the fixed volume elements
      // size: num_itf
      std::vector<int> num_fixed_node;

      // stores the pressure and velocity info of the nodes from the fixed volume elements
      // size: num_itf x (4 x num_fixed_node[ii])
      std::vector<std::vector<double>> fixed_node_sol;
    
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
      SI_quad_point(const ALocal_Interface * const &itf, const FEType &in_type, const int &in_nqp_sur);

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
          const ALocal_Interface * const &itf_part,
          const SI_solution * const &SI_sol );

      void search_opposite_rotated_point(
          const Vector_3 &fixed_pt,
          const ALocal_Interface * const &itf_part,
          const SI_solution * const &SI_sol,
          const int &itf_id,
          int &tag,
          int &rotated_ee );

      void search_opposite_fixed_point(
          const Vector_3 &rotated_pt,
          const ALocal_Interface * const &itf_part,
          const SI_solution * const &SI_sol,
          const int &itf_id,
          int &tag,
          int &fixed_ee );

    private:
      const int nqp_sur;

      // For the interface integral
      const std::unique_ptr<FEAElement> anchor_elementv;

      // Defined with only one quadrature point,
      // given by the found closest point
      const std::unique_ptr<FEAElement> opposite_elementv;

      const std::unique_ptr<FEAElement> elements;

      // For the interface integral
      const std::unique_ptr<const IQuadPts> quad_s;

      // Free parametric surface quadrature point
      std::unique_ptr<IQuadPts> free_quad;

      // stores the current rotated element number for each quadrature point
      // size: num_itf x num_fixed_ele[ii] x numQuadPts(surface)
      std::vector<std::vector<int>> fixed_qp_curr_rotated_ee;

      // stores the current rotated element xi = (r, s) for each quadrature point
      // size: num_itf x num_fixed_ele[ii] x numQuadPts(surface) x 2
      std::vector<std::vector<double>> fixed_qp_curr_rotated_xi;
      std::vector<std::vector<double>> fixed_qp_curr_rotated_eta;

      // rotated counterpart
      std::vector<std::vector<int>> rotated_qp_curr_fixed_ee;
      std::vector<std::vector<double>> rotated_qp_curr_fixed_xi;
      std::vector<std::vector<double>> rotated_qp_curr_fixed_eta;
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
