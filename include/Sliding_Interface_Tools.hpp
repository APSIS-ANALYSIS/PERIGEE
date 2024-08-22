#ifndef SLIDING_INTERFACE_TOOLS
#define SLIDING_INTERFACE_TOOLS
// ============================================================================
// FE_Tools.hpp
// This file defines help functions for the sliding-interface technique.
//
// Date Created: Aug. 16 2024
// ============================================================================
#include "ALocal_Interface.hpp"
#include "FE_Tools.hpp"

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

      void get_rotated_local(const ALocal_Interface * const &itf,
        const int &ii, const int &tag, const int &ee,
        int * const &local_ien, double * const &local_sol) const
      {
        for(int nn = 0; nn < nLocBas; ++nn)
        {
          local_ien[nn] = itf->get_rotated_lien(ii, tag, ee * nLocBas + nn);

          for(int dd = 0; dd < dof_sol; ++dd)
            local_sol[dof_sol * nn + dd] = rotated_node_sol[ii][dof_sol * local_ien[nn] + dd];
        }
      }

      void update_node_sol(const PDNSolution * const &sol);

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

      void set_curr(const double &itf_id, const int &fixed_ee, const int &qua,
        const int &ele_tag, const int &rotated_ee, const std::vector<double> &xi);

      void get_curr(const double &itf_id, const int &fixed_ee, const int &qua,
        int &ele_tag, int &rotated_ee, std::vector<double> &xi) const;

      void search_all_opposite_point(
        const double &curr_time,
        FEAElement * const &fixed_elementv,
        FEAElement * const &rotated_elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        IQuadPts * const &free_quad,
        ALocal_Interface * const &itf_part );

      void search_opposite_point(
        const double &curr_time,
        const Vector_3 &fixed_pt,
        const ALocal_Interface * const &itf_part,
        const int &itf_id,
        FEAElement * rotated_elementv,
        FEAElement * elements,
        int &tag,
        int &rotated_ee,
        IQuadPts * const &rotated_xi );

      ~SI_quad_point() = default;

    private:
      int nqp_sur;

      // stores the current rotated element tag for each quadrature point
      // size: num_itf x num_fixed_ele[ii] x numQuadPts(surface)
      std::vector<std::vector<int>> curr_tag;

      // stores the current rotated element number for each quadrature point
      // size: num_itf x num_fixed_ele[ii] x numQuadPts(surface)
      std::vector<std::vector<int>> curr_ee;

      // stores the current rotated element xi = (r, s) for each quadrature point
      // size: num_itf x num_fixed_ele[ii] x numQuadPts(surface) x 2
      std::vector<std::vector<std::vector<double>>> curr_xi;

  };
}

#endif
