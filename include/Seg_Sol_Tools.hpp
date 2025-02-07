#ifndef SEG_SOL_TOOLS_HPP
#define SEG_SOL_TOOLS_HPP
// ==================================================================
// Seg_Sol_Tools.hpp
// ------------------------------------------------------------------
// SEG_SOL_T namespace contains a suite of frequently used tools for
// the manipulation of solution vectors in the segregated algorithms.
// 
// Date Created: May 21 2017
// ==================================================================
#include "PDNSolution.hpp"
#include "Math_Tools.hpp"
#include "ALocal_InflowBC.hpp"

namespace SEG_SOL_T
{
  // ===============================================================
  // Part 1: Update solution values
  // ----------------------------------------------------------------
  // ! UpdateU: The solution vector sol is assumed to have 7 dofs
  //            [disp_x, disp_y, disp_z, pres, vel_x, vel_y, vel_z].
  //            The solution vector usol is assumed to have 3 dofs
  //            [dd_x, dd_y, dd_z].
  //            The updateU will update the disp components in sol by
  //            the values of the usol:
  //            disp_x += val * dd_x
  //            disp_y += val * dd_y
  //            disp_z += val * dd_z.
  // ----------------------------------------------------------------
  void UpdateU( const double &val, const PDNSolution * const &usol,
      PDNSolution * const &sol );

  // ----------------------------------------------------------------  
  // ! UpdateV: The mesh displacement at time step n+1 is used to get
  //            the mesh velocity using the Generalized-alpha kinematic
  //            relation. The mesh velocity is stored in the first
  //            three slots in the dot_sol vector.
  //     mesh_velo_n+1 = (mesh_disp_n+1 - mesh_disp_n) / (gamma * dt)
  //                    + (1 - 1 / gamma) * mesh_velo_n
  // ----------------------------------------------------------------  
  void UpdateV( const double &dt, const double &gamma,
      const PDNSolution * const &pre_dot_sol,
      const PDNSolution * const &pre_sol,
      const PDNSolution * const &sol,
      PDNSolution * const &dot_sol );

  // ----------------------------------------------------------------
  // ! PlusAiPV: The solution vector sol is assumed to have 7 dofs
  //               [disp_x, disp_y, disp_z, pres, vel_x, vel_y, vel_z].
  //             The PV vector is assumed to have 4 dofs:
  //               [dp, dv_x, dv_y, dv_z].
  //             The updated solution takes the results:
  //                        disp_x += a dv_x
  //                        disp_y += a dv_y
  //                        disp_z += a dv_z
  //                        pres   += b dp
  //                        vel_x  += c dv_x
  //                        vel_y  += c dv_y
  //                        vel_z  += c dv_z.
  //             E.G., taking a = 0, b = -1, c = -1, one gets the
  //             update in the initialization of dot_sol vector using
  //             consistent mass matrix for the Generalized-alpha.
  // ----------------------------------------------------------------
  void PlusAiPV(const double &aa, const double &bb, const double &cc,
      const PDNSolution * const &step, PDNSolution * const &sol ); 


  // ----------------------------------------------------------------
  // ! PlusAiPV: The update is performed for the nodal indices given
  //             by pnode->node_loc_solid. The operation is identical
  //             to the original PlusAiPV.
  // ----------------------------------------------------------------
  void PlusAiPV(const double &aa, const double &bb, const double &cc,
      const APart_Node * const &pnode,
      const PDNSolution * const &step, PDNSolution * const &sol ); 


  // ----------------------------------------------------------------
  // ! PlusAiUPV: The solution vector sol is assumed to have 7 dofs
  //               [disp_x, disp_y, disp_z, pres, vel_x, vel_y, vel_z].
  //             The UPV vector is assumed to have 7 dofs:
  //               [dd_x, dd_y, dd_z, dp, dv_x, dv_y, dv_z].
  //             The updated solution takes the results:
  //                        disp_x += a dd_x
  //                        disp_y += a dd_y
  //                        disp_z += a dd_z
  //                        pres   += b dp
  //                        vel_x  += c dv_x
  //                        vel_y  += c dv_y
  //                        vel_z  += c dv_z.
  // ----------------------------------------------------------------
  void PlusAiUPV( const double &aa, const double &bb, const double &cc,
      const PDNSolution * const &step, PDNSolution * const &sol );


  // ----------------------------------------------------------------
  // ! PlusAiUPV: The solution is updated for the node given in the 
  //              list by node_loc_solid. The operation is identical 
  //              to the original PlusAiUPV function.
  // ----------------------------------------------------------------
  void PlusAiUPV(const double &aa, const double &bb, const double &cc,
      const APart_Node * const &pnode,
      const PDNSolution * const &step, PDNSolution * const &sol ); 


  // ----------------------------------------------------------------
  // ! PlusAiVPV: The solution vector sol is assumed to have 7 dofs
  //               [disp_x, disp_y, disp_z, pres, vel_x, vel_y, vel_z].
  //             The UPV vector is assumed to have 7 dofs:
  //               [dd_x, dd_y, dd_z, dp, dv_x, dv_y, dv_z].
  //             The updated solution takes the results:
  //                        disp_x += a dv_x
  //                        disp_y += a dv_y
  //                        disp_z += a dv_z
  //                        pres   += b dp
  //                        vel_x  += c dv_x
  //                        vel_y  += c dv_y
  //                        vel_z  += c dv_z.
  // ----------------------------------------------------------------
  void PlusAiVPV( const double &aa, const double &bb, const double &cc,
      const PDNSolution * const &step, PDNSolution * const &sol );
  
 
  // ---------------------------------------------------------------- 
  // ! PlusAiVPV: The solution is updated for the node givne in the
  //              list of node_loc_solid. The operatioin is identical
  //              to the original PlusAiVPV function.
  // ---------------------------------------------------------------- 
  void PlusAiVPV( const double &aa, const double &bb, const double &cc,
      const APart_Node * const &pnode,
      const PDNSolution * const &step, PDNSolution * const &sol );

  // ================================================================
  // Part 2: Insert values
  // ----------------------------------------------------------------
  // ! Insert_plug_inflow_UPV: Ths U_P_V 7-dof solution vector's Inflow 
  //                      nodes gets the inflow values inserted.
  //   This one insert plug inflow profile, at the inflow nodes,
  //          vel_x = val * out_normal_x       
  //          vel_y = val * out_normal_y       
  //          vel_z = val * out_normal_z       
  //   \para val: the inserted speed, i.e., velocity magnitude
  //   \para infnbc: the inflow nodal bc class
  //   \para sol : the solution class
  // ----------------------------------------------------------------
  void Insert_plug_inflow_UPV(const double &val,
     const ALocal_InflowBC * const &infnbc,
     PDNSolution * const &sol );

  void Insert_zero_solid_UV( const APart_Node * const &pnode,
      PDNSolution * const &sol );

  // ================================================================
  // Part 3: Extract values
  // ----------------------------------------------------------------
  // ! extract_solid_U : The U-P-V 7-dof solution's solid displacement
  //                     gets extracted to a 3-dof solution vector.
  //   This function is designed for assessing the convergence in the 
  //   prestress generation algorithm.
  // ================================================================
  void Extract_solid_U( const APart_Node * const &pnode,
      const PDNSolution * const &sol, PDNSolution * const &output );

  // ----------------------------------------------------------------
  // ! extract_solid_P : The U-P-V 7-dof solution's solid pressure
  //                     gets extracted to a single dof solution vector.
  //   Note: the pressure on the solid-fluid interface is set to be
  //   zero in the output. This function is designed for assessing the
  //   convergence in the prestress generation algorithm.
  // ----------------------------------------------------------------
  void Extract_solid_P( const APart_Node * const &pnode,
      const PDNSolution * const &sol, PDNSolution * const &output );

  // ================================================================
  // Part 4: Check solution
  // ----------------------------------------------------------------
  // ! CheckUV: This routine checks the dot{U}_n+alpha_m - V_n+alpha_f
  //            = 0 for all the consistent Newton-Raphson iterations.
  //            This means that, using the update procedure proposed
  //            in the segregated algorithm, the kinematic equation
  //            is guaranteed to be satisfied.
  //            If the equation is not satisfied up to a tolerance,
  //            an error message will be sent.
  // ----------------------------------------------------------------
  void CheckUV(const PDNSolution * const &dotU, 
      const PDNSolution * const &V );

  // ----------------------------------------------------------------
  // ! Check_dotSol_Sol: This routine checks the satisfaction of the
  //            relation
  //    y_n+1 = y_n + dt doty_n + gamma dt (doty_n+1 - doty_n)
  // ----------------------------------------------------------------
  void Check_dotSol_Sol( const double &dt, const double &gamma,
      const PDNSolution * const &dotSol_n,
      const PDNSolution * const &dotSol_np1,
      const PDNSolution * const &Sol_n,
      const PDNSolution * const &Sol_np1 );

}

#endif
