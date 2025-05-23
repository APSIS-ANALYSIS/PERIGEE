// ============================================================================
// post_PLA_tet4.cpp
//
// The quadrature routine for Power Loss Analysis.
// It calculates P_dir, P_prod, P_mod, P_diss.
// 
// This routine works only for linear tetrahedral element.
//
// Author: Xuanming Huang
// Date: May 22 2025
// ============================================================================
#include "HDF5_Reader.hpp"
#include "Tet_Tools.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Tet4.hpp"
#include "Tensor2_3D.hpp"

std::vector<int> ReadNodeMapping( const char * const &node_mapping_file,
    const char * const &mapping_type, const int &node_size );

std::vector<double> ReadPETSc_Vec( const std::string &solution_file_name,
    const std::vector<int> &nodemap,
    const int &vec_size, const int &in_dof );

int main( int argc, char * argv[] )
{
  std::string sol_bname("SOL_");

  // visualization sampling pattern over time
  int time_start = 0;
  int time_step = 1;
  int time_end = 1;

  constexpr int dof = 4;

  constexpr int v_nqp = 5; 

#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif
  
  SYS_T::print_fatal_if(SYS_T::get_MPI_size() != 1, "ERROR: post_PLA_tet4 is a serial program! \n");

  // Directly read in the volumetric and wall file from the file
  // that record the preprocessor command lines.
  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );
  const std::string geo_file  = cmd_h5r -> read_string("/", "geo_file");
  const std::string elemType_str = cmd_h5r -> read_string("/", "elemType");
  const FEType elemType = FE_T::to_FEType(elemType_str);

  delete cmd_h5r; H5Fclose(prepcmd_file);

  // Now read the material properties from the solver cmd h5 file
  prepcmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  
  cmd_h5r = new HDF5_Reader( prepcmd_file );
  
  const double rho = cmd_h5r -> read_doubleScalar("/", "fl_density");
  const double fluid_mu = cmd_h5r -> read_doubleScalar("/", "fl_mu");

  delete cmd_h5r; H5Fclose(prepcmd_file);

  // enforce this code is for linear element only
  SYS_T::print_fatal_if( elemType != FEType::Tet4, "Error: element type should be linear tet element.\n");

  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);  

  std::string out_bname = sol_bname;

  // Print the key data on screen
  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -sol_bname: "<<sol_bname<<endl;
  cout<<" -time_start: "<<time_start<<endl;
  cout<<" -time_step: "<<time_step<<endl;
  cout<<" -time_end: "<<time_end<<endl;
  cout<<"----------------------------------\n";
  cout<<" geo_file: "<<geo_file<<endl;
  cout<<" elemType: "<<elemType_str<<endl;
  cout<<" fl_density: "<<rho<<endl;
  cout<<" fl_mu: "<<fluid_mu<<endl;
  cout<<"==== Command Line Arguments ===="<<endl;

  // Make sure the files exist on disk
  SYS_T::file_check( geo_file );

  // Now read in the volumetric mesh info
  int v_nFunc, v_nElem;
  std::vector<int> v_vecIEN;
  std::vector<double> v_ctrlPts;

  VTK_T::read_vtu_grid(geo_file, v_nFunc, v_nElem, v_ctrlPts, v_vecIEN);

  cout<<endl<<"Volumetric mesh contains "<<v_nElem<<" elements and "<<v_nFunc<<" vertices.\n";

  const std::vector<int> global_node_idx = VTK_T::read_int_PointData(geo_file, "GlobalNodeID");
  const std::vector<int> global_ele_idx  = VTK_T::read_int_CellData(geo_file, "GlobalElementID");

  // Quadrature rule
  IQuadPts * quad = new QuadPts_Gauss_Tet( v_nqp );

  quad -> print_info();

  // Element container
  FEAElement * element = new FEAElement_Tet4( quad-> get_num_quadPts() );

  // Read the node mappings
  const auto analysis_new2old = ReadNodeMapping("node_mapping.h5", "new_2_old", v_nFunc );

  // Read solutions
  std::ostringstream time_index;

  const double inv_T = 1.0 / ( static_cast<double>((time_end - time_start)/time_step) + 1.0 );

  // Container for time-averaged velocity on nodes
  std::vector<Vector_3> ave_velo(v_nFunc, Vector_3(0, 0, 0));

  // Container for time-averaged velocity gradient on quadrature points
  std::vector<Tensor2_3D> ave_velo_grad(v_nqp * v_nElem, Tensor2_3D(0, 0, 0, 0, 0, 0, 0, 0, 0));

  // Container for time-averged Renold stress on quadrature points
  std::vector<Tensor2_3D> ave_Renold_stress(v_nqp * v_nElem, Tensor2_3D(0, 0, 0, 0, 0, 0, 0, 0, 0));

  // Container for P_dir, P_prod
  double P_dir = 0, P_prod = 0;

  // Container for P_mod, P_diss
  double P_mod = 0, P_diss = 0;

  // Calculate ave_velo and ave_velo_grad
  std::cout<<"==== Calculate ave_velo and ave_velo_grad ===="<<std::endl;
  for(int time = time_start; time <= time_end; time += time_step)
  {
    // Generate the file name
    std::string name_to_read(sol_bname);
    std::string name_to_write(out_bname);
    time_index.str("");
    time_index << 900000000 + time;
    name_to_read.append(time_index.str());

    std::cout<<"Time "<<time<<": Read "<<name_to_read<<std::endl;

    // Read in the solution vector and arrange them into the natural numbering
    const auto sol = ReadPETSc_Vec( name_to_read, analysis_new2old, v_nFunc*dof, dof );

    for(int nn=0; nn<v_nFunc; ++nn)
    {
      // Use nn to store, and use global_node_idx to get, like sol
      Vector_3 node_velo { sol[ nn * dof + 1 ],
                           sol[ nn * dof + 2 ],
                           sol[ nn * dof + 3 ] };
      ave_velo[nn] += inv_T * node_velo;
    }

    for(int ee=0; ee<v_nElem; ++ee)
    {
      const int tet_node[4] { v_vecIEN[4*ee+0], v_vecIEN[4*ee+1],
                              v_vecIEN[4*ee+2], v_vecIEN[4*ee+3] };

      const double ectrl_x[4] { v_ctrlPts[3*tet_node[0] + 0], v_ctrlPts[3*tet_node[1] + 0],
                                v_ctrlPts[3*tet_node[2] + 0], v_ctrlPts[3*tet_node[3] + 0] };

      const double ectrl_y[4] { v_ctrlPts[3*tet_node[0] + 1], v_ctrlPts[3*tet_node[1] + 1],
                                v_ctrlPts[3*tet_node[2] + 1], v_ctrlPts[3*tet_node[3] + 1] };

      const double ectrl_z[4] { v_ctrlPts[3*tet_node[0] + 2], v_ctrlPts[3*tet_node[1] + 2],
                                v_ctrlPts[3*tet_node[2] + 2], v_ctrlPts[3*tet_node[3] + 2] };

      element -> buildBasis(quad, ectrl_x, ectrl_y, ectrl_z);

      const double esol_u[4] { sol[ global_node_idx[tet_node[0]] * dof + 1 ],
                               sol[ global_node_idx[tet_node[1]] * dof + 1 ],
                               sol[ global_node_idx[tet_node[2]] * dof + 1 ],
                               sol[ global_node_idx[tet_node[3]] * dof + 1 ] };

      const double esol_v[4] { sol[ global_node_idx[tet_node[0]] * dof + 2 ],
                               sol[ global_node_idx[tet_node[1]] * dof + 2 ],
                               sol[ global_node_idx[tet_node[2]] * dof + 2 ],
                               sol[ global_node_idx[tet_node[3]] * dof + 2 ] };

      const double esol_w[4] { sol[ global_node_idx[tet_node[0]] * dof + 3 ],
                               sol[ global_node_idx[tet_node[1]] * dof + 3 ],
                               sol[ global_node_idx[tet_node[2]] * dof + 3 ],
                               sol[ global_node_idx[tet_node[3]] * dof + 3 ] };

      for(int qua=0; qua<v_nqp; ++qua)
      {

        double Rx[4], Ry[4], Rz[4];
        element -> get_gradR(qua, Rx, Ry, Rz);

        // Velocity gradient
        const double ux = esol_u[0] * Rx[0] + esol_u[1] * Rx[1] + esol_u[2] * Rx[2] + esol_u[3] * Rx[3];
        const double vx = esol_v[0] * Rx[0] + esol_v[1] * Rx[1] + esol_v[2] * Rx[2] + esol_v[3] * Rx[3];
        const double wx = esol_w[0] * Rx[0] + esol_w[1] * Rx[1] + esol_w[2] * Rx[2] + esol_w[3] * Rx[3];

        const double uy = esol_u[0] * Ry[0] + esol_u[1] * Ry[1] + esol_u[2] * Ry[2] + esol_u[3] * Ry[3];
        const double vy = esol_v[0] * Ry[0] + esol_v[1] * Ry[1] + esol_v[2] * Ry[2] + esol_v[3] * Ry[3];
        const double wy = esol_w[0] * Ry[0] + esol_w[1] * Ry[1] + esol_w[2] * Ry[2] + esol_w[3] * Ry[3];

        const double uz = esol_u[0] * Rz[0] + esol_u[1] * Rz[1] + esol_u[2] * Rz[2] + esol_u[3] * Rz[3];
        const double vz = esol_v[0] * Rz[0] + esol_v[1] * Rz[1] + esol_v[2] * Rz[2] + esol_v[3] * Rz[3];
        const double wz = esol_w[0] * Rz[0] + esol_w[1] * Rz[1] + esol_w[2] * Rz[2] + esol_w[3] * Rz[3];

        ave_velo_grad[v_nqp * ee + qua] += inv_T * Tensor2_3D{ ux, uy, uz,
                                                               vx, vy, vz,
                                                               wx, wy, wz };
      }
    }
  }// Loop over each time instance

  // Calculate ave_Renold_stress, P_mod, P_diss
  std::cout<<"==== Calculate ave_Renold_stress, P_dir, P_mod, P_diss ===="<<std::endl;
  for(int time = time_start; time <= time_end; time += time_step)
  {
    // Generate the file name
    std::string name_to_read(sol_bname);
    std::string name_to_write(out_bname);
    time_index.str("");
    time_index << 900000000 + time;
    name_to_read.append(time_index.str());

    std::cout<<"Time "<<time<<": Read "<<name_to_read<<std::endl;

    // Read in the solution vector and arrange them into the natural numbering
    const auto sol = ReadPETSc_Vec( name_to_read, analysis_new2old, v_nFunc*dof, dof );

    for(int ee=0; ee<v_nElem; ++ee)
    {
      const int tet_node[4] { v_vecIEN[4*ee+0], v_vecIEN[4*ee+1],
                              v_vecIEN[4*ee+2], v_vecIEN[4*ee+3] };

      const double ectrl_x[4] { v_ctrlPts[3*tet_node[0] + 0], v_ctrlPts[3*tet_node[1] + 0],
                                v_ctrlPts[3*tet_node[2] + 0], v_ctrlPts[3*tet_node[3] + 0] };

      const double ectrl_y[4] { v_ctrlPts[3*tet_node[0] + 1], v_ctrlPts[3*tet_node[1] + 1],
                                v_ctrlPts[3*tet_node[2] + 1], v_ctrlPts[3*tet_node[3] + 1] };

      const double ectrl_z[4] { v_ctrlPts[3*tet_node[0] + 2], v_ctrlPts[3*tet_node[1] + 2],
                                v_ctrlPts[3*tet_node[2] + 2], v_ctrlPts[3*tet_node[3] + 2] };

      element -> buildBasis(quad, ectrl_x, ectrl_y, ectrl_z);

      const double esol_u[4] { sol[ global_node_idx[tet_node[0]] * dof + 1 ],
                               sol[ global_node_idx[tet_node[1]] * dof + 1 ],
                               sol[ global_node_idx[tet_node[2]] * dof + 1 ],
                               sol[ global_node_idx[tet_node[3]] * dof + 1 ] };

      const double esol_v[4] { sol[ global_node_idx[tet_node[0]] * dof + 2 ],
                               sol[ global_node_idx[tet_node[1]] * dof + 2 ],
                               sol[ global_node_idx[tet_node[2]] * dof + 2 ],
                               sol[ global_node_idx[tet_node[3]] * dof + 2 ] };

      const double esol_w[4] { sol[ global_node_idx[tet_node[0]] * dof + 3 ],
                               sol[ global_node_idx[tet_node[1]] * dof + 3 ],
                               sol[ global_node_idx[tet_node[2]] * dof + 3 ],
                               sol[ global_node_idx[tet_node[3]] * dof + 3 ] };

      // Fluctuation on nodes
      const double fluc_u[4] { esol_u[0] - ave_velo[global_node_idx[tet_node[0]]].x(),
                               esol_u[1] - ave_velo[global_node_idx[tet_node[1]]].x(),
                               esol_u[2] - ave_velo[global_node_idx[tet_node[2]]].x(),
                               esol_u[3] - ave_velo[global_node_idx[tet_node[3]]].x() };

      const double fluc_v[4] { esol_v[0] - ave_velo[global_node_idx[tet_node[0]]].y(),
                               esol_v[1] - ave_velo[global_node_idx[tet_node[1]]].y(),
                               esol_v[2] - ave_velo[global_node_idx[tet_node[2]]].y(),
                               esol_v[3] - ave_velo[global_node_idx[tet_node[3]]].y() };

      const double fluc_w[4] { esol_w[0] - ave_velo[global_node_idx[tet_node[0]]].z(),
                               esol_w[1] - ave_velo[global_node_idx[tet_node[1]]].z(),
                               esol_w[2] - ave_velo[global_node_idx[tet_node[2]]].z(),
                               esol_w[3] - ave_velo[global_node_idx[tet_node[3]]].z() };
      
      std::array<Vector_3, 4> node_coor { Vector_3(ectrl_x[0], ectrl_y[0], ectrl_z[0]),
                                          Vector_3(ectrl_x[1], ectrl_y[1], ectrl_z[1]),
                                          Vector_3(ectrl_x[2], ectrl_y[2], ectrl_z[2]),
                                          Vector_3(ectrl_x[3], ectrl_y[3], ectrl_z[3]) };
      TET_T::Tet4 tet_elem (node_coor);

      const double e_delta = std::pow(tet_elem.get_volume(), 1.0/3.0);
      // Here we use a constant C_s = 0.17
      const double CC = std::pow(0.17 * e_delta, 2);

      for(int qua=0; qua<v_nqp; ++qua)
      {
        double R[4], Rx[4], Ry[4], Rz[4];
        element -> get_R_gradR(qua, R, Rx, Ry, Rz);

        // Velocity gradient
        const double ux = esol_u[0] * Rx[0] + esol_u[1] * Rx[1] + esol_u[2] * Rx[2] + esol_u[3] * Rx[3];
        const double vx = esol_v[0] * Rx[0] + esol_v[1] * Rx[1] + esol_v[2] * Rx[2] + esol_v[3] * Rx[3];
        const double wx = esol_w[0] * Rx[0] + esol_w[1] * Rx[1] + esol_w[2] * Rx[2] + esol_w[3] * Rx[3];

        const double uy = esol_u[0] * Ry[0] + esol_u[1] * Ry[1] + esol_u[2] * Ry[2] + esol_u[3] * Ry[3];
        const double vy = esol_v[0] * Ry[0] + esol_v[1] * Ry[1] + esol_v[2] * Ry[2] + esol_v[3] * Ry[3];
        const double wy = esol_w[0] * Ry[0] + esol_w[1] * Ry[1] + esol_w[2] * Ry[2] + esol_w[3] * Ry[3];

        const double uz = esol_u[0] * Rz[0] + esol_u[1] * Rz[1] + esol_u[2] * Rz[2] + esol_u[3] * Rz[3];
        const double vz = esol_v[0] * Rz[0] + esol_v[1] * Rz[1] + esol_v[2] * Rz[2] + esol_v[3] * Rz[3];
        const double wz = esol_w[0] * Rz[0] + esol_w[1] * Rz[1] + esol_w[2] * Rz[2] + esol_w[3] * Rz[3];

        // Velocity fluctuation
        const double fu = fluc_u[0] * R[0] + fluc_u[1] * R[1] + fluc_u[2] * R[2] + fluc_u[3] * R[3];
        const double fv = fluc_v[0] * R[0] + fluc_v[1] * R[1] + fluc_v[2] * R[2] + fluc_v[3] * R[3];
        const double fw = fluc_w[0] * R[0] + fluc_w[1] * R[1] + fluc_w[2] * R[2] + fluc_w[3] * R[3];

        // Fluctuation gradient
        const double fux = fluc_u[0] * Rx[0] + fluc_u[1] * Rx[1] + fluc_u[2] * Rx[2] + fluc_u[3] * Rx[3];
        const double fvx = fluc_v[0] * Rx[0] + fluc_v[1] * Rx[1] + fluc_v[2] * Rx[2] + fluc_v[3] * Rx[3];
        const double fwx = fluc_w[0] * Rx[0] + fluc_w[1] * Rx[1] + fluc_w[2] * Rx[2] + fluc_w[3] * Rx[3];

        const double fuy = fluc_u[0] * Ry[0] + fluc_u[1] * Ry[1] + fluc_u[2] * Ry[2] + fluc_u[3] * Ry[3];
        const double fvy = fluc_v[0] * Ry[0] + fluc_v[1] * Ry[1] + fluc_v[2] * Ry[2] + fluc_v[3] * Ry[3];
        const double fwy = fluc_w[0] * Ry[0] + fluc_w[1] * Ry[1] + fluc_w[2] * Ry[2] + fluc_w[3] * Ry[3];

        const double fuz = fluc_u[0] * Rz[0] + fluc_u[1] * Rz[1] + fluc_u[2] * Rz[2] + fluc_u[3] * Rz[3];
        const double fvz = fluc_v[0] * Rz[0] + fluc_v[1] * Rz[1] + fluc_v[2] * Rz[2] + fluc_v[3] * Rz[3];
        const double fwz = fluc_w[0] * Rz[0] + fluc_w[1] * Rz[1] + fluc_w[2] * Rz[2] + fluc_w[3] * Rz[3];

        ave_Renold_stress[v_nqp * ee +qua] -= inv_T * Tensor2_3D{ fu*fu, fu*fv, fu*fw,
                                                                  fv*fu, fv*fv, fv*fw,
                                                                  fw*fu, fw*fv, fw*fw };

        Tensor2_3D fg { fux, fuy, fuz, 
                        fvx, fvy, fvz,
                        fwx, fwy, fwz };

        Tensor2_3D vg { ux, uy, uz,
                        vx, vy, vz,
                        wx, wy, wz };

        Tensor2_3D avg = ave_velo_grad[v_nqp * ee + qua];
        Tensor2_3D avg_T = avg;
        avg_T.transpose();

        Tensor2_3D rate_of_strain = 0.5 * (avg + avg_T);
        double S_bar = std::sqrt(2 * (rate_of_strain.MatContraction()));

        double mu_t = CC * S_bar;

        double detJ = element -> get_detJac(qua);
        double qw = quad -> get_qw(qua);

        P_mod += inv_T * mu_t * detJ * qw * (vg.MatContraction() + vg.MatTContraction(vg));

        P_diss += inv_T * fluid_mu * detJ * qw * (fg.MatContraction() + fg.MatTContraction(fg));
      }
    }
  }

  // Calculate P_dir, P_prod
  std::cout<<"==== Calculate P_dir, P_prod ===="<<std::endl;
  for(int ee=0; ee<v_nElem; ++ee)
  {
    const int tet_node[4] { v_vecIEN[4*ee+0], v_vecIEN[4*ee+1],
                            v_vecIEN[4*ee+2], v_vecIEN[4*ee+3] };

    const double ectrl_x[4] { v_ctrlPts[3*tet_node[0] + 0], v_ctrlPts[3*tet_node[1] + 0],
                              v_ctrlPts[3*tet_node[2] + 0], v_ctrlPts[3*tet_node[3] + 0] };

    const double ectrl_y[4] { v_ctrlPts[3*tet_node[0] + 1], v_ctrlPts[3*tet_node[1] + 1],
                              v_ctrlPts[3*tet_node[2] + 1], v_ctrlPts[3*tet_node[3] + 1] };

    const double ectrl_z[4] { v_ctrlPts[3*tet_node[0] + 2], v_ctrlPts[3*tet_node[1] + 2],
                              v_ctrlPts[3*tet_node[2] + 2], v_ctrlPts[3*tet_node[3] + 2] };

    element -> buildBasis(quad, ectrl_x, ectrl_y, ectrl_z);

    for(int qua=0; qua<v_nqp; ++qua)
    {
      Tensor2_3D avg = ave_velo_grad[v_nqp * ee + qua];
      Tensor2_3D ars = ave_Renold_stress[v_nqp * ee + qua];

      double detJ = element -> get_detJac(qua);
      double qw = quad -> get_qw(qua);

      P_dir += fluid_mu * detJ * qw * (avg.MatContraction() + avg.MatTContraction(avg));

      P_prod += rho * detJ * qw * ars.MatContraction(avg);
    }
  }
  
  const double P_tot_2 = P_dir + P_prod + P_mod;
  const double P_tot_3 = P_dir + P_diss + P_mod;

  std::cout<<"==== Result ===="<<std::endl;
  std::cout<<"  P_dir = "<<P_dir<<std::endl;
  std::cout<<"  P_prod = "<<P_prod<<std::endl;
  std::cout<<"  P_mod = "<<P_mod<<std::endl;
  std::cout<<"  P_diss = "<<P_diss<<std::endl;
  std::cout<<"  P_tot_2 = P_dir + P_prod + P_mod = "<<P_tot_2<<std::endl;
  std::cout<<"  P_tot_3 = P_dir + P_diss + P_mod = "<<P_tot_3<<std::endl;

  delete quad; delete element;
  PetscFinalize();
  return EXIT_SUCCESS;
}
// END of MAIN Function


// Read the node mappings
std::vector<int> ReadNodeMapping( const char * const &node_mapping_file,
    const char * const &mapping_type, const int &node_size )
{
  hid_t file_id = H5Fopen(node_mapping_file, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t data_id = H5Dopen(file_id, mapping_type, H5P_DEFAULT);

  hid_t data_space = H5Dget_space( data_id );
  hid_t data_rank = H5Sget_simple_extent_ndims( data_space );

  if( data_rank != 1)
  {
    SYS_T::commPrint("Error: the node mapping file has wrong format. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  hsize_t * data_dims = new hsize_t [1];

  H5Sget_simple_extent_dims( data_space, data_dims, NULL );

  hid_t mem_space = H5Screate_simple(data_rank, data_dims, NULL);

  hsize_t dSize = data_dims[0];

  if( int(dSize) != node_size )
  {
    SYS_T::commPrint("Error: the allocated array has wrong size! \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  std::vector<int> out(node_size, -1);

  H5Dread( data_id, H5T_NATIVE_INT, mem_space, data_space, H5P_DEFAULT, &out[0] );

  delete [] data_dims;
  H5Sclose( mem_space );
  H5Sclose(data_space);
  H5Dclose(data_id);
  H5Fclose(file_id);

  return out;
}

std::vector<double> ReadPETSc_Vec( const std::string &solution_file_name,
    const std::vector<int> &nodemap,
    const int &vec_size, const int &in_dof )
{
  Vec sol_temp;
  VecCreate(PETSC_COMM_SELF, &sol_temp);
  VecSetType(sol_temp, VECSEQ);

  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_SELF, solution_file_name.c_str(),
      FILE_MODE_READ, &viewer);
  VecLoad(sol_temp, viewer);
  PetscViewerDestroy(&viewer);

  // Check the solution length
  PetscInt get_sol_temp_size;
  VecGetSize(sol_temp, &get_sol_temp_size);
  if( get_sol_temp_size != vec_size )
  {
    SYS_T::commPrint("The solution size %d is not compatible with the size %d given by partition file! \n",
        get_sol_temp_size, vec_size);
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  std::vector<double> veccopy(vec_size, 0.0);
  double * array_temp;
  VecGetArray(sol_temp, &array_temp);

  for(int ii=0; ii<vec_size; ++ii) veccopy[ii] = array_temp[ii];

  VecRestoreArray(sol_temp, &array_temp);
  VecDestroy(&sol_temp);

  // copy the solution varibles to the correct location
  std::vector<double> sol(vec_size, 0.0);

  // check the nodemap size
  if( (int)nodemap.size() * in_dof != vec_size ) SYS_T::print_fatal("Error: node map size is incompatible with the solution length. \n");

  for(unsigned int ii=0; ii<nodemap.size(); ++ii)
  {
    const int index = nodemap[ii];
    for(int jj=0; jj<in_dof; ++jj)
      sol[in_dof*index+jj] = veccopy[in_dof*ii+jj];
  }

  return sol;
}

// EOF