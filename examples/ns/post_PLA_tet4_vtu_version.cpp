// ============================================================================
// post_PLA_tet4.cpp
//
// The quadrature routine for Power Loss Analysis.
// It calculates P_dir, P_prod, P_mod, P_diss.
// 
// This routine works only for linear tetrahedral element.
//
// This versin reads the solution from vtu files
//
// Author: Xuanming Huang
// Date: May 22 2025
// ============================================================================
#include "Tet_Tools.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Tet4.hpp"
#include "Tensor2_3D.hpp"

std::string get_vtu_result_name(const std::string &pre_fix, const int &time,
  const int &is_pvtu = 0);

int main( int argc, char * argv[] )
{
  std::string sol_bname("50M_c2_");

  // visualization sampling pattern over time
  int time_start = 0;
  int time_step = 1;
  int time_end = 1;

  // constexpr int dof = 4;

  constexpr int v_nqp = 5;

  std::string geo_file;

  std::string elemType_str("Tet4");

  double rho = 1050.0;
  double fluid_mu = 0.0035;

  int is_pvtu = 0;

  // double aspect_ratio_threshold = 100;

#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif
  
  SYS_T::print_fatal_if(SYS_T::get_MPI_size() != 1, "ERROR: post_PLA_tet4_vtu_version is a serial program! \n");

  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionString("-elemType", elemType_str);
  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionInt("-is_pvtu", is_pvtu);
  SYS_T::GetOptionReal("-rho", rho);
  SYS_T::GetOptionReal("-fluid_mu", fluid_mu);

  // Print the key data on screen
  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -sol_bname: "<<sol_bname<<endl;
  cout<<" -time_start: "<<time_start<<endl;
  cout<<" -time_step: "<<time_step<<endl;
  cout<<" -time_end: "<<time_end<<endl;
  cout<<" -is_pvtu: "<<is_pvtu<<endl;
  cout<<"----------------------------------\n";
  cout<<" elemType: "<<elemType_str<<endl;
  cout<<" fl_density: "<<rho<<endl;
  cout<<" fl_mu: "<<fluid_mu<<endl;
  cout<<"==== Command Line Arguments ===="<<endl;

  const FEType elemType = FE_T::to_FEType(elemType_str);

  // enforce this code is for linear element only
  SYS_T::print_fatal_if( elemType != FEType::Tet4, "Error: element type should be linear tet element.\n");

  // #pragma omp parallel
  // {
  //     std::cout << "Hello from thread " << omp_get_thread_num() 
  //               << " of " << omp_get_num_threads() << std::endl;
  // }

  // Now read in the volumetric mesh info
  int v_nFunc, v_nElem;
  std::vector<int> v_vecIEN;
  std::vector<double> v_ctrlPts;

  // Quadrature rule
  IQuadPts * quad = new QuadPts_Gauss_Tet( v_nqp );

  quad -> print_info();

  const double inv_T = 1.0 / ( static_cast<double>((time_end - time_start)/time_step) + 1.0 );

  std::vector<Vector_3> velo_sol {};

  // Container for time-averaged velocity on nodes
  std::vector<Vector_3> ave_velo {};

  // Container for time-averaged velocity gradient on quadrature points
  std::vector<Tensor2_3D> ave_velo_grad {};

  // Container for time-averged Renold stress on quadrature points
  std::vector<Tensor2_3D> ave_Renold_stress {};

  std::vector<double> elem_vol {};

  // std::vector<double> asp_ratio {};

  // std::vector<int> GlobalElementID {};

  // Total value for omp
  double total_P_dir = 0, total_P_prod = 0, total_P_mod = 0, total_P_diss = 0;

  // int test_begin = 552126;
  // int test_num = 552127;

  // Calculate ave_velo and ave_velo_grad
  std::cout<<"==== Calculate ave_velo and ave_velo_grad ===="<<std::endl;
  for(int time = time_start; time <= time_end; time += time_step)
  {
    // Generate the file name
    geo_file = get_vtu_result_name(sol_bname, time, is_pvtu);

    std::cout<<"Time "<<time<<": Read "<<geo_file<<std::endl;

    if(time == time_start)
    {
      VTK_T::read_vtu_grid(geo_file, v_nFunc, v_nElem, v_ctrlPts, v_vecIEN, elem_vol);

      cout<<endl<<"Volumetric mesh contains "<<v_nElem<<" elements and "<<v_nFunc<<" vertices.\n";

      ave_velo = std::vector<Vector_3>(v_nFunc, Vector_3(0, 0, 0));

      ave_velo_grad = std::vector<Tensor2_3D>(v_nqp * v_nElem, Tensor2_3D(0, 0, 0, 0, 0, 0, 0, 0, 0));

      ave_Renold_stress = std::vector<Tensor2_3D>(v_nqp * v_nElem, Tensor2_3D(0, 0, 0, 0, 0, 0, 0, 0, 0));

      // asp_ratio = VTK_T::read_double_CellData(geo_file, "Aspect_ratio");

      // GlobalElementID = VTK_T::read_int_CellData(geo_file, "GlobalElementID");
    }

    velo_sol = VTK_T::read_Vector3_PointData(geo_file, "Velocity");

    PERIGEE_OMP_PARALLEL_FOR
    for(int nn=0; nn<v_nFunc; ++nn)
    {
      Vector_3 * node_velo = &velo_sol[nn];

      ave_velo[nn] += inv_T * *node_velo;
    }
    
    PERIGEE_OMP_PARALLEL_FOR
    for(int ee=0; ee<v_nElem; ++ee)
    {
      // // Filter
      // if(time == time_start && asp_ratio[ee] > aspect_ratio_threshold)
      // {
      //   GlobalElementID[ee] = -1;
      //   // std::cout<<"Warning: Element "<<ee<<" has a very high aspect ratio: "
      //   //          <<asp_ratio[ee]<<" > 100.0.\n";
      // }
      
      // if(GlobalElementID[ee] != -1)
      // {
        // Element container
        FEAElement * element = new FEAElement_Tet4( quad-> get_num_quadPts() );

        const int tet_node[4] { v_vecIEN[4*ee+0], v_vecIEN[4*ee+1],
                                v_vecIEN[4*ee+2], v_vecIEN[4*ee+3] };

        const double ectrl_x[4] { v_ctrlPts[3*tet_node[0] + 0], v_ctrlPts[3*tet_node[1] + 0],
                                  v_ctrlPts[3*tet_node[2] + 0], v_ctrlPts[3*tet_node[3] + 0] };

        const double ectrl_y[4] { v_ctrlPts[3*tet_node[0] + 1], v_ctrlPts[3*tet_node[1] + 1],
                                  v_ctrlPts[3*tet_node[2] + 1], v_ctrlPts[3*tet_node[3] + 1] };

        const double ectrl_z[4] { v_ctrlPts[3*tet_node[0] + 2], v_ctrlPts[3*tet_node[1] + 2],
                                  v_ctrlPts[3*tet_node[2] + 2], v_ctrlPts[3*tet_node[3] + 2] };

        element -> buildBasis(quad, ectrl_x, ectrl_y, ectrl_z);

        const double esol_u[4] { velo_sol[ tet_node[0] ].x(),
                                velo_sol[ tet_node[1] ].x(),
                                velo_sol[ tet_node[2] ].x(),
                                velo_sol[ tet_node[3] ].x() };

        const double esol_v[4] { velo_sol[ tet_node[0] ].y(),
                                velo_sol[ tet_node[1] ].y(),
                                velo_sol[ tet_node[2] ].y(),
                                velo_sol[ tet_node[3] ].y() };

        const double esol_w[4] { velo_sol[ tet_node[0] ].z(),
                                velo_sol[ tet_node[1] ].z(),
                                velo_sol[ tet_node[2] ].z(),
                                velo_sol[ tet_node[3] ].z() };

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

        delete element;
      }
    // }
  }// Loop over each time instance

  // Calculate ave_Renold_stress, P_mod, P_diss
  std::cout<<"==== Calculate ave_Renold_stress, P_dir, P_mod, P_diss ===="<<std::endl;
  for(int time = time_start; time <= time_end; time += time_step)
  {
    // Generate the file name
    geo_file = get_vtu_result_name(sol_bname, time, is_pvtu);

    std::cout<<"Time "<<time<<": Read "<<geo_file<<std::endl;

    velo_sol = VTK_T::read_Vector3_PointData(geo_file, "Velocity");

    // PERIGEE_OMP_PARALLEL
    // {
    //   // Container for P_mod, P_diss
    //   double P_mod = 0, P_diss = 0;

    #pragma omp parallel for reduction(+:total_P_mod) reduction(+:total_P_diss)
    for(int ee=0; ee<v_nElem; ++ee)
    {
      // if(GlobalElementID[ee] != -1)
      // {
        FEAElement * element = new FEAElement_Tet4( quad-> get_num_quadPts() );

        const int tet_node[4] { v_vecIEN[4*ee+0], v_vecIEN[4*ee+1],
                                v_vecIEN[4*ee+2], v_vecIEN[4*ee+3] };

        const double ectrl_x[4] { v_ctrlPts[3*tet_node[0] + 0], v_ctrlPts[3*tet_node[1] + 0],
                                  v_ctrlPts[3*tet_node[2] + 0], v_ctrlPts[3*tet_node[3] + 0] };

        const double ectrl_y[4] { v_ctrlPts[3*tet_node[0] + 1], v_ctrlPts[3*tet_node[1] + 1],
                                  v_ctrlPts[3*tet_node[2] + 1], v_ctrlPts[3*tet_node[3] + 1] };

        const double ectrl_z[4] { v_ctrlPts[3*tet_node[0] + 2], v_ctrlPts[3*tet_node[1] + 2],
                                  v_ctrlPts[3*tet_node[2] + 2], v_ctrlPts[3*tet_node[3] + 2] };

        element -> buildBasis(quad, ectrl_x, ectrl_y, ectrl_z);

        const double esol_u[4] { velo_sol[ tet_node[0] ].x(),
                                  velo_sol[ tet_node[1] ].x(),
                                  velo_sol[ tet_node[2] ].x(),
                                  velo_sol[ tet_node[3] ].x() };

        const double esol_v[4] { velo_sol[ tet_node[0] ].y(),
                                  velo_sol[ tet_node[1] ].y(),
                                  velo_sol[ tet_node[2] ].y(),
                                  velo_sol[ tet_node[3] ].y() };

        const double esol_w[4] { velo_sol[ tet_node[0] ].z(),
                                  velo_sol[ tet_node[1] ].z(),
                                  velo_sol[ tet_node[2] ].z(),
                                  velo_sol[ tet_node[3] ].z() };

        // Fluctuation on nodes
        const double fluc_u[4] { esol_u[0] - ave_velo[tet_node[0]].x(),
                                  esol_u[1] - ave_velo[tet_node[1]].x(),
                                  esol_u[2] - ave_velo[tet_node[2]].x(),
                                  esol_u[3] - ave_velo[tet_node[3]].x() };

        const double fluc_v[4] { esol_v[0] - ave_velo[tet_node[0]].y(),
                                  esol_v[1] - ave_velo[tet_node[1]].y(),
                                  esol_v[2] - ave_velo[tet_node[2]].y(),
                                  esol_v[3] - ave_velo[tet_node[3]].y() };

        const double fluc_w[4] { esol_w[0] - ave_velo[tet_node[0]].z(),
                                  esol_w[1] - ave_velo[tet_node[1]].z(),
                                  esol_w[2] - ave_velo[tet_node[2]].z(),
                                  esol_w[3] - ave_velo[tet_node[3]].z() };
        
        std::array<Vector_3, 4> node_coor { Vector_3(ectrl_x[0], ectrl_y[0], ectrl_z[0]),
                                            Vector_3(ectrl_x[1], ectrl_y[1], ectrl_z[1]),
                                            Vector_3(ectrl_x[2], ectrl_y[2], ectrl_z[2]),
                                            Vector_3(ectrl_x[3], ectrl_y[3], ectrl_z[3]) };
        TET_T::Tet4 tet_elem (node_coor);

        // node_coor[0].print();
        // node_coor[1].print();
        // node_coor[2].print();
        // node_coor[3].print();

        double e_volume = tet_elem.get_volume();
        double e_delta = 0.0;
        // if(e_volume > 1e-12)
          e_delta = std::pow(e_volume, 1.0/3.0);
        // else if(e_volume <= 1e-12 && e_volume > 0.0)
        //   e_delta = std::pow(e_volume * 1e+6, 1.0/3.0) / 100.0; // avoid zero volume
        // else
        //   SYS_T::print_fatal("Error: Element %d has zero volume.\n", ee);

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

          // cout<< "ee = "<<ee<<", qua = "<<qua<<"\n";
          // cout<< "  detJ = "<<detJ<<"\n";
          // cout<< "  volume = "<<tet_elem.get_volume()<<"\n";
          // cout<< "e_delta = "<<e_delta<<", CC = "<<CC<<"\n";
          // cout<< "ars = "<<"\n";
          // ave_Renold_stress[v_nqp * ee + qua].print();
          // cout<< "fg = "<<"\n";
          // fg.print();
          // cout<< "vg = "<<"\n";
          // vg.print();
          // cout<< "avg = "<<"\n";
          // avg.print();
          // cout<< "rate_of_strain = "<<"\n";
          // rate_of_strain.print();
          // cout<< "  S_bar = "<<S_bar<<", mu_t = "<<mu_t<<"\n";

          total_P_mod += inv_T * mu_t * detJ * qw * (vg.MatContraction() + vg.MatTContraction(vg));

          total_P_diss += inv_T * fluid_mu * detJ * qw * (fg.MatContraction() + fg.MatTContraction(fg));
        }

        delete element;
      // }
    }

    //   PERIGEE_OMP_CRITICAL
    //   {
    //     total_P_mod += P_mod;
    //     total_P_diss += P_diss;
    //   }
    // }
  }

  // Calculate P_dir, P_prod
  std::cout<<"==== Calculate P_dir, P_prod ===="<<std::endl;

  // Generate the file name
  geo_file = get_vtu_result_name(sol_bname, time_start, is_pvtu);

  // PERIGEE_OMP_PARALLEL
  // {
  //   // Container for P_dir, P_prod
  //   double P_dir = 0, P_prod = 0;
  #pragma omp parallel for reduction(+:total_P_dir) reduction(+:total_P_prod)
  for(int ee=0; ee<v_nElem; ++ee)
  {
    // if(GlobalElementID[ee] != -1)
    // {
      FEAElement * element = new FEAElement_Tet4( quad-> get_num_quadPts() );
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

        total_P_dir += fluid_mu * detJ * qw * (avg.MatContraction() + avg.MatTContraction(avg));

        total_P_prod += rho * detJ * qw * ars.MatContraction(avg);
      }

      delete element;
    // }

    // PERIGEE_OMP_CRITICAL
    // {
    //   total_P_dir += P_dir;
    //   total_P_prod += P_prod;
    // }
  }
  
  const double P_tot_2 = total_P_dir + total_P_prod + total_P_mod;
  const double P_tot_3 = total_P_dir + total_P_diss + total_P_mod;

  std::cout<<"==== Result ===="<<std::endl;
  std::cout<<"  P_dir = "<<total_P_dir<<std::endl;
  std::cout<<"  P_prod = "<<total_P_prod<<std::endl;
  std::cout<<"  P_mod = "<<total_P_mod<<std::endl;
  std::cout<<"  P_diss = "<<total_P_diss<<std::endl;
  std::cout<<"  P_tot_2 = P_dir + P_prod + P_mod = "<<P_tot_2<<std::endl;
  std::cout<<"  P_tot_3 = P_dir + P_diss + P_mod = "<<P_tot_3<<std::endl;

  delete quad;
  PetscFinalize();
  return EXIT_SUCCESS;
}
// END of MAIN Function

std::string get_vtu_result_name(const std::string &pre_fix, 
    const int &time, const int &is_pvtu)
{
  std::string name = pre_fix;

  std::ostringstream time_index;

  std::string pre_time;

  time_index.str("");
  time_index<<time;
  pre_time += time_index.str();
   
  name += pre_time;

  name += "_0_0";

  if (is_pvtu)
    name += ".pvtu";
  else
    name += ".vtu";

  return name;
}

// EOF