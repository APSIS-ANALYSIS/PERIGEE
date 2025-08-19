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
#include "HDF5_Writer.hpp"
#include "HDF5_Reader.hpp"

std::string get_vtu_result_name(const std::string &pre_fix, const int &time,
  const int &is_pvtu = 0);

template<typename T> std::vector<T> get_inherited(const std::vector<T> original, const std::vector<int> inherit)
{
  int length = VEC_T::get_size(inherit);
  std::vector<T> output (length);
  PERIGEE_OMP_PARALLEL_FOR
  for(int ii=0; ii<length; ++ii)
    output[ii] = original[inherit[ii]];

  return output;
}

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

  int read_h5f = 0;
  std::string h5f_name = "IEN_inherit.h5";

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
  SYS_T::GetOptionInt("-read_h5f", read_h5f);
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

  std::vector<Vector_3> original_velo_sol {};
  std::vector<Vector_3> velo_sol {};

  std::vector<double> original_du_dx {};
  std::vector<double> du_dx {};

  std::vector<double> original_du_dy {};
  std::vector<double> du_dy {};

  std::vector<double> original_du_dz {};
  std::vector<double> du_dz {};

  std::vector<double> original_dv_dx {};
  std::vector<double> dv_dx {};

  std::vector<double> original_dv_dy {};
  std::vector<double> dv_dy {};

  std::vector<double> original_dv_dz {};
  std::vector<double> dv_dz {};

  std::vector<double> original_dw_dx {};
  std::vector<double> dw_dx {};

  std::vector<double> original_dw_dy {};
  std::vector<double> dw_dy {};

  std::vector<double> original_dw_dz {};
  std::vector<double> dw_dz {};

  // Container for time-averaged velocity on nodes
  std::vector<Vector_3> ave_velo {};

  // Container for time-averaged velocity gradient on quadrature points
  std::vector<Tensor2_3D> ave_velo_grad {};

  // Container for time-averged Renold stress on quadrature points
  std::vector<Tensor2_3D> ave_Renold_stress {};

  std::vector<double> elem_vol {};

  std::vector<int> inherit {};

  std::vector<int> is_considered {};

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
      if(read_h5f)
      {
        VTK_T::read_vtu_grid(geo_file, v_nFunc, v_ctrlPts);

        hid_t file_id = H5Fopen(h5f_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

        auto h5r = SYS_T::make_unique<HDF5_Reader>( file_id );

        v_vecIEN = h5r -> read_intVector("/", "ien");

        inherit = h5r -> read_intVector("/", "inherit");

        v_nElem = VEC_T::get_size(inherit);

        H5Fclose( file_id );
      }
      else
      {
        VTK_T::read_vtu_grid(geo_file, v_nFunc, v_nElem, v_ctrlPts, v_vecIEN, elem_vol, inherit);

        SYS_T::execute("rm -rf IEN_inherit.h5");

        hid_t file_id = H5Fcreate(h5f_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT, H5P_DEFAULT);

        auto h5w = SYS_T::make_unique<HDF5_Writer>( file_id );

        h5w -> write_intVector("ien", v_vecIEN);

        h5w -> write_intVector("inherit", inherit);

        H5Fclose( file_id );
      }

      cout<<endl<<"Volumetric mesh contains "<<v_nElem<<" elements and "<<v_nFunc<<" vertices.\n";

      ave_velo = std::vector<Vector_3>(v_nElem, Vector_3(0, 0, 0));

      ave_velo_grad = std::vector<Tensor2_3D>(v_nElem, Tensor2_3D(0, 0, 0, 0, 0, 0, 0, 0, 0));

      ave_Renold_stress = std::vector<Tensor2_3D>(v_nElem, Tensor2_3D(0, 0, 0, 0, 0, 0, 0, 0, 0));

      // asp_ratio = VTK_T::read_double_CellData(geo_file, "Aspect_ratio");

      // GlobalElementID = VTK_T::read_int_CellData(geo_file, "GlobalElementID");

      is_considered = std::vector<int>(v_nElem, 0);

      PERIGEE_OMP_PARALLEL_FOR
      for(int ee=0; ee<v_nElem; ++ee)
      {
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
        
        std::array<Vector_3, 4> node_coor { Vector_3(ectrl_x[0], ectrl_y[0], ectrl_z[0]),
                                            Vector_3(ectrl_x[1], ectrl_y[1], ectrl_z[1]),
                                            Vector_3(ectrl_x[2], ectrl_y[2], ectrl_z[2]),
                                            Vector_3(ectrl_x[3], ectrl_y[3], ectrl_z[3]) };
        TET_T::Tet4 tet_elem (node_coor);

        Vector_3 centroid = tet_elem.get_centroid();

        if(centroid.y()<=0.03)
          is_considered[ee] = 1;
      }
    }

    original_velo_sol = VTK_T::read_Vector3_CellData(geo_file, "Velocity");
    velo_sol = get_inherited(original_velo_sol, inherit);

    original_du_dx = VTK_T::read_double_CellData(geo_file, "dX-Velocity.dx");
    du_dx = get_inherited(original_du_dx, inherit);

    original_du_dy = VTK_T::read_double_CellData(geo_file, "dX-Velocity.dy");
    du_dy = get_inherited(original_du_dy, inherit);

    original_du_dz = VTK_T::read_double_CellData(geo_file, "dX-Velocity.dz");
    du_dz = get_inherited(original_du_dz, inherit);

    original_dv_dx = VTK_T::read_double_CellData(geo_file, "dY-Velocity.dx");
    dv_dx = get_inherited(original_dv_dx, inherit);

    original_dv_dy = VTK_T::read_double_CellData(geo_file, "dY-Velocity.dy");
    dv_dy = get_inherited(original_dv_dy, inherit);

    original_dv_dz = VTK_T::read_double_CellData(geo_file, "dY-Velocity.dz");
    dv_dz = get_inherited(original_dv_dz, inherit);

    original_dw_dx = VTK_T::read_double_CellData(geo_file, "dZ-Velocity.dx");
    dw_dx = get_inherited(original_dw_dx, inherit);

    original_dw_dy = VTK_T::read_double_CellData(geo_file, "dZ-Velocity.dy");
    dw_dy = get_inherited(original_dw_dy, inherit);

    original_dw_dz = VTK_T::read_double_CellData(geo_file, "dZ-Velocity.dz");
    dw_dz = get_inherited(original_dw_dz, inherit);

    PERIGEE_OMP_PARALLEL_FOR
    for(int ee=0; ee<v_nElem; ++ee)
    {
      if(is_considered[ee])
      {
        Vector_3 * node_velo = &velo_sol[ee];

        ave_velo[ee] += inv_T * *node_velo;

        Tensor2_3D velo_grad (du_dx[ee], du_dy[ee], du_dz[ee],
                              dv_dx[ee], dv_dy[ee], dv_dz[ee],
                              dw_dx[ee], dw_dy[ee], dw_dz[ee]);

        ave_velo_grad[ee] += inv_T * velo_grad;
      }
    }
  }

  // Calculate ave_Renold_stress, P_mod, P_diss
  std::cout<<"==== Calculate ave_Renold_stress, P_dir, P_mod, P_diss ===="<<std::endl;
  for(int time = time_start; time <= time_end; time += time_step)
  {
    // Generate the file name
    geo_file = get_vtu_result_name(sol_bname, time, is_pvtu);

    std::cout<<"Time "<<time<<": Read "<<geo_file<<std::endl;

    original_velo_sol = VTK_T::read_Vector3_CellData(geo_file, "Velocity");
    velo_sol = get_inherited(original_velo_sol, inherit);

    original_du_dx = VTK_T::read_double_CellData(geo_file, "dX-Velocity.dx");
    du_dx = get_inherited(original_du_dx, inherit);

    original_du_dy = VTK_T::read_double_CellData(geo_file, "dX-Velocity.dy");
    du_dy = get_inherited(original_du_dy, inherit);

    original_du_dz = VTK_T::read_double_CellData(geo_file, "dX-Velocity.dz");
    du_dz = get_inherited(original_du_dz, inherit);

    original_dv_dx = VTK_T::read_double_CellData(geo_file, "dY-Velocity.dx");
    dv_dx = get_inherited(original_dv_dx, inherit);

    original_dv_dy = VTK_T::read_double_CellData(geo_file, "dY-Velocity.dy");
    dv_dy = get_inherited(original_dv_dy, inherit);

    original_dv_dz = VTK_T::read_double_CellData(geo_file, "dY-Velocity.dz");
    dv_dz = get_inherited(original_dv_dz, inherit);

    original_dw_dx = VTK_T::read_double_CellData(geo_file, "dZ-Velocity.dx");
    dw_dx = get_inherited(original_dw_dx, inherit);

    original_dw_dy = VTK_T::read_double_CellData(geo_file, "dZ-Velocity.dy");
    dw_dy = get_inherited(original_dw_dy, inherit);

    original_dw_dz = VTK_T::read_double_CellData(geo_file, "dZ-Velocity.dz");
    dw_dz = get_inherited(original_dw_dz, inherit);

    // PERIGEE_OMP_PARALLEL
    // {
    //   // Container for P_mod, P_diss
    //   double P_mod = 0, P_diss = 0;

    #pragma omp parallel for reduction(+:total_P_mod) reduction(+:total_P_diss)
    for(int ee=0; ee<v_nElem; ++ee)
    {
      if(is_considered[ee])
      {
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
        
        std::array<Vector_3, 4> node_coor { Vector_3(ectrl_x[0], ectrl_y[0], ectrl_z[0]),
                                            Vector_3(ectrl_x[1], ectrl_y[1], ectrl_z[1]),
                                            Vector_3(ectrl_x[2], ectrl_y[2], ectrl_z[2]),
                                            Vector_3(ectrl_x[3], ectrl_y[3], ectrl_z[3]) };
        TET_T::Tet4 tet_elem (node_coor);

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

        Vector_3 fluc = velo_sol[ee] - ave_velo[ee];
        double fu = fluc.x();
        double fv = fluc.y();
        double fw = fluc.z();

        ave_Renold_stress[ee] -= inv_T * Tensor2_3D{ fu*fu, fu*fv, fu*fw,
                                                      fv*fu, fv*fv, fv*fw,
                                                      fw*fu, fw*fv, fw*fw };

        
        Tensor2_3D vg (du_dx[ee], du_dy[ee], du_dz[ee],
                      dv_dx[ee], dv_dy[ee], dv_dz[ee],
                      dw_dx[ee], dw_dy[ee], dw_dz[ee]);

        Tensor2_3D avg = ave_velo_grad[ee];
        Tensor2_3D avg_T = avg;
        avg_T.transpose();

        Tensor2_3D rate_of_strain = 0.5 * (avg + avg_T);
        double S_bar = std::sqrt(2 * (rate_of_strain.MatContraction()));

        double mu_t = CC * S_bar;
        
        Tensor2_3D fg = vg - avg;

        for(int qua=0; qua<v_nqp; ++qua)
        {
          double detJ = element -> get_detJac(qua);
          double qw = quad -> get_qw(qua);

          total_P_mod += inv_T * mu_t * detJ * qw * (vg.MatContraction() + vg.MatTContraction(vg));

          total_P_diss += inv_T * fluid_mu * detJ * qw * (fg.MatContraction() + fg.MatTContraction(fg));
        }

        delete element;
      }
    }
  }

  // Calculate P_dir, P_prod
  std::cout<<"==== Calculate P_dir, P_prod ===="<<std::endl;

  // Generate the file name
  // geo_file = get_vtu_result_name(sol_bname, time_start, is_pvtu);

  #pragma omp parallel for reduction(+:total_P_dir) reduction(+:total_P_prod)
  for(int ee=0; ee<v_nElem; ++ee)
  {
    if(is_considered[ee])
    {
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

      Tensor2_3D avg = ave_velo_grad[ee];
      Tensor2_3D ars = ave_Renold_stress[ee];

      for(int qua=0; qua<v_nqp; ++qua)
      {
        double detJ = element -> get_detJac(qua);
        double qw = quad -> get_qw(qua);

        total_P_dir += fluid_mu * detJ * qw * (avg.MatContraction() + avg.MatTContraction(avg));

        total_P_prod += rho * detJ * qw * ars.MatContraction(avg);
      }

      delete element;
    }
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