// ============================================================================
// vis_fsi_smooth.cpp
//
// This routine provides visualization of the solid sub-domain in FSI problems.
//
// Date: Mar 22 2024
// ============================================================================

#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "APart_Node_FSI.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "QuadPts_vis_hex8.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Hex8.hpp"
#include "VisDataPrep_Stress_Recovery.hpp"  
#include "VTK_Writer_Stress_Recovery.hpp"
#include "MaterialModel_GOH06_ST91_Mixed.hpp"
#include "MaterialModel_GOH06_Incompressible_Mixed.hpp"
#include "MaterialModel_GOH14_ST91_Mixed.hpp"
#include "Tissue_property.hpp"

int main ( int argc , char * argv[] )
{
  const std::string element_part_file = "epart.h5";
  const std::string an_v_mapping_file = "node_mapping_v.h5";
  const std::string an_p_mapping_file = "node_mapping_p.h5";
  const std::string pn_v_mapping_file = "post_node_mapping_v.h5";
  const std::string pn_p_mapping_file = "post_node_mapping_p.h5";
  const std::string part_v_file="./ppart/postpart_v";
  const std::string part_p_file="./ppart/postpart_p";

  std::string grad_sol_bname("SOL_Grad_disp_");
  
  std::string disp_sol_bname("SOL_disp_");
  std::string pres_sol_bname("SOL_pres_");
  std::string out_bname("VIS_Stress_");

  int time_start = 0;
  int time_step = 1;
  int time_end = 1;
  bool isXML = true;
  bool isRef = false;
  bool isClean = true;

  // Load analysis code parameter from solver_cmd.h5 file
  hid_t prepcmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

  double dt = cmd_h5r -> read_doubleScalar("/","init_step");

  const int sol_rec_freq = cmd_h5r -> read_intScalar("/", "sol_record_freq");

  delete cmd_h5r; H5Fclose(prepcmd_file);

  // ===== PETSc Initialization =====
#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  // Material parameters
  bool is_read_material = true;    // bool flag to decide if one wants to read material model from h5 file
  double solid_density = 1.0e-3;
  double solid_E = 121716;
  double solid_nu = 0.47;
  double solid_mu = 4.14e4;
  double solid_f1the = 47.0;
  double solid_f1phi = 90.0;
  double solid_f2the = -47.0;
  double solid_f2phi = 90.0;
  double solid_fk1 = 2.0188e6;
  double solid_fk2 = 20.0;
  double solid_fkd = 0.226;

  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionReal("-dt", dt);
  SYS_T::GetOptionString("-grad_sol_bname", grad_sol_bname);
  SYS_T::GetOptionString("-disp_sol_bname", disp_sol_bname);
  SYS_T::GetOptionString("-pres_sol_bname", pres_sol_bname);
  SYS_T::GetOptionString("-out_bname", out_bname);
  SYS_T::GetOptionBool("-xml", isXML);
  SYS_T::GetOptionBool("-ref", isRef);
  SYS_T::GetOptionBool("-clean", isClean);
  SYS_T::GetOptionBool("-is_read_material", is_read_material);
  SYS_T::GetOptionReal( "-sl_density", solid_density);
  SYS_T::GetOptionReal( "-sl_E", solid_E);
  SYS_T::GetOptionReal( "-sl_nu", solid_nu);
  SYS_T::GetOptionReal( "-sl_mu", solid_mu);
  SYS_T::GetOptionReal( "-sl_f1the", solid_f1the);
  SYS_T::GetOptionReal( "-sl_f1phi", solid_f1phi);
  SYS_T::GetOptionReal( "-sl_f2the", solid_f2the);
  SYS_T::GetOptionReal( "-sl_f2phi", solid_f2phi);
  SYS_T::GetOptionReal( "-sl_fk1", solid_fk1);
  SYS_T::GetOptionReal( "-sl_fk2", solid_fk2);
  SYS_T::GetOptionReal( "-sl_fkd", solid_fkd);
 

  // Correct time_step if it does not match with sol_rec_freq
  if( time_step % sol_rec_freq != 0 ) time_step = sol_rec_freq;

  SYS_T::cmdPrint("-grad_sol_bname:", grad_sol_bname);
  SYS_T::cmdPrint("-disp_sol_bname",  disp_sol_bname);
  SYS_T::cmdPrint("-pres_sol_bname:", pres_sol_bname);
  SYS_T::cmdPrint("-out_bname:", out_bname);
  SYS_T::cmdPrint("-time_start:", time_start);
  SYS_T::cmdPrint("-time_step:", time_step);
  SYS_T::cmdPrint("-time_end:", time_end);
  SYS_T::cmdPrint("-dt:",dt);
  if(isXML) SYS_T::commPrint("-xml: true \n");
  else SYS_T::commPrint("-xml: false \n");

  if(isRef) SYS_T::commPrint("-ref: true \n");
  else SYS_T::commPrint("-ref: false \n");

  if(isClean) SYS_T::commPrint("-clean: true \n");
  else SYS_T::commPrint("-clean: false \n");

  // If demand cleaning, remove all previous visualization files
  if( isClean )
  {
    SYS_T::execute("rm -rf VIS_S_*_p*.vtu");
    SYS_T::execute("rm -rf VIS_S_*.pvtu");
    SYS_T::execute("rm -rf VIS_S_.pvd");
  }

  if( is_read_material )
  {
    SYS_T::commPrint("-is_read_material: true \n");
    SYS_T::file_check("material_model.h5" );
    SYS_T::commPrint("material_model.h5 found. \n");
  }
  else
  {
    SYS_T::commPrint("-is_read_material: false \n");
    SYS_T::cmdPrint("-sl_density:", solid_density);
    SYS_T::cmdPrint("-sl_E:", solid_E);
    SYS_T::cmdPrint("-sl_nu:", solid_nu);
    SYS_T::cmdPrint("-sl_mu:", solid_mu);
    SYS_T::cmdPrint("-sl_f1the:", solid_f1the);
    SYS_T::cmdPrint("-sl_f1phi:", solid_f1phi);
    SYS_T::cmdPrint("-sl_f2the:", solid_f2the);
    SYS_T::cmdPrint("-sl_f2phi:", solid_f2phi);
    SYS_T::cmdPrint("-sl_fk1:", solid_fk1);
    SYS_T::cmdPrint("-sl_fk2:", solid_fk2);
    SYS_T::cmdPrint("-sl_fkd:", solid_fkd);
  }

  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_v_file, 0);

  SYS_T::print_fatal_if(size != PartBasic->get_cpu_size(), "Error: number of processors does not match with prepost! \n");

  SYS_T::commPrint("===> %d processor(s) are assigned.", size);

  FEANode * fNode = new FEANode(part_v_file, rank);

  ALocal_IEN * locIEN_v = new ALocal_IEN(part_v_file, rank);
  ALocal_IEN * locIEN_p = new ALocal_IEN(part_p_file, rank);

  IAGlobal_Mesh_Info * GMIptr_v = new AGlobal_Mesh_Info_FEM_3D(part_v_file, rank);
  IAGlobal_Mesh_Info * GMIptr_p = new AGlobal_Mesh_Info_FEM_3D(part_p_file, rank);

  ALocal_Elem * locElem = new ALocal_Elem(part_v_file, rank);

  APart_Node * pNode_v = new APart_Node(part_v_file, rank);
  APart_Node * pNode_p = new APart_Node(part_p_file, rank);

  Tissue_property * tp_data = new Tissue_property(part_v_file, rank);

  // Allocate the quadrature rule and element container
  IQuadPts * quad = nullptr;
  FEAElement * element = nullptr; 

  // We assume that the same element type is used for pressure and velocity
  if( GMIptr_v->get_elemType() == 501 )
  {
    quad = new QuadPts_vis_tet4();
    element = new FEAElement_Tet4( quad-> get_num_quadPts() );
  }
  else if( GMIptr_v->get_elemType() == 601 )
  {
    quad = new QuadPts_vis_hex8();
    element = new FEAElement_Hex8( quad-> get_num_quadPts() );
  }
  else SYS_T::print_fatal( "Error: unsupported element type \n" );

  quad -> print_info();

  // material model
  IMaterialModel * matmodel;
  if( is_read_material )
  {
    // load solid_nu from h5 file
    hid_t model_file = H5Fopen("material_model.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

    HDF5_Reader * model_h5r = new HDF5_Reader( model_file );

    solid_nu = model_h5r -> read_doubleScalar("/", "nu");

    delete model_h5r; H5Fclose(model_file);

    if( solid_nu == 0.5 )
    {
      matmodel = new MaterialModel_GOH06_Incompressible_Mixed( "material_model.h5" );
    }
    else
    {
      //matmodel = new MaterialModel_GOH06_ST91_Mixed( "material_model.h5" );
      matmodel = new MaterialModel_GOH14_ST91_Mixed( "material_model.h5" );
    }
  }
  else
  {
    if( solid_nu == 0.5 )
    {
      matmodel = new MaterialModel_GOH06_Incompressible_Mixed( solid_density, solid_mu,
        solid_f1the, solid_f1phi, solid_f2the, solid_f2phi, solid_fk1, solid_fk2, solid_fkd );
    }
    else
    {
      //matmodel = new MaterialModel_GOH06_ST91_Mixed( solid_density, solid_E, solid_nu,
      //  solid_f1the, solid_f1phi, solid_f2the, solid_f2phi, solid_fk1, solid_fk2, solid_fkd );
      matmodel = new MaterialModel_GOH14_ST91_Mixed( solid_density, solid_E, solid_nu,
        solid_f1the, solid_f1phi, solid_f2the, solid_f2phi, solid_fk1, solid_fk2, solid_fkd );
    }
  }

  matmodel -> print_info();

  // For the solid subdomain, we need to prepare a mapping from the
  // FSI nodal index to the solid subdomain nodal index
  std::vector<int> subdomain_nodes; subdomain_nodes.clear();
  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    if( locElem -> get_elem_tag(ee) == 1 )
      VEC_T::insert_end( subdomain_nodes, locIEN_v -> get_LIEN(ee) );
  }
  VEC_T::sort_unique_resize( subdomain_nodes );
  const int num_subdomain_nodes = static_cast<int>( subdomain_nodes.size() );

  const int nelem = locElem -> get_nlocalele();
  const int nlocbas = locIEN_v -> get_stride();
  std::vector<int> sIEN( nelem * nlocbas );
  for( int ee=0; ee<nelem; ++ee )
    for(int ii=0; ii<nlocbas; ++ii)
      sIEN[ee*nlocbas+ii] = VEC_T::get_pos( subdomain_nodes, locIEN_v->get_LIEN(ee, ii) );

  // Visualization prepration
  IVisDataPrep * visprep = new VisDataPrep_Smooth_FSI();
  visprep->print_info();
  double ** pointArrays = new double * [3];
  pointArrays[0] = new double [pNode_v->get_nlocghonode() * 9]; // stress, grad_disp or F
  pointArrays[1] = new double [pNode_p->get_nlocghonode() * 1]; // vonMise or pres
  pointArrays[2] = new double [pNode_p->get_nlocghonode() * 3]; // disp

  VTK_Writer_Smooth_FSI * vtk_w = new VTK_Writer_Smooth_FSI( GMIptr_v->get_nElem(),
      element->get_nLocBas(), element_part_file );
  std::ostringstream time_index;
  for(int time = time_start; time<=time_end; time += time_step)
  {
    std::string gradsol_name_to_read(grad_sol_bname);
    std::string disp_name_to_read(disp_sol_bname);  
    std::string pres_name_to_read(pres_sol_bname);
    std::string name_to_write(out_bname);
    time_index.str("");
    time_index<< 900000000 + time;
    gradsol_name_to_read.append(time_index.str());
    disp_name_to_read.append(time_index.str());
    pres_name_to_read.append(time_index.str());
    name_to_write.append(time_index.str());

    SYS_T::commPrint("Time %d: Read %s %s and Write %s \n",
        time, gradsol_name_to_read.c_str(), disp_name_to_read.c_str(), pres_name_to_read.c_str(),
        name_to_write.c_str() );
    visprep->get_pointArray(gradsol_name_to_read, disp_name_to_read, pres_name_to_read,
        an_v_mapping_file, an_p_mapping_file,
        pn_v_mapping_file, pn_p_mapping_file,
        pNode_v, pNode_p, GMIptr_v->get_nFunc(), GMIptr_p->get_nFunc(),
        pointArrays);
    vtk_w->writeOutput( fNode, locIEN_v, locIEN_p, sIEN, locElem,
        visprep, matmodel, element, quad, tp_data, pointArrays, rank, size,
        num_subdomain_nodes, time * dt, out_bname, name_to_write, isXML );
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // Clean up memory
  delete quad; delete element; delete visprep; delete matmodel; delete tp_data;
  delete fNode; delete locIEN_v; delete locIEN_p; delete GMIptr_v; delete GMIptr_p;
  delete PartBasic; delete locElem; delete pNode_v; delete pNode_p;
  delete [] pointArrays[0]; delete [] pointArrays[1];
  delete [] pointArrays; delete vtk_w;
  PetscFinalize();

  return EXIT_SUCCESS ;
}

//EOF
