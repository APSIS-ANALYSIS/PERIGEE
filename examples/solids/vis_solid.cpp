#include "AGlobal_Mesh_Info.hpp"
#include "ANL_Tools.hpp"
#include "FEAElementFactory.hpp"
#include "HDF5_Reader.hpp"
#include "MaterialModel_ich_NeoHookean.hpp"
#include "MaterialModel_vol_Incompressible.hpp"
#include "MaterialModel_Mixed_Elasticity.hpp"
#include "QuadPtsFactory.hpp"
#include "VisDataPrep_Solid.hpp"
#include "VTK_Writer_Solid.hpp"

int main( int argc, char * argv[] )
{
  const std::string element_part_file = "epart.h5";
  const std::string anode_mapping_file = "node_mapping.h5";
  const std::string pnode_mapping_file = "post_node_mapping.h5";
  const std::string part_file = "./ppart/part";

  std::string disp_sol_bname("SOL_disp_");
  std::string velo_sol_bname("SOL_velo_");
  std::string pres_sol_bname("SOL_pres_");
  std::string out_bname("SOL_");

  int time_start = 0;
  int time_step = 1;
  int time_end = 1;
  bool isXML = true;
  bool isRef = false;
  bool isClean = true;

  // Read analysis code parameter if the solver_cmd.h5 exists
  hid_t prepcmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

  const int init_index = cmd_h5r -> read_intScalar("/","init_index");
  const int final_index = cmd_h5r -> read_intScalar("/","final_index");
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

  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  time_start = init_index;
  time_end = final_index;
  time_step = sol_rec_freq;

  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionReal("-dt", dt);
  SYS_T::GetOptionString("-disp_sol_bname", disp_sol_bname);
  SYS_T::GetOptionString("-velo_sol_bname", velo_sol_bname);
  SYS_T::GetOptionString("-pres_sol_bname", pres_sol_bname);
  SYS_T::GetOptionString("-out_bname", out_bname);
  SYS_T::GetOptionBool("-xml", isXML);
  SYS_T::GetOptionBool("-ref", isRef);
  SYS_T::GetOptionBool("-clean", isClean);

  // Correct time_step if it does not match with sol_rec_freq
  if( time_step % sol_rec_freq != 0 ) time_step = sol_rec_freq;

  SYS_T::cmdPrint("-disp_sol_bname:", disp_sol_bname);
  SYS_T::cmdPrint("-velo_sol_bname:", velo_sol_bname);
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

  if( isClean )
  {
    SYS_T::execute( ("rm -rf " + out_bname + "*_p*.vtu").c_str() );
    SYS_T::execute( ("rm -rf " + out_bname + "*.pvtu").c_str() );
    SYS_T::execute( ("rm -rf " + out_bname + ".pvd").c_str() );
  }

  SYS_T::print_fatal_if(size != ANL_T::get_cpu_size(part_file, rank),
      "Error: number of processors does not match with prepost! \n");

  SYS_T::commPrint("===> %d processor(s) are assigned.", size);

  auto fNode = SYS_T::make_unique<FEANode>(part_file, rank);
  auto locIEN = SYS_T::make_unique<ALocal_IEN>(part_file, rank);
  auto GMIptr = SYS_T::make_unique<AGlobal_Mesh_Info>(part_file, rank);
  auto locElem = SYS_T::make_unique<ALocal_Elem>(part_file, rank);
  auto pNode = SYS_T::make_unique<APart_Node>(part_file, rank);

  const auto elemType = GMIptr->get_elemType();
  std::unique_ptr<IQuadPts> quad = QuadPtsFactory::createVisQuadrature(elemType);
  std::unique_ptr<FEAElement> element = ElementFactory::createVolElement(elemType,
      quad->get_num_quadPts());

  quad -> print_info();

  // Material model for Cauchy stress
  const double solid_mu = 6.666666666e4;
  const double solid_rho0 = 1.0;

  std::unique_ptr<IMaterialModel_ich> imodel =
    SYS_T::make_unique<MaterialModel_ich_NeoHookean>(solid_mu);

  std::unique_ptr<IMaterialModel_vol> vmodel =
    SYS_T::make_unique<MaterialModel_vol_Incompressible>(solid_rho0);

  auto matmodel = SYS_T::make_unique<MaterialModel_Mixed_Elasticity>(
      std::move(vmodel), std::move(imodel));

  std::unique_ptr<IVisDataPrep> visprep = SYS_T::make_unique<VisDataPrep_Solid>();
  visprep->print_info();

  double ** solArrays = new double * [visprep->get_ptarray_size()];
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    solArrays[ii] = new double [pNode->get_nlocghonode() * visprep->get_ptarray_comp_length(ii)];

  auto vtk_w = SYS_T::make_unique<VTK_Writer_Solid>( GMIptr->get_nElem(),
      GMIptr->get_nLocBas(), element_part_file, std::move(matmodel) );

  std::ostringstream time_index;

  for(int time = time_start; time<=time_end; time += time_step)
  {
    std::string disp_name_to_read(disp_sol_bname);
    std::string velo_name_to_read(velo_sol_bname);
    std::string pres_name_to_read(pres_sol_bname);
    std::string name_to_write(out_bname);
    time_index.str("");
    time_index << 900000000 + time;
    disp_name_to_read.append(time_index.str());
    velo_name_to_read.append(time_index.str());
    pres_name_to_read.append(time_index.str());
    name_to_write.append(time_index.str());

    SYS_T::commPrint("Time %d: Read %s %s %s and Write %s \n",
        time, disp_name_to_read.c_str(), pres_name_to_read.c_str(),
        velo_name_to_read.c_str(), name_to_write.c_str() );

    std::vector<std::string> sol_names;
    sol_names.push_back(disp_name_to_read);
    sol_names.push_back(pres_name_to_read);
    sol_names.push_back(velo_name_to_read);

    visprep->get_pointArray(sol_names, anode_mapping_file, pnode_mapping_file,
        pNode.get(), GMIptr->get_nFunc(), solArrays);

    vtk_w->writeOutput( fNode.get(), locIEN.get(), locElem.get(), visprep.get(),
        element.get(), quad.get(), solArrays, rank, size,
        time * dt, out_bname, name_to_write, isXML, isRef );
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    delete [] solArrays[ii];
  delete [] solArrays;

  PetscFinalize();
  return EXIT_SUCCESS;
}
