// ============================================================================
// vis_ns.cpp
//
// Visualization driver for the NS equation solver.
// 
// Author: Ju Liu, liujuy@gmail.com
// Date Created: Feb. 12 2020
// ============================================================================
#include "AGlobal_Mesh_Info.hpp"
#include "ANL_Tools.hpp"
#include "APart_Node.hpp"
#include "QuadPtsFactory.hpp"
#include "FEAElementFactory.hpp"
#include "VisDataPrep_NS.hpp"
#include "VTK_Writer_NS.hpp"

int main( int argc, char * argv[] )
{
  const std::string element_part_file = "epart.h5";
  const std::string part_file="postpart";
  
  std::string sol_bname("SOL_");
  std::string out_bname = sol_bname;
  int time_start = 0, time_step = 1, time_end = 1;
  bool isXML = true, isRestart = false;

  // Read analysis code parameter if the solver_cmd.h5 exists
  double dt = HDF5_T::read_doubleScalar("solver_cmd.h5", "/", "init_step"); 

  // ===== Initialize the MPI run =====
#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  PetscMPIInt rank = SYS_T::get_MPI_rank();
  PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionReal("-dt", dt);
  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionString("-out_bname", out_bname);
  SYS_T::GetOptionBool("-xml", isXML);
  SYS_T::GetOptionBool("-restart", isRestart);
  
  SYS_T::commPrint("=== Command line arguments ===\n");
  SYS_T::cmdPrint("-sol_bname:", sol_bname);
  SYS_T::cmdPrint("-out_bname:", out_bname);
  SYS_T::cmdPrint("-time_start:", time_start);
  SYS_T::cmdPrint("-time_step:", time_step);
  SYS_T::cmdPrint("-time_end:", time_end);
  SYS_T::cmdPrint("-dt:",dt);
  if(isXML) SYS_T::commPrint("-xml: true \n");
  else SYS_T::commPrint("-xml: false \n");

  if(isRestart) SYS_T::commPrint("-restart: true \n");
  else SYS_T::commPrint("-restart: false \n");
  SYS_T::commPrint("==============================\n");
  
  // Clean the visualization files if not restart
  if( !isRestart )
  {
    SYS_T::execute("rm -rf *_p*.vtu");
    SYS_T::execute("rm -rf *.pvtu");
    SYS_T::execute("rm -rf *_.pvd");
  }
  
  auto fNode = SYS_T::make_unique<FEANode>(part_file, rank);
  
  auto locIEN = SYS_T::make_unique<ALocal_IEN>(part_file, rank);
  
  auto GMIptr = SYS_T::make_unique<AGlobal_Mesh_Info>(part_file,rank);
  
  auto locElem = SYS_T::make_unique<ALocal_Elem>(part_file, rank);
  
  auto pNode = SYS_T::make_unique<APart_Node>(part_file, rank);
  
  SYS_T::print_fatal_if(size != ANL_T::get_cpu_size(part_file, rank), "Error: number of processors does not match with prepost! \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for:", size);

  // Allocate the quadrature rule and element container
  const auto elemType = GMIptr->get_elemType();
  std::unique_ptr<IQuadPts> quad = QuadPtsFactory::createVisQuadrature(elemType);
  std::unique_ptr<FEAElement> element = ElementFactory::createVolElement(elemType,
    quad->get_num_quadPts());

  // Print the sampling points on screen
  quad -> print_info();

  // Create the visualization data object
  std::unique_ptr<IVisDataPrep> visprep = SYS_T::make_unique<VisDataPrep_NS>();

  visprep->print_info();
 
  // Allocate the container to store the solution values into physical fields.
  double ** solArrays = new double * [visprep->get_ptarray_size()];
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    solArrays[ii] = new double [pNode->get_nlocghonode() * visprep->get_ptarray_comp_length(ii)];

  const auto epart_map = VIS_T::read_epart( element_part_file, GMIptr->get_nElem() );

  const auto anode_mapping = HDF5_T::read_intVector("node_mapping.h5", "/", "old_2_new");
  const auto pnode_mapping = HDF5_T::read_intVector("post_node_mapping.h5", "/", "new_2_old");

  for(int time = time_start; time<=time_end; time+= time_step)
  {
    const std::string suffix = std::to_string(900000000 + time);
    const std::string name_to_read  = sol_bname + suffix;
    const std::string name_to_write = out_bname + suffix;

    SYS_T::commPrint("Time %d: Read %s and Write %s \n",
        time, name_to_read.c_str(), name_to_write.c_str() );

    visprep->get_pointArray(name_to_read, anode_mapping, pnode_mapping,
        pNode.get(), solArrays);

    VTK_Writer_NS::writeOutput( fNode.get(), locIEN.get(), locElem.get(),
        visprep.get(), element.get(), quad.get(), solArrays, epart_map,
        GMIptr->get_nLocBas(),
        rank, size, time * dt, sol_bname, out_bname, name_to_write, isXML );
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== Clean the memory =====
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    delete [] solArrays[ii];
  delete [] solArrays;
  
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
