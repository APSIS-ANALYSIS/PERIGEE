#include "Global_Part_METIS.hpp"

Global_Part_METIS::Global_Part_METIS( const int &cpu_size,
        const int &in_ncommon, const bool &isDualGraph,
        const IMesh * const &mesh,
        const IIEN * const &IEN,
        const std::string &element_part_name,
        const std::string &node_part_name )
: isMETIS(true), isDual(isDualGraph), dual_edge_ncommon(in_ncommon)
{
  const idx_t nElem = mesh->get_nElem();
  const idx_t nFunc = mesh->get_nFunc();
  const idx_t nLocBas = mesh->get_nLocBas();
  
  if(cpu_size <= 1)
  {
    std::cerr<<"ERROR: METIS cannot handle partition graph into one subdomain. \n";
    exit(1);
  }

  if(isDualGraph)
  {
    if( in_ncommon < 1 || in_ncommon > nLocBas - 1 )
    {
      std::cerr<<"ERROR: ncommon: "<<in_ncommon<<" is wrong! \n";
      exit(1);
    }
  }

  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING] = 0;

  std::cout<<"-- allocating eptr and eind arrays... \n";

  const idx_t eptr_size = nElem + 1;
  const idx_t eind_size = nElem * nLocBas; 

  idx_t * eptr = new idx_t [eptr_size];
  
  std::cout<<"---- eptr allocated with length "<<eptr_size<<" and size: ";
  SYS_T::print_mem_size(double(eptr_size) * sizeof(idx_t));
  std::cout<<"\n";

  idx_t * eind = new idx_t [eind_size];

  std::cout<<"---- eind allocated with length "<<eind_size<<" and size: ";
  SYS_T::print_mem_size(double(eind_size) * sizeof(idx_t));
  std::cout<<"\n";

  clock_t time_tracker = clock();
  
  for( idx_t ee=0; ee<nElem; ++ee )
  {
    eptr[ee] = ee * nLocBas;
    for(int ii=0; ii<nLocBas; ++ii)
      eind[ee*nLocBas + ii] = IEN->get_IEN(ee, ii);
  }
  eptr[nElem] = nElem * nLocBas;
  
  time_tracker = clock() - time_tracker;

  std::cout<<"---- eptr eind generated, taking ";
  std::cout<<((double) time_tracker)/CLOCKS_PER_SEC<<" seconds. \n";

  epart = new idx_t [nElem];
  npart = new idx_t [nFunc];
  std::cout<<"---- epart and npart vector has been allocated. \n";

  idx_t ne = nElem;
  idx_t nn = nFunc;
  idx_t nparts = cpu_size;
  idx_t ncommon = in_ncommon;

  int metis_result;
  idx_t objval;
  
  time_tracker = clock();
  if( isDualGraph )
  {
    std::cout<<"---- calling METIS_PartMeshDual ... \n";
    metis_result = METIS_PartMeshDual( &ne, &nn, eptr, eind, NULL,
        NULL, &ncommon, &nparts, NULL, options, &objval, epart, npart );
  }
  else
  {
    std::cout<<"---- calling METIS_PartMeshNodal ... \n";
    metis_result = METIS_PartMeshNodal( &ne, &nn, eptr, eind, NULL,
        NULL, &nparts, NULL, NULL, &objval, epart, npart );
  }

  if(metis_result != METIS_OK)
  {
    std::cerr<<"ERROR: PARTITION FAILED: "<<std::endl;
    switch( metis_result)
    {
      case METIS_ERROR_INPUT:
        std::cout << "METIS_ERROR_INPUT" << std::endl;
        break;
      case METIS_ERROR_MEMORY:
        std::cout << "METIS_ERROR_MEMORY" << std::endl;
        break;
      case METIS_ERROR:
        std::cout << "METIS_ERROR" << std::endl;
        break;
      default:
        break;
    }
    exit(1);
  }

  delete [] eptr; eptr = nullptr;
  delete [] eind; eind = nullptr;
  
  time_tracker = clock() - time_tracker;

  std::cout<<"-- METIS partition successfully completed, taking ";
  std::cout<<((double) time_tracker)/CLOCKS_PER_SEC<<" seconds. \n";
 
  time_tracker = clock();
  std::cout<<"-- writing epart file takes "; 
  write_part_hdf5(element_part_name, epart, nElem, cpu_size, isDualGraph, in_ncommon, true );
  time_tracker = clock() - time_tracker;
  std::cout<<((double) time_tracker)/CLOCKS_PER_SEC<<" seconds. \n";
  
  time_tracker = clock();
  std::cout<<"-- writing npart file takes "; 
  write_part_hdf5(node_part_name, npart, nFunc, cpu_size, isDualGraph, in_ncommon, true );
  time_tracker = clock() - time_tracker;
  std::cout<<((double) time_tracker)/CLOCKS_PER_SEC<<" seconds. \n";

  std::cout<<"=== Global partition generated. \n";
}

Global_Part_METIS::~Global_Part_METIS()
{
  delete [] epart; delete [] npart; epart = nullptr; npart = nullptr;
}

void Global_Part_METIS::write_part_hdf5( const std::string &fileName,
    const idx_t * const &part_in,
    const int &part_size, const int &cpu_size,
    const bool &part_isdual, const int &in_ncommon,
    const bool &isMETIS ) const
{
  std::string fName( fileName ); fName.append(".h5");

  // file creation
  hid_t file_id = H5Fcreate( fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  h5w -> write_intScalar("part_size", part_size);
  h5w -> write_intScalar("cpu_size", cpu_size);

  h5w->write_intScalar("part_isdual", ( part_isdual ? 1 : 0 ) );
  h5w->write_intScalar("in_ncommon", in_ncommon);

  h5w->write_intScalar("isMETIS", ( isMETIS ? 1 : 0 ) );
  h5w->write_intVector( "part", part_in, part_size );

  h5w->write_intScalar("isSerial", ( is_serial() ? 1 : 0 ) );

  delete h5w; H5Fclose(file_id);
}

void Global_Part_METIS::write_part_hdf5_64bit( const std::string &fileName,
    const int64_t * const &part_in,
    const int64_t &part_size, const int &cpu_size,
    const bool &part_isdual, const int &in_ncommon,
    const bool &isMETIS ) const
{
  std::string fName( fileName ); fName.append( ".h5" );

  hid_t file_id = H5Fcreate( fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  h5w->write_int64Scalar("part_size", part_size);
  h5w->write_intScalar("cpu_size", cpu_size);

  h5w->write_intScalar("part_isdual", ( part_isdual ? 1 : 0 ) );
  h5w->write_intScalar("in_ncommon", in_ncommon);

  h5w->write_intScalar("isMETIS", ( isMETIS ? 1 : 0 ) );
  h5w->write_int64Vector( "part", part_in, part_size );

  h5w->write_intScalar("isSerial", ( is_serial() ? 1 : 0 ) );

  delete h5w; H5Fclose(file_id);
}

// EOF
