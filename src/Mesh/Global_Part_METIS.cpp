#include "Global_Part_METIS.hpp"

Global_Part_METIS::Global_Part_METIS( const int &cpu_size,
        const int &in_ncommon, const bool &isDualGraph,
        const IMesh * const &mesh,
        const IIEN * const &IEN,
        const char * const &element_part_name,
        const char * const &node_part_name )
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

  idx_t nparts = cpu_size;
  idx_t objval;
  idx_t ncommon = in_ncommon;

  std::cout<<"-- allocating eptr and eind arrays... \n";

  const idx_t eptr_size = nElem + 1;
  const idx_t eind_size = nElem * nLocBas; 

  idx_t * eptr = new idx_t [eptr_size];
  
  std::cout<<"---- eptr allocated, with length "<<eptr_size<<" and size: ";
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

  time_tracker = clock();
  int metis_result;
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
  delete [] epart; delete [] npart;
  epart = nullptr; npart = nullptr;
  std::cout<<"-- Global_Part_METIS deleted. \n";
}

void Global_Part_METIS::write_part_hdf5( const char * const &fileName,
    const idx_t * const &part_in,
    const int &part_size, const int &cpu_size,
    const bool &part_isdual, const int &in_ncommon,
    const bool &isMETIS ) const
{
  std::string fName( fileName );
  fName.append(".h5");

  hid_t file_id, dataspace_id_n, dataspace_id_1;
  hsize_t dim_n[1], dim_1[1];

  // file creation
  file_id = H5Fcreate( fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // dataspace creation
  dim_n[0] = part_size; dim_1[0] = 1;

  dataspace_id_n = H5Screate_simple(1, dim_n, NULL);
  dataspace_id_1 = H5Screate_simple(1, dim_1, NULL);

  // dataset
  hid_t setid_part_size, setid_cpu_size, setid_part_isdual, setid_in_ncommon;
  hid_t setid_isMETIS, setid_part;

  setid_part_size = H5Dcreate( file_id, "part_size", H5T_NATIVE_INT, dataspace_id_1,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  setid_cpu_size = H5Dcreate( file_id, "cpu_size", H5T_NATIVE_INT, dataspace_id_1,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  setid_part_isdual = H5Dcreate( file_id, "part_isdual", H5T_NATIVE_INT, dataspace_id_1,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  setid_in_ncommon = H5Dcreate( file_id, "in_ncommon", H5T_NATIVE_INT, dataspace_id_1,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  setid_isMETIS = H5Dcreate( file_id, "isMETIS", H5T_NATIVE_INT, dataspace_id_1,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  setid_part = H5Dcreate( file_id, "part", H5T_NATIVE_INT, dataspace_id_n,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );


  H5Dwrite( setid_part_size, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &part_size);
  H5Dwrite( setid_cpu_size, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &cpu_size);

  int intbool;
  if(part_isdual)
    intbool = 1;
  else
    intbool = 0; 
  H5Dwrite( setid_part_isdual, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &intbool);

  H5Dwrite( setid_in_ncommon, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &in_ncommon);

  if(isMETIS)
    intbool = 1;
  else
    intbool = 0; 
  H5Dwrite( setid_isMETIS, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &intbool);

  H5Dwrite( setid_part, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      part_in );

  H5Dclose( setid_part_size );   H5Dclose( setid_cpu_size );
  H5Dclose( setid_part_isdual ); H5Dclose( setid_in_ncommon );
  H5Dclose( setid_isMETIS );     H5Dclose( setid_part );
  H5Sclose( dataspace_id_n );
  H5Sclose( dataspace_id_1 );
  H5Fclose(file_id);
}

void Global_Part_METIS::write_part_hdf5_64bit( const char * const &fileName,
    const int64_t * const &part_in,
    const int64_t &part_size, const int &cpu_size,
    const bool &part_isdual, const int &in_ncommon,
    const bool &isMETIS ) const
{
  std::string fName( fileName );
  fName.append( ".h5" );

  hid_t file_id;

  file_id = H5Fcreate( fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  h5w->write_int64Scalar("part_size", part_size);
  h5w->write_intScalar("cpu_size", cpu_size);

  int intbool;
  if(part_isdual)
    intbool = 1;
  else
    intbool = 0;

  h5w->write_intScalar("part_isdual", intbool);
  h5w->write_intScalar("in_ncommon", in_ncommon);

  if(isMETIS)
    intbool = 1;
  else
    intbool = 0;

  h5w->write_intScalar("isMETIS", intbool);

  h5w->write_int64Vector( "part", part_in, part_size );

  delete h5w;

  H5Fclose(file_id);
}

// EOF
