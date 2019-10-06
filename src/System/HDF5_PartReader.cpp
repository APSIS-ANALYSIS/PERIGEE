#include "HDF5_PartReader.hpp"

HDF5_PartReader::HDF5_PartReader( std::string partFileBaseName, const int &cpu_rank )
{
  // add sufix to base name
  std::string fName;
  fName.assign(partFileBaseName);
  fName.append("_p");

  if( cpu_rank / 10 == 0 )
    fName.append("0000");
  else if( cpu_rank / 100 == 0 )
    fName.append("000");
  else if( cpu_rank / 1000 == 0 )
    fName.append("00");
  else if( cpu_rank / 10000 == 0 )
    fName.append("0");

  std::stringstream sstrm;
  sstrm<<cpu_rank;
  fName.append(sstrm.str());
  fName.append(".h5");

  // open the file
  file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );  
}

HDF5_PartReader::~HDF5_PartReader()
{
  H5Fclose( file_id );
}

int HDF5_PartReader::get_intScalar( const char * group_name, const char * data_name ) const
{
  hid_t data_rank;
  hsize_t * data_dims;
  int * temp;

  get_intData(group_name, data_name, data_rank, data_dims, temp);

  if( data_rank != 1 )
  {
    std::cerr<<"ERROR: "<<data_name<<" in "<<group_name<<" is not a scalar. \n";
    exit(1);
  }
  if( data_dims[0] != 1 )
  {
    std::cerr<<"ERROR: "<<data_name<<" in "<<group_name<<" is not a scalar. \n";
    exit(1);
  }

  int value = temp[0];
  
  delete [] data_dims; delete [] temp;
  
  return value;
}

void HDF5_PartReader::get_intData( const char * group_name, const char * data_name,
           hid_t &data_rank, hsize_t * &data_dims, int * &data ) const
{
  // open data file in group
  hid_t group_id = H5Gopen( file_id, group_name, H5P_DEFAULT );
  hid_t data_id = H5Dopen( group_id, data_name, H5P_DEFAULT );

  // retrive dataspace of the dataset
  hid_t data_space = H5Dget_space( data_id );
  
  data_rank = H5Sget_simple_extent_ndims( data_space );
  data_dims = new hsize_t [data_rank];
  herr_t status = H5Sget_simple_extent_dims( data_space, data_dims, NULL );

  check_error(status, "get_intData");

  // deine the memory allocation for reading
  hid_t mem_space = H5Screate_simple(data_rank, data_dims, NULL);
  
  // setup the buffer
  hsize_t dSize = 1;
  for(hid_t ii=0; ii<data_rank; ++ii)
    dSize = dSize * data_dims[ii];

  data = new int [dSize];

  status = H5Dread( data_id, H5T_NATIVE_INT, mem_space, data_space, 
     H5P_DEFAULT, data );

  check_error(status, "get_intData");

  // clean out
  //delete [] data_dims; 
  H5Sclose( mem_space );
  H5Sclose( data_space );
  H5Dclose( data_id );
  H5Gclose( group_id );
}

void HDF5_PartReader::get_doubleData( const char * group_name, const char * data_name,
           hid_t &data_rank, hsize_t * &data_dims, double * &data ) const
{
  // open data file in group
  hid_t group_id = H5Gopen( file_id, group_name, H5P_DEFAULT );
  hid_t data_id = H5Dopen( group_id, data_name, H5P_DEFAULT );

  // retrive dataspace of the dataset
  hid_t data_space = H5Dget_space( data_id );
  
  data_rank = H5Sget_simple_extent_ndims( data_space );
  data_dims = new hsize_t [data_rank];
  herr_t status = H5Sget_simple_extent_dims( data_space, data_dims, NULL );

  check_error(status, "get_doubleData");

  // deine the memory allocation for reading
  hid_t mem_space = H5Screate_simple(data_rank, data_dims, NULL);
  
  // setup the buffer
  hsize_t dSize = 1;
  for(hid_t ii=0; ii<data_rank; ++ii)
    dSize = dSize * data_dims[ii];

  data = new double [dSize];

  status = H5Dread( data_id, H5T_NATIVE_DOUBLE, mem_space, data_space, 
     H5P_DEFAULT, data );

  check_error(status, "get_doubleData");
  
  // clean out
  //delete [] data_dims; 
  H5Sclose( mem_space );
  H5Sclose( data_space );
  H5Dclose( data_id );
  H5Gclose( group_id );
}

void HDF5_PartReader::get_intRowData( const char * group_name, const char * data_name,
    hid_t &data_rank, hsize_t * &data_dims, const int &row, int * &row_data ) const
{
  // open data file in group
  hid_t group_id = H5Gopen( file_id, group_name, H5P_DEFAULT );
  hid_t data_id  = H5Dopen( group_id, data_name, H5P_DEFAULT );

  // retrive dataspace of the data_id
  hid_t data_space = H5Dget_space( data_id );
  data_rank = H5Sget_simple_extent_ndims( data_space );

  if(data_rank != 2)
  {
    std::cerr<<"ERROR: "<<data_name<<" in "<<group_name<<" is not a two dimensional array. \n";
    exit(1);
  }

  data_dims = new hsize_t [data_rank];
  herr_t status = H5Sget_simple_extent_dims( data_space, data_dims, NULL );
  
  check_error(status, "get_intRowData");
  
  // setup hyperslab
  hsize_t offset[2]; offset[0] = row; offset[1] = 0;
  hsize_t count[2];  count[0]  = 1;   count[1]  = data_dims[1];
  hsize_t stride[2]; stride[0] = 1;   stride[1] = 1;
  hsize_t block[2];  block[0]  = 1;   block[1]  = 1;

  status = H5Sselect_hyperslab(data_space, H5S_SELECT_SET, offset, stride, count, block);
  
  check_error(status, "get_intRowData");
  
  // assign memory layout
  hid_t row_rank = 1;
  hsize_t row_dim[1]; 
  row_dim[0] = data_dims[1];
  hid_t row_space = H5Screate_simple(row_rank, row_dim, NULL);
 
  // assign row_data 
  row_data = new int [data_dims[1]];
  status = H5Dread( data_id, H5T_NATIVE_INT, row_space, data_space, H5P_DEFAULT, row_data);

  // clean up
  //delete [] data_dims;
  H5Sclose(row_space);
  H5Sclose(data_space);
  H5Gclose(group_id);
  H5Dclose(data_id);
}


void HDF5_PartReader::get_doubleRowData( const char * group_name, const char * data_name,
    hid_t &data_rank, hsize_t * &data_dims, const int &row, double * &row_data ) const
{
  // open data file in group
  hid_t group_id = H5Gopen( file_id, group_name, H5P_DEFAULT );
  hid_t data_id  = H5Dopen( group_id, data_name, H5P_DEFAULT );

  // retrive dataspace of the data_id
  hid_t data_space = H5Dget_space( data_id );
  data_rank = H5Sget_simple_extent_ndims( data_space );

  if(data_rank != 2)
  {
    std::cerr<<"ERROR: "<<data_name<<" in "<<group_name<<" is not a two dimensional array. \n";
    exit(1);
  }

  data_dims = new hsize_t [data_rank];
  herr_t status = H5Sget_simple_extent_dims( data_space, data_dims, NULL );
  
  check_error(status, "get_doubleRowData");
  
  // setup hyperslab
  hsize_t offset[2]; offset[0] = row; offset[1] = 0;
  hsize_t count[2];  count[0]  = 1;   count[1]  = data_dims[1];
  hsize_t stride[2]; stride[0] = 1;   stride[1] = 1;
  hsize_t block[2];  block[0]  = 1;   block[1]  = 1;

  status = H5Sselect_hyperslab(data_space, H5S_SELECT_SET, offset, stride, count, block);
  
  check_error(status, "get_doubleRowData");
  
  // assign memory layout
  hid_t row_rank = 1;
  hsize_t row_dim[1]; 
  row_dim[0] = data_dims[1];
  hid_t row_space = H5Screate_simple(row_rank, row_dim, NULL);
 
  // assign row_data 
  row_data = new double [data_dims[1]];
  status = H5Dread( data_id, H5T_NATIVE_DOUBLE, row_space, data_space, H5P_DEFAULT, row_data);

  // clean up
  //delete [] data_dims;
  H5Sclose(row_space);
  H5Sclose(data_space);
  H5Gclose(group_id);
  H5Dclose(data_id);
}

void HDF5_PartReader::get_GMI_degree( int &sDegree, int &tDegree, int &uDegree ) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "degree";
  hid_t data_rank; 
  hsize_t * data_dims;
  int * degree;
  
  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, degree);
  sDegree = degree[0];
  tDegree = degree[1];
  uDegree = degree[2];

  delete [] degree; delete [] data_dims;
}


void HDF5_PartReader::get_GMI_degree( int &sDegree, int &tDegree ) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "degree";
  hid_t data_rank; 
  hsize_t * data_dims;
  int * degree;
  
  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, degree);
  sDegree = degree[0];
  tDegree = degree[1];

  delete [] degree; delete [] data_dims;
}

void HDF5_PartReader::get_GMI_h_max( double &hx_max, double &hy_max, double &hz_max) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "h_max";
  hid_t data_rank;
  hsize_t * data_dims;
  double * hmax;

  get_doubleData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, hmax);
  hx_max = hmax[0];
  hy_max = hmax[1];
  hz_max = hmax[2];

  delete [] hmax; delete [] data_dims;
}

void HDF5_PartReader::get_GMI_h_max( double &hx_max, double &hy_max) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "h_max";
  hid_t data_rank;
  hsize_t * data_dims;
  double * hmax;

  get_doubleData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, hmax);
  hx_max = hmax[0];
  hy_max = hmax[1];

  delete [] hmax; delete [] data_dims;
}

void HDF5_PartReader::get_GMI_h_min( double &hx_min, double &hy_min, double &hz_min) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "h_min";
  hid_t data_rank;
  hsize_t * data_dims;
  double * hmin;

  get_doubleData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, hmin);

  hx_min = hmin[0];
  hy_min = hmin[1];
  hz_min = hmin[2];

  delete [] hmin; delete [] data_dims;
}

void HDF5_PartReader::get_GMI_h_min( double &hx_min, double &hy_min ) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "h_min";
  hid_t data_rank;
  hsize_t * data_dims;
  double * hmin;

  get_doubleData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, hmin);

  hx_min = hmin[0];
  hy_min = hmin[1];

  delete [] hmin; delete [] data_dims;
}

void HDF5_PartReader::get_GMI_nElem( int &nElem, int &nElem_x, int &nElem_y, int &nElem_z ) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "nElem";
  hid_t data_rank; 
  hsize_t * data_dims;
  int * temp;
  
  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, temp);

  nElem   = temp[0];
  nElem_x = temp[1];
  nElem_y = temp[2];
  nElem_z = temp[3];
  
  delete [] temp; delete [] data_dims;
}

void HDF5_PartReader::get_GMI_nElem( int &nElem, int &nElem_x, int &nElem_y ) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "nElem";
  hid_t data_rank; 
  hsize_t * data_dims;
  int * temp;
  
  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, temp);

  nElem   = temp[0];
  nElem_x = temp[1];
  nElem_y = temp[2];
  
  delete [] temp; delete [] data_dims;
}


void HDF5_PartReader::get_GMI_nElem( int &nElem ) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "nElem";
  hid_t data_rank; 
  hsize_t * data_dims;
  int * temp;
  
  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, temp);

  nElem   = temp[0];
  
  delete [] temp; delete [] data_dims;
}

void HDF5_PartReader::get_GMI_nFunc( int &nFunc, int &nFunc_x, int &nFunc_y, int &nFunc_z ) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "nFunc";
  hid_t data_rank; 
  hsize_t * data_dims;
  int * temp;
  
  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, temp);

  nFunc   = temp[0];
  nFunc_x = temp[1];
  nFunc_y = temp[2];
  nFunc_z = temp[3];

  delete [] temp; delete [] data_dims;
}

void HDF5_PartReader::get_GMI_nFunc( int &nFunc, int &nFunc_x, int &nFunc_y ) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "nFunc";
  hid_t data_rank; 
  hsize_t * data_dims;
  int * temp;
  
  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, temp);

  nFunc   = temp[0];
  nFunc_x = temp[1];
  nFunc_y = temp[2];

  delete [] temp; delete [] data_dims;
}


void HDF5_PartReader::get_GMI_nFunc( int &nFunc ) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "nFunc";
  hid_t data_rank; 
  hsize_t * data_dims;
  int * temp;
  
  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, temp);

  nFunc   = temp[0];

  delete [] temp; delete [] data_dims;
}


void HDF5_PartReader::get_GMI_nLocBas( int &nLocBas ) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "nLocBas";
  hid_t data_rank; 
  hsize_t * data_dims;
  int * temp;
  
  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, temp);

  nLocBas = temp[0];
  delete [] temp; delete [] data_dims;
}


void HDF5_PartReader::get_GMI_nLocBas( std::vector<int> &nLocBas ) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "nLocBas";
  hid_t data_rank; 
  hsize_t * data_dims;
  int * temp;
  
  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, temp);

  const int len = data_dims[0];

  VEC_T::fillArray(nLocBas, temp, len);
  VEC_T::shrink2fit(nLocBas);

  delete [] temp; delete [] data_dims;
}


void HDF5_PartReader::get_GMI_probDim( int &probDim ) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "probDim";
  hid_t data_rank;
  hsize_t * data_dims;
  int * temp;

  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, temp);

  probDim = temp[0];
  delete [] temp; delete [] data_dims;
}

void HDF5_PartReader::get_GMI_dofNum( int &dofNum ) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "dofNum";
  hid_t data_rank;
  hsize_t * data_dims;
  int * temp;

  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, temp);

  dofNum = temp[0];
  delete [] temp; delete [] data_dims;
}

void HDF5_PartReader::get_GMI_dofMat( int &dofMat ) const
{
  std::string group_name = "Global_Mesh_Info";
  std::string data_name  = "dofMat";
 
  // If the preprocessor does not specify the dofMat data, read dofNum instead 
  if( !H5Lexists(file_id, "/Global_Mesh_Info/dofMat", H5P_DEFAULT) ) 
    data_name = "dofNum";

  hid_t data_rank;
  hsize_t * data_dims;
  int * temp;
  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, temp);
  dofMat = temp[0];
  delete [] temp; delete [] data_dims;
}


void HDF5_PartReader::get_EXT_elemType( int &elemType ) const
{
  std::string group_name = "Extraction";
  std::string data_name  = "elemType";
  hid_t data_rank;
  hsize_t * data_dims;
  int * temp;

  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, temp);

  elemType = temp[0];
  delete [] temp; delete [] data_dims;
}



void HDF5_PartReader::get_LIEN( const int &e, int * &LIEN, int &length ) const
{
  std::string group_name = "LIEN";
  std::string data_name = "LIEN";

  hid_t data_rank; hsize_t * data_dims;

  get_intRowData(group_name.c_str(), data_name.c_str(), data_rank,
      data_dims, e, LIEN );

  length = data_dims[1]; delete [] data_dims;
}

void HDF5_PartReader::get_LIEN( std::vector<int> &LIEN ) const
{
  hid_t data_rank;
  hsize_t * data_dims;
  int * temp;

  get_intData("LIEN", "LIEN", data_rank, data_dims, temp );

  if(data_rank != 2)
  {
    SYS_T::commPrint("Error: LIEN is not a two-dimensional array on disk. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  int len = data_dims[0] * data_dims[1];

  VEC_T::fillArray(LIEN, temp, len);
  VEC_T::shrink2fit(LIEN);

  delete [] temp; delete [] data_dims;
}


void HDF5_PartReader::get_1D_LIEN( std::vector<int> &LIEN ) const
{
  hid_t data_rank;
  hsize_t * data_dims;
  int * temp;

  get_intData("LIEN", "LIEN", data_rank, data_dims, temp);

  const int len = data_dims[0];

  VEC_T::fillArray(LIEN, temp, len);
  VEC_T::shrink2fit(LIEN);

  delete [] temp; delete [] data_dims;
}



void HDF5_PartReader::get_LE( int * &elem_loc, int &nlocalele ) const
{
  std::string group_name = "Local_Elem";
  std::string data_name = "elem_loc";

  hid_t data_rank; hsize_t * data_dims;

  get_intData(group_name.c_str(), data_name.c_str(), data_rank, data_dims, elem_loc);

  std::string data_name_2 = "nlocalele";
  hid_t data_rank_2; hsize_t * data_dims_2;
  int * temp;

  get_intData(group_name.c_str(), data_name_2.c_str(), data_rank_2, data_dims_2, temp );

  nlocalele = temp[0];

  delete [] temp; delete [] data_dims; delete [] data_dims_2;
}


void HDF5_PartReader::get_L2GN( int &nlocghonode, std::vector<int> &local_to_global ) const
{
  const std::string group_name = "Local_Node";
  nlocghonode = get_intScalar(group_name.c_str(), "nlocghonode");

  int * l2g_temp;

  hid_t drank_1; hsize_t * ddims_1;
  get_intData(group_name.c_str(), "local_to_global", drank_1, ddims_1, l2g_temp);

  if(drank_1 != 1)
  {
    std::cerr<<"ERROR: local_to_global is not in correct format. \n";
    exit(1);
  }
  if(int(ddims_1[0]) != nlocghonode)
  {
    std::cerr<<"ERROR: local_to_global is not in correct format. \n";
    exit(1);
  }
  delete [] ddims_1;

  VEC_T::fillArray(local_to_global, l2g_temp, nlocghonode);

  delete [] l2g_temp; l2g_temp = NULL;
}


void HDF5_PartReader::get_LN( int &nlocalnode, int &nghostnode, int &nbadnode, 
    int &nlocghonode, int &ntotalnode, int * &local_to_global, 
    int * &node_ghost, int * &node_loc, int * &node_loc_original ) const
{
  std::string group_name = "Local_Node";
  nlocalnode = get_intScalar(group_name.c_str(), "nlocalnode");
  nghostnode = get_intScalar(group_name.c_str(), "nghostnode");
  nbadnode = get_intScalar(group_name.c_str(), "nbadnode");
  nlocghonode = get_intScalar(group_name.c_str(), "nlocghonode");
  ntotalnode = get_intScalar(group_name.c_str(), "ntotalnode");

  // local_to_global
  hid_t drank_1; hsize_t * ddims_1;
  get_intData(group_name.c_str(), "local_to_global", drank_1, ddims_1, local_to_global);

  if(drank_1 != 1)
  {
    std::cerr<<"ERROR: local_to_global is not in correct format. \n";
    exit(1);
  }
  if(int(ddims_1[0]) != nlocghonode)
  {
    std::cerr<<"ERROR: local_to_global is not in correct format. \n";
    exit(1);
  }
  delete [] ddims_1;

  // node_ghost
  if(nghostnode > 0)
  {
    hid_t drank_2; hsize_t * ddims_2;
    get_intData(group_name.c_str(), "node_ghost", drank_2, ddims_2, node_ghost);

    if(drank_2 != 1)
    {
      std::cerr<<"ERROR: node_ghost is not in correct format. \n";
      exit(1);
    }
    if( int(ddims_2[0]) != nghostnode )
    {
      std::cerr<<"ERROR: node_ghost is not in correct format. \n";
      exit(1);
    }
    delete [] ddims_2;
  }
  else
    node_ghost = NULL;

  // node_loc
  hid_t drank_3; hsize_t * ddims_3;
  get_intData(group_name.c_str(), "node_loc", drank_3, ddims_3, node_loc);
  if(drank_3 != 1)
  {
    std::cerr<<"ERROR: node_loc is not in correct format. \n";
    exit(1);
  }
  if( int(ddims_3[0]) != nlocalnode )
  {
    std::cerr<<"ERROR: node_loc is not in correct format. \n";
    exit(1);
  }
  delete [] ddims_3;

  // node_loc_original
  hid_t drank_4; hsize_t * ddims_4;
  get_intData(group_name.c_str(), "node_loc_original", drank_4, ddims_4,
      node_loc_original);
  if(drank_4 != 1)
  {
    std::cerr<<"ERROR: node_loc_original is not in correct format. \n";
    exit(1);
  }
  if( int(ddims_4[0]) != nlocalnode )
  {
    std::cerr<<"ERROR: node_loc_original is not in correct format. \n";
    exit(1);
  }
  delete [] ddims_4;
}

void HDF5_PartReader::get_PI( int &cpu_rank, int &cpu_size, 
    int &dual_edge_ncommon ) const
{
  std::string group_name = "Part_Info";
  cpu_rank = get_intScalar(group_name.c_str(), "cpu_rank");
  cpu_size = get_intScalar(group_name.c_str(), "cpu_size");
  dual_edge_ncommon = get_intScalar(group_name.c_str(), "dual_edge_ncommon");
}

void HDF5_PartReader::get_CPL( double * &ctrlPts_x_loc, double * &ctrlPts_y_loc,
    double * &ctrlPts_z_loc, double * &ctrlPts_w_loc, int &len ) const
{
  const std::string group_name = "ctrlPts_loc";

  // ctrlPts_x_loc
  hid_t rank_x; hsize_t * dims_x;
  get_doubleData(group_name.c_str(), "ctrlPts_x_loc", rank_x, dims_x, ctrlPts_x_loc);
  if(rank_x != 1)
  {
    std::cerr<<"ERROR: ctrlPts_x_loc is in wrong format. \n";
    exit(1);
  }
  len = int(dims_x[0]); delete [] dims_x;

  // ctrlPts_y_loc
  hid_t rank_y; hsize_t * dims_y;
  get_doubleData(group_name.c_str(), "ctrlPts_y_loc", rank_y, dims_y, ctrlPts_y_loc);
  if(rank_y != 1)
  {
    std::cerr<<"ERROR: ctrlPts_y_loc is in wrong format. \n";
    exit(1);
  }
  if( len != int(dims_y[0]) )
  {
    std::cerr<<"ERROR: ctrlPts_y_loc is in wrong format. \n";
    exit(1);
  }
  delete [] dims_y;

  // ctrlPts_z_loc
  hid_t rank_z; hsize_t * dims_z;
  get_doubleData(group_name.c_str(), "ctrlPts_z_loc", rank_z, dims_z, ctrlPts_z_loc);
  if(rank_z != 1)
  {
    std::cerr<<"ERROR: ctrlPts_z_loc is in wrong format. \n";
    exit(1);
  }
  if( len != int(dims_z[0]) )
  {
    std::cerr<<"ERROR: ctrlPts_z_loc is in wrong format. \n";
    exit(1);
  }
  delete [] dims_z;

  // ctrlPts_w_loc
  hid_t rank_w; hsize_t * dims_w;
  get_doubleData(group_name.c_str(), "ctrlPts_w_loc", rank_w, dims_w, ctrlPts_w_loc);
  if(rank_w != 1)
  {
    std::cerr<<"ERROR: ctrlPts_w_loc is in wrong format. \n";
    exit(1);
  }
  if( len != int(dims_w[0]) )
  {
    std::cerr<<"ERROR: ctrlPts_w_loc is in wrong format. \n";
    exit(1);
  }
  delete [] dims_w;
}

void HDF5_PartReader::get_CPL( std::vector<double> &ctrlPts_x_loc,
    std::vector<double> &ctrlPts_y_loc, std::vector<double> &ctrlPts_z_loc,
    std::vector<double> &ctrlPts_w_loc ) const
{
  double * ctrl_x; double * ctrl_y; double * ctrl_z; double * ctrl_w;
  int len;
  get_CPL(ctrl_x, ctrl_y, ctrl_z, ctrl_w, len);

  VEC_T::fillArray(ctrlPts_x_loc, ctrl_x, len);
  VEC_T::fillArray(ctrlPts_y_loc, ctrl_y, len);
  VEC_T::fillArray(ctrlPts_z_loc, ctrl_z, len);
  VEC_T::fillArray(ctrlPts_w_loc, ctrl_w, len);

  delete [] ctrl_x; delete [] ctrl_y; delete [] ctrl_z; delete [] ctrl_w;
}

void HDF5_PartReader::get_BC_LID_dof( std::vector<int> &LID, int &dof ) const
{
  int nlocghonode = get_intScalar("Local_Node", "nlocghonode");

  int * lid_array; hid_t rank_lid; hsize_t * dims_lid;
  get_intData( "bc", "LID", rank_lid, dims_lid, lid_array );

  if(rank_lid != 1)
  {
    std::cerr<<"ERROR: LID array is not in correct dimension. \n";
    exit(1);
  }
  if( dims_lid[0] % nlocghonode != 0 )
  {
    std::cerr<<"ERROR: LID length is not dof * nlocghonode. \n";
    exit(1);
  }

  dof = dims_lid[0] / nlocghonode;

  VEC_T::fillArray(LID, lid_array, dims_lid[0]);

  delete [] lid_array; delete [] dims_lid;
}

void HDF5_PartReader::get_BC_LD( std::vector<int> &LDN, std::vector<int> &Num_LD ) const
{
  int * num_array; hid_t rank_num; hsize_t * dims_num;
  LDN.clear();

  get_intData("bc", "Num_LD", rank_num, dims_num, num_array);
  if( rank_num != 1)
  {
    std::cerr<<"ERROR: Num_LD is in wrong dimension. \n";
    exit(1);
  }
  VEC_T::fillArray(Num_LD, num_array, dims_num[0]);
  delete [] num_array; delete [] dims_num;

  int size_ldn = 0;
  for(unsigned int ii=0; ii<Num_LD.size(); ++ii)
    size_ldn += Num_LD[ii];

  if(size_ldn > 0)
  {
    int * ldn_array; hid_t rank_ldn; hsize_t * dims_ldn;
    get_intData("bc", "LDN", rank_ldn, dims_ldn, ldn_array);
    if(rank_ldn != 1)
    {
      std::cerr<<"ERROR: LDN is ont in one dimension. \n";
      exit(1);
    }
    VEC_T::fillArray(LDN, ldn_array, dims_ldn[0]);
    delete [] ldn_array; delete [] dims_ldn;
  }
}

void HDF5_PartReader::get_BC_LP(std::vector<int> &LPSN, std::vector<int> &LPMN,
    std::vector<int> &Num_LP ) const
{
  int * num_array; hid_t rank_num; hsize_t * dims_num;
  LPSN.clear(); LPMN.clear();

  get_intData("bc", "Num_LP", rank_num, dims_num, num_array);
  if(rank_num != 1)
  {
    std::cerr<<"ERROR: Num_LP is in wrong dimension. \n";
    exit(1);
  }
  VEC_T::fillArray(Num_LP, num_array, dims_num[0]);
  delete [] dims_num; delete [] num_array;

  int size_lp = 0;
  for(unsigned int ii=0; ii<Num_LP.size(); ++ii)
    size_lp += Num_LP[ii];

  if(size_lp > 0)
  {
    int * lpsn; int * lpmn;
    hid_t rank_s; hid_t rank_m;
    hsize_t * dims_s; hsize_t * dims_m;
    get_intData("bc", "LPSN", rank_s, dims_s, lpsn);
    get_intData("bc", "LPMN", rank_m, dims_m, lpmn);
    if(rank_s != 1)
    {
      std::cerr<<"ERROR: LPSN is not in 1-dim array. \n";
      exit(1);
    }
    if(rank_m != 1)
    {
      std::cerr<<"ERROR: LPMN is not in 1-dim array. \n";
      exit(1);
    }
    VEC_T::fillArray(LPSN, lpsn, dims_s[0]);
    VEC_T::fillArray(LPMN, lpmn, dims_m[0]);
    delete [] lpsn; delete [] lpmn; delete [] dims_s; delete [] dims_m;
  }
}

void HDF5_PartReader::get_BC_LBCE( std::vector<int> &Num_LBCElem,
    std::vector<int> &LFront_Elem, std::vector<int> &LBack_Elem,
    std::vector<int> &LLeft_Elem,  std::vector<int> &LRight_Elem, 
    std::vector<int> &LTop_Elem, std::vector<int> &LBottom_Elem ) const
{
  int * num; hid_t rank_num; hsize_t * dims_num;
  LFront_Elem.clear(); LBack_Elem.clear(); LLeft_Elem.clear();
  LRight_Elem.clear(); LTop_Elem.clear(); LBottom_Elem.clear();

  get_intData("bc", "Num_LBCElem", rank_num, dims_num, num);
  if(rank_num != 1)
  {
    std::cerr<<"ERROR: Num_LBCElem is in wrong dimension. \n";
    exit(1);
  }
  int dof;
  if(dims_num[0] % 6 != 0)
  {
    std::cerr<<"ERROR: Num_LBCElem is in wrong size. \n";
    exit(1);
  }
  dof = dims_num[0] / 6;
  VEC_T::fillArray(Num_LBCElem, num, dims_num[0]);
  delete [] num; delete [] dims_num;

  int num_front = 0; int num_back = 0;
  int num_left = 0; int num_right = 0;
  int num_top = 0; int num_bottom = 0;
  for(int ii=0; ii<dof; ++ii)
  {
    num_front  += Num_LBCElem[6*ii + 0];
    num_back   += Num_LBCElem[6*ii + 1];
    num_left   += Num_LBCElem[6*ii + 2];
    num_right  += Num_LBCElem[6*ii + 3];
    num_top    += Num_LBCElem[6*ii + 4];
    num_bottom += Num_LBCElem[6*ii + 5];
  }

  if(num_front > 0)
  {
    int * front; hid_t rank_front; hsize_t * dims_front;
    get_intData("bc","LFront_Elem", rank_front, dims_front, front);
    assert(rank_front == 1); assert(int(dims_front[0]) == num_front);
    VEC_T::fillArray(LFront_Elem, front, dims_front[0]);
    delete [] front; delete [] dims_front;
  }

  if(num_back > 0)
  {
    int * back; hid_t rank_back; hsize_t * dims_back;
    get_intData("bc", "LBack_Elem", rank_back, dims_back, back);
    assert(rank_back == 1); assert(int(dims_back[0]) == num_back);
    VEC_T::fillArray(LBack_Elem, back, dims_back[0]);
    delete [] back; delete [] dims_back;
  }

  if(num_left > 0)
  {
    int * left; hid_t rank_left; hsize_t * dims_left;
    get_intData("bc", "LLeft_Elem", rank_left, dims_left, left);
    assert(rank_left == 1); assert(int(dims_left[0]) == num_left);
    VEC_T::fillArray(LLeft_Elem, left, dims_left[0]);
    delete [] left; delete [] dims_left;
  }

  if(num_right >0)
  {
    int * right; hid_t rank_right; hsize_t * dims_right;
    get_intData("bc", "LRight_Elem", rank_right, dims_right, right);
    assert(rank_right == 1); assert(int(dims_right[0]) == num_right);
    VEC_T::fillArray(LRight_Elem, right, dims_right[0]);
    delete [] right; delete [] dims_right;
  }

  if(num_top >0)
  {
    int * top; hid_t rank_top; hsize_t * dims_top;
    get_intData("bc", "LTop_Elem", rank_top, dims_top, top);
    assert(rank_top == 1); assert(int(dims_top[0]) == num_top);
    VEC_T::fillArray(LTop_Elem, top, dims_top[0]);
    delete [] top; delete [] dims_top;
  }

  if(num_bottom > 0)
  {
    int * bottom; hid_t rank_bottom; hsize_t * dims_bottom;
    get_intData("bc", "LBottom_Elem", rank_bottom, dims_bottom, bottom);
    assert(rank_bottom == 1); assert(int(dims_bottom[0]) == num_bottom);
    VEC_T::fillArray(LBottom_Elem, bottom, dims_bottom[0]);
    delete [] bottom; delete [] dims_bottom;
  }
}


void HDF5_PartReader::get_NBC_LID_dof( std::vector<int> &LID, int &dof ) const
{
  const int nlocghonode = get_intScalar("Local_Node", "nlocghonode");

  int * lid_array; hid_t rank_lid; hsize_t * dims_lid;
  get_intData( "nbc", "LID", rank_lid, dims_lid, lid_array );

  if(rank_lid != 1)
  {
    SYS_T::commPrint("Error: LID array is not in correct dimension. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  if( dims_lid[0] % nlocghonode != 0 )
  {
    SYS_T::commPrint("Error: LID array length is not dof * nlocghonode. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  dof = dims_lid[0] / nlocghonode;

  VEC_T::fillArray(LID, lid_array, dims_lid[0]);

  delete [] lid_array; delete [] dims_lid;
}


void HDF5_PartReader::get_NBC_LD( std::vector<int> &LDN, std::vector<int> &Num_LD ) const
{
  int * num_array; hid_t rank_num; hsize_t * dims_num;

  get_intData("nbc", "Num_LD", rank_num, dims_num, num_array);

  if( rank_num != 1 )
  {
    SYS_T::commPrint("Error: Num_LD is in wrong dimension. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  VEC_T::fillArray(Num_LD, num_array, dims_num[0]);
  delete [] num_array; delete [] dims_num;

  int size_ldn = 0;
  for(unsigned int ii=0; ii<Num_LD.size(); ++ii)
    size_ldn += Num_LD[ii];

  LDN.clear();

  if(size_ldn >0)
  {
    int * ldn_array; hid_t rank_ldn; hsize_t * dims_ldn;
    get_intData("nbc", "LDN", rank_ldn, dims_ldn, ldn_array);
    if(rank_ldn != 1)
    {
      SYS_T::commPrint("Error: LDN is not in one dimension. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    if(dims_ldn[0] != hsize_t(size_ldn))
    {
      SYS_T::commPrint("Error: LDN length does not match Num_LD. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    VEC_T::fillArray(LDN, ldn_array, dims_ldn[0]);
    delete [] ldn_array; delete [] dims_ldn;
  }
}


void HDF5_PartReader::get_NBC_LPS(std::vector<int> &LPSN, std::vector<int> &LPMN, 
    std::vector<int> &Num_LPS) const
{
  int * num_array; hid_t rank_num; hsize_t * dims_num;
  LPSN.clear(); LPMN.clear();

  get_intData("nbc", "Num_LPS", rank_num, dims_num, num_array);

  if(rank_num != 1)
  {
    SYS_T::commPrint("Error: Num_LPS is not one-dimensional. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  VEC_T::fillArray(Num_LPS, num_array, dims_num[0]);
  delete [] dims_num; delete [] num_array;

  int size_lp = 0;
  for(unsigned int ii=0; ii<Num_LPS.size(); ++ii)
    size_lp += Num_LPS[ii];

  if(size_lp > 0)
  {
    int * lpsn; int * lpmn;
    hid_t rank_s; hid_t rank_m;
    hsize_t * dims_s; hsize_t * dims_m;
    get_intData("nbc", "LPSN", rank_s, dims_s, lpsn);
    get_intData("nbc", "LPMN", rank_m, dims_m, lpmn);
    if(rank_s != 1)
    {
      SYS_T::commPrint("Error: LPSN is not 1-dimensional. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    if(rank_m != 1)
    {
      SYS_T::commPrint("Error: LPMN is not 1-dimensional. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    if(dims_s[0] != hsize_t(size_lp) || dims_m[0] != hsize_t(size_lp))
    {
      SYS_T::commPrint("Error: LPMN or LPMN length does not match Num_LPS. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    VEC_T::fillArray(LPSN, lpsn, dims_s[0]);
    VEC_T::fillArray(LPMN, lpmn, dims_m[0]);
    delete [] lpsn; delete [] lpmn;
    delete [] dims_s; delete [] dims_m;
  }
}


void HDF5_PartReader::get_NBC_LPM(std::vector<int> &LocalMaster, 
    std::vector<int> &LocalMasterSlave, std::vector<int> &Num_LPM) const
{
  int * num_array; hid_t rank_num; hsize_t * dims_num;
  LocalMaster.clear(); LocalMasterSlave.clear();

  get_intData("nbc", "Num_LPM", rank_num, dims_num, num_array);

  if(rank_num != 1)
  {
    SYS_T::commPrint("Error: Num_LPM is not 1-dimensional. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  VEC_T::fillArray(Num_LPM, num_array, dims_num[0]);
  delete [] num_array; delete [] dims_num;

  unsigned int size_lp = 0;
  for(unsigned int ii=0; ii<Num_LPM.size(); ++ii)
    size_lp += Num_LPM[ii];

  if(size_lp > 0)
  {
    int * lm; int * lms;
    hid_t rank_m; hid_t rank_s;
    hsize_t * dims_m; hsize_t * dims_s;
    get_intData("nbc", "LocalMaster", rank_m, dims_m, lm);
    get_intData("nbc", "LocalMasterSlave", rank_s, dims_s, lms);

    if(rank_s != 1 || rank_m != 1)
    {
      SYS_T::commPrint("Error: LocalMaster or LocalMasterSlave is not 1-dimensional. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    if(dims_m[0] != size_lp || dims_s[0] != size_lp )
    {
      SYS_T::commPrint("Error: LocalMaster or LocalMasterSlave does not match Num_LPM. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    VEC_T::fillArray(LocalMaster, lm, dims_m[0]);
    VEC_T::fillArray(LocalMasterSlave, lms, dims_s[0]);
    delete [] lm; delete [] lms;
    delete [] dims_m; delete [] dims_s;
  }
}



void HDF5_PartReader::get_EBC( std::vector<int> &num_lbcelem, std::vector<int> &lfront_elem,
    std::vector<int> &lback_elem, std::vector<int> &lleft_elem,
    std::vector<int> &lright_elem, std::vector<int> &ltop_elem,
    std::vector<int> &lbottom_elem  ) const
{
  num_lbcelem.clear();
  lfront_elem.clear(); lback_elem.clear();
  lleft_elem.clear(); lright_elem.clear();
  ltop_elem.clear(); lbottom_elem.clear();

  int * num; hid_t rank_num; hsize_t * dims_num;

  get_intData("ebc", "Num_LBCElem", rank_num, dims_num, num);
  if(rank_num != 1 || dims_num[0] != 6)
  {
    SYS_T::commPrint("Error: Num_LBCElem is not in correct format. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  VEC_T::fillArray(num_lbcelem, num, dims_num[0]);
  delete [] num; delete [] dims_num;

  if(num_lbcelem[0] > 0)
  {
    int * front; hid_t rank_front; hsize_t * dims_front;
    get_intData("ebc", "LFront_Elem", rank_front, dims_front, front);
    if(rank_front != 1 || (int) dims_front[0] != num_lbcelem[0])
    {
      SYS_T::commPrint("Error: ebc LFront_Elem is not in correct format. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    VEC_T::fillArray(lfront_elem, front, dims_front[0]);
    delete [] front; delete [] dims_front;
  }

  if(num_lbcelem[1] > 0)
  {
    int * back; hid_t rank_back; hsize_t * dims_back;
    get_intData("ebc","LBack_Elem", rank_back, dims_back, back);
    if(rank_back != 1 || (int) dims_back[0] != num_lbcelem[1])
    {
      SYS_T::commPrint("Error: ebc LBack_Elem is not in correct format. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    VEC_T::fillArray(lback_elem, back, dims_back[0]);
    delete [] back; delete [] dims_back;
  }

  if(num_lbcelem[2] > 0)
  {
    int * left; hid_t rank_left; hsize_t * dims_left;
    get_intData("ebc", "LLeft_Elem", rank_left, dims_left, left);
    if(rank_left != 1 || (int) dims_left[0] != num_lbcelem[2])
    {
      SYS_T::commPrint("Error: ebc LLeft_Elem is not in correct format. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    VEC_T::fillArray(lleft_elem, left, dims_left[0]);
    delete [] left; delete [] dims_left;
  }

  if(num_lbcelem[3] > 0)
  {
    int * right; hid_t rank_right; hsize_t * dims_right;
    get_intData("ebc","LRight_Elem", rank_right, dims_right, right);
    if(rank_right != 1 || (int) dims_right[0] != num_lbcelem[3])
    {
      SYS_T::commPrint("Error: ebc LRight_Elem is not in correct format. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    VEC_T::fillArray(lright_elem, right, dims_right[0]);
    delete [] right; delete [] dims_right;
  }

  if(num_lbcelem[4] > 0)
  {
    int * top; hid_t rank_top; hsize_t * dims_top;
    get_intData("ebc", "LTop_Elem", rank_top, dims_top, top);
    if(rank_top != 1 || (int) dims_top[0] != num_lbcelem[4])
    {
      SYS_T::commPrint("Error: ebc LTop_Elem is not in correct format. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    VEC_T::fillArray(ltop_elem, top, dims_top[0]);
    delete [] top; delete [] dims_top;
  }

  if(num_lbcelem[5] > 0)
  {
    int * bottom; hid_t rank_bottom; hsize_t * dims_bottom;
    get_intData("ebc", "LBottom_Elem", rank_bottom, dims_bottom, bottom);
    if(rank_bottom != 1 || (int) dims_bottom[0] != num_lbcelem[5])
    {
      SYS_T::commPrint("Error: ebc LBottom_Elem is not in correct format. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    VEC_T::fillArray(lbottom_elem, bottom, dims_bottom[0]);
    delete [] bottom; delete [] dims_bottom;
  }
}



void HDF5_PartReader::get_EXT_x( std::vector<double> &ext_x ) const
{
  std::string group_name = "Extraction";
  std::string data_name  = "extractor_x";
  hid_t data_rank; hsize_t * data_dims;

  double * temp;

  get_doubleData(group_name.c_str(), data_name.c_str(), data_rank, 
      data_dims, temp );

  int len = data_dims[0] * data_dims[1];

  VEC_T::fillArray(ext_x, temp, len);

  delete [] temp; delete [] data_dims;
}


void HDF5_PartReader::get_EXT_y( std::vector<double> &ext_y ) const
{
  std::string group_name = "Extraction";
  std::string data_name  = "extractor_y";
  hid_t data_rank; hsize_t * data_dims;

  double * temp;

  get_doubleData(group_name.c_str(), data_name.c_str(), data_rank, 
      data_dims, temp );

  int len = data_dims[0] * data_dims[1];

  VEC_T::fillArray(ext_y, temp, len);

  delete [] temp; delete [] data_dims;
}


void HDF5_PartReader::get_EXT_z( std::vector<double> &ext_z ) const
{
  std::string group_name = "Extraction";
  std::string data_name  = "extractor_z";
  hid_t data_rank; hsize_t * data_dims;

  double * temp;

  get_doubleData(group_name.c_str(), data_name.c_str(), data_rank, 
      data_dims, temp );

  int len = data_dims[0] * data_dims[1];

  VEC_T::fillArray(ext_z, temp, len);

  delete [] temp; delete [] data_dims;
}


void HDF5_PartReader::get_EXT_TS_full( std::vector<double> &ext ) const
{
  const std::string group_name = "Extraction";
  const std::string data_name  = "extractor";
  hid_t data_rank;
  hsize_t * data_dims;
  double * temp;

  get_doubleData(group_name.c_str(), data_name.c_str(), data_rank,
      data_dims, temp );

  const int len = data_dims[0];

  VEC_T::fillArray(ext, temp, len);

  delete [] temp; delete [] data_dims;
}



void HDF5_PartReader::get_LIEN( const int &e, std::vector<int> &LIEN ) const
{
  int * temp; int len;
  get_LIEN(e, temp, len);
  VEC_T::fillArray(LIEN, temp, len);
  delete [] temp;
}

void HDF5_PartReader::get_LE( std::vector<int> &elem_loc, int &nlocalele ) const
{
  int * temp;
  get_LE(temp, nlocalele);
  VEC_T::fillArray(elem_loc, temp, nlocalele);
  delete [] temp;
}

void HDF5_PartReader::get_LN( int &nlocalnode, int &nghostnode, int &nbadnode,
    int &nlocghonode, int &ntotalnode, std::vector<int> &local_to_global,
    std::vector<int> &node_ghost, std::vector<int> &node_loc,
    std::vector<int> &node_loc_original ) const
{
  int * temp_ltg; int * temp_ghost; int * temp_loc; int * temp_original;

  get_LN(nlocalnode, nghostnode, nbadnode, nlocghonode, ntotalnode, temp_ltg,
      temp_ghost, temp_loc, temp_original );

  VEC_T::fillArray(local_to_global, temp_ltg, nlocghonode);

  if( nghostnode == 0 )
    node_ghost.clear();
  else
    VEC_T::fillArray(node_ghost, temp_ghost, nghostnode);

  VEC_T::fillArray(node_loc, temp_loc, nlocalnode);
  VEC_T::fillArray(node_loc_original, temp_original, nlocalnode);

  delete [] temp_ltg; delete [] temp_ghost; delete [] temp_loc;
  delete [] temp_original;
}

void HDF5_PartReader::get_hxyz( std::vector<double> &hx, std::vector<double> &hy,
    std::vector<double> &hz ) const
{
  double * temp_hx; double * temp_hy; double * temp_hz;
  hid_t temp_rank; 
  hsize_t * temp_dims_x;
  hsize_t * temp_dims_y;
  hsize_t * temp_dims_z;
  get_doubleData("Mesh_Size", "hx", temp_rank, temp_dims_x, temp_hx );
  get_doubleData("Mesh_Size", "hy", temp_rank, temp_dims_y, temp_hy );
  get_doubleData("Mesh_Size", "hz", temp_rank, temp_dims_z, temp_hz );
  VEC_T::fillArray(hx, temp_hx, temp_dims_x[0]);
  VEC_T::fillArray(hy, temp_hy, temp_dims_y[0]);
  VEC_T::fillArray(hz, temp_hz, temp_dims_z[0]);

  delete [] temp_hx; delete [] temp_hy; delete [] temp_hz;
  delete [] temp_dims_x; delete [] temp_dims_y; delete [] temp_dims_z;
}


void HDF5_PartReader::get_hxy( std::vector<double> &hx, std::vector<double> &hy
    ) const
{
  double * temp_hx; double * temp_hy;
  hid_t temp_rank; 
  hsize_t * temp_dims_x;
  hsize_t * temp_dims_y;
  get_doubleData("Mesh_Size", "hx", temp_rank, temp_dims_x, temp_hx );
  get_doubleData("Mesh_Size", "hy", temp_rank, temp_dims_y, temp_hy );
  VEC_T::fillArray(hx, temp_hx, temp_dims_x[0]);
  VEC_T::fillArray(hy, temp_hy, temp_dims_y[0]);

  delete [] temp_hx; delete [] temp_hy;
  delete [] temp_dims_x; delete [] temp_dims_y;
}

// EOF
