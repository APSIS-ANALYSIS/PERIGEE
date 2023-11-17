#include "AGlobal_Mesh_Info_FEM_3D.hpp"

AGlobal_Mesh_Info_FEM_3D::AGlobal_Mesh_Info_FEM_3D( 
    const std::string &fileBaseName,
    const int &cpu_rank )
{
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );
  
  std::vector<int> vdeg = h5r -> read_intVector("Global_Mesh_Info", "degree");
  
  xdegree = vdeg[0]; ydegree = vdeg[1]; zdegree = vdeg[2];

  nElem    = h5r -> read_intScalar("Global_Mesh_Info", "nElem");
  nFunc    = h5r -> read_intScalar("Global_Mesh_Info", "nFunc");
  nLocBas  = h5r -> read_intScalar("Global_Mesh_Info", "nLocBas");
  probDim  = h5r -> read_intScalar("Global_Mesh_Info", "probDim");
  elemType = h5r -> read_intScalar("Global_Mesh_Info", "elemType");

  delete h5r;
  H5Fclose( file_id );
}

AGlobal_Mesh_Info_FEM_3D::AGlobal_Mesh_Info_FEM_3D( 
    const std::string &fileBaseName,
    const int &cpu_rank, const int &num_enrich_node )
{
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );
  
  std::vector<int> vdeg = h5r -> read_intVector("Global_Mesh_Info", "degree");
  
  xdegree = vdeg[0]; ydegree = vdeg[1]; zdegree = vdeg[2];

  nElem    = h5r -> read_intScalar("Global_Mesh_Info", "nElem");
  nFunc    = h5r -> read_intScalar("Global_Mesh_Info", "nFunc");
  nLocBas  = h5r -> read_intScalar("Global_Mesh_Info", "nLocBas");
  probDim  = h5r -> read_intScalar("Global_Mesh_Info", "probDim");
  elemType = h5r -> read_intScalar("Global_Mesh_Info", "elemType");

  delete h5r; H5Fclose( file_id );
  
  SYS_T::print_fatal_if( num_enrich_node < 0, "Error: AGlobal_Mesh_Info_FEM_3D::AGlobal_Mesh_Info_FEM_3D input number of enrichment nodes is negative. \n");

  nFunc = nFunc + nElem * num_enrich_node;
  nLocBas = nLocBas + num_enrich_node;
  elemType = elemType + 10;
}

void AGlobal_Mesh_Info_FEM_3D::print_info() const
{
  std::cout<<"AGlobal_Mesh_Info_FEM_3D:"<<std::endl;
  std::cout<<"degree: "<<xdegree<<'\t'<<ydegree<<'\t'<<zdegree<<'\n';
  std::cout<<"nElem: "<<nElem<<'\n';
  std::cout<<"nFunc: "<<nFunc<<'\n';
  std::cout<<"nLocBas: "<<nLocBas<<std::endl;
  std::cout<<"probDim: "<<probDim<<std::endl;
  std::cout<<"elemType: "<<get_elemType()<<std::endl;
}

// EOF
