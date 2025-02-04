#include "ALocal_IEN.hpp"

ALocal_IEN::ALocal_IEN( const std::string &fileBaseName, const int &cpu_rank )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  auto h5r = SYS_T::make_unique<HDF5_Reader>(file_id);

  nlocalele = h5r -> read_intScalar("Local_Elem", "nlocalele");

  nLocBas = h5r -> read_intScalar("Global_Mesh_Info", "nLocBas");

  int num_row, num_col;
  LIEN = h5r -> read_intMatrix("LIEN", "LIEN", num_row, num_col);
  
  H5Fclose( file_id );

  SYS_T::print_fatal_if( num_row != nlocalele, "Error: ALocal_IEN::LIEN size does not match the number of element. \n");

  SYS_T::print_fatal_if( num_col != nLocBas, "Error: ALocal_IEN::LIEN size does not match the value of nLocBas. \n");
}

ALocal_IEN::ALocal_IEN( const HDF5_Reader * const &h5r ) 
{
  nlocalele = h5r -> read_intScalar("Local_Elem", "nlocalele");

  nLocBas = h5r -> read_intScalar("Global_Mesh_Info", "nLocBas");

  int num_row, num_col;
  LIEN = h5r -> read_intMatrix("LIEN", "LIEN", num_row, num_col);
  
  SYS_T::print_fatal_if( num_row != nlocalele, "Error: ALocal_IEN::LIEN size does not match the number of element. \n");

  SYS_T::print_fatal_if( num_col != nLocBas, "Error: ALocal_IEN::LIEN size does not match the value of nLocBas. \n");
}

void ALocal_IEN::print_info() const
{
  std::cout<<"ALocal_IEN: \n";

  std::cout<<"nlocalele = "<<nlocalele<<'\n';
  std::cout<<"nLocBas = "<<nLocBas<<'\n';

  for(int ee = 0; ee<nlocalele; ++ee)
  {
    const std::vector<int> temp = get_LIEN( ee );
    std::cout<<ee<<" : \t";
    VEC_T::print(temp);
    std::cout<<'\n';
  }
}

// EOF
