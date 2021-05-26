#include "ALocal_IEN.hpp"

ALocal_IEN::ALocal_IEN( const std::string &fileBaseName, const int &cpu_rank )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  nlocalele = h5r -> read_intScalar("Local_Elem", "nlocalele");

  nLocBas = h5r -> read_intScalar("Global_Mesh_Info", "nLocBas");

  h5r -> read_intMatrix("LIEN", "LIEN", LIEN, nlocalele, nLocBas);
  
  delete h5r; H5Fclose( file_id );

  SYS_T::print_fatal_if( int(LIEN.size()) != nLocBas * nlocalele, 
    "Error: ALocal_IEN::LIEN is in wrong format.\n" );
}

ALocal_IEN::~ALocal_IEN()
{}

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
