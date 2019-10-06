#include "ALocal_IEN_wBubble.hpp"

ALocal_IEN_wBubble::ALocal_IEN_wBubble( 
    const std::string &fileBaseName, const int &cpu_rank,
    const int &nbubble_per_cell )
: ALocal_IEN( fileBaseName, cpu_rank )
{
  // ----------------------------------------------------------------
  // Read nlocghonode from part.h5 file
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const int nlgnode = h5r->read_intScalar("Local_Node", "nlocghonode");

  delete h5r; H5Fclose( file_id );
  // ----------------------------------------------------------------

  // the un-enriched mesh nLocBas defines the geometry
  const int nLocBas_ori = nLocBas;

  // Define physical nLocBas, here we add local bubble number
  nLocBas = nLocBas_ori + nbubble_per_cell;

  std::vector<int> temp = LIEN;

  // Define LIEN
  LIEN.clear(); LIEN.resize( nLocBas * nlocalele );

  for( int ee=0; ee<nlocalele; ++ee )
  {
    for( int ii=0; ii<nLocBas_ori; ++ii)
      LIEN[ee*nLocBas + ii] = temp[ee*nLocBas_ori + ii];

    const int offset = nlgnode + ee * nbubble_per_cell;

    for( int ii=0; ii<nbubble_per_cell; ++ii )
      LIEN[ee*nLocBas + nLocBas_ori + ii] = offset + ii; 
  }
}


ALocal_IEN_wBubble::~ALocal_IEN_wBubble()
{}


void ALocal_IEN_wBubble::print_info() const
{
  std::cout<<"ALocal_IEN_wBubble: \n";

  std::cout<<"nlocalele = "<<nlocalele<<'\n'; 
  std::cout<<"nLocBas = "<<nLocBas<<'\n';

  std::vector<int> temp; temp.resize(nLocBas);
  for(int ee=0; ee < nlocalele; ++ee)
  {
    get_LIEN_e(ee, &temp[0]);
    std::cout<<ee<<" : \t";
    VEC_T::print(temp);
    std::cout<<'\n';
  }
}

// EOF
