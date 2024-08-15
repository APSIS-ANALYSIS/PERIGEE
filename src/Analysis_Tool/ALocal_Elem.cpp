#include "ALocal_Elem.hpp"

ALocal_Elem::ALocal_Elem(const std::string &fileBaseName, const int &cpu_rank)
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  elem_loc = h5r->read_intVector( "/Local_Elem", "elem_loc" );

  nlocalele = h5r->read_intScalar( "/Local_Elem", "nlocalele" );

  isTagged = h5r -> check_data("/Local_Elem/elem_phy_tag");

  if( isTagged )
  {
    elem_tag = h5r->read_intVector("/Local_Elem", "elem_phy_tag");
    SYS_T::print_fatal_if( VEC_T::get_size( elem_tag ) != nlocalele, "Error: ALocal_Elem::ALocal_Elem function elem_tag_length is not equal to nlocalele.\n");
  }
  else
    elem_tag.clear();

  // isRotated = h5r -> check_data("/Local_Elem/elem_rotated_tag");

  // if( isRotated )
  //   elem_rotated = h5r->read_intVector("/Local_Elem", "elem_rotated_tag");
  // else
  //   elem_rotated.clear();
    
  delete h5r; H5Fclose( file_id );
}

int ALocal_Elem::get_nlocalele( const int &tag_val ) const
{
  if( isTagged )
  {
    int counter = 0;
    for(int ee=0; ee<nlocalele; ++ee)
    {
      if( elem_tag[ee] == tag_val) counter += 1;
    }
    return counter;
  }
  else
  {
    SYS_T::print_fatal("Error: ALocal_Elem::get_nlocalele with input tag value should be called when isTagged = true. \n");
    return -1;
  }
}

// int ALocal_Elem::get_nlocalele_rotated( const int &rotated_val ) const
// {
//   if( isRotated )
//   {
//     int counter = 0;
//     for(int ee=0; ee<nlocalele; ++ee)
//     {
//       if( elem_rotated[ee] == rotated_val) counter += 1;
//     }
//     return counter;
//   }
//   else
//   {
//     SYS_T::print_fatal("Error: ALocal_Elem::get_nlocalele_rotated with input rotated tag value should be called when isRotated = true. \n");
//     return -1;
//   }
// }

void ALocal_Elem::print_info() const
{
  std::cout<<"nlocalelem: "<<nlocalele<<std::endl;
  if( isTagged )
  {
    std::cout<<"elem_loc(elem_tag): \n";
    for(int ii=0; ii<get_nlocalele(); ++ii)
      std::cout<<get_elem_loc(ii)<<"("<<get_elem_tag(ii)<<")"<<'\t';
    std::cout<<std::endl;
  }
  // if ( isRotated )
  // {
  //   std::cout<<"elem_loc(elem_rotated): \n";
  //   for(int ii=0; ii<get_nlocalele(); ++ii)
  //     std::cout<<get_elem_loc(ii)<<"("<<get_elem_rotated(ii)<<")"<<'\t';
  //   std::cout<<std::endl;
  // }
  else
  {
    std::cout<<"elem_loc: \n";
    VEC_T::print(elem_loc);
  }
}

// EOF
