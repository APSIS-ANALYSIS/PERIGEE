#include "Tissue_property.hpp"

Tissue_property::Tissue_property( const std::string &fileBaseName, const int &rank )
:cpu_rank( rank )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  nlocghonode_s = h5r -> read_intScalar("directionBasis_loc", "nlocghonode_s");
  if(nlocghonode_s>0)
  {
    node_locgho_solid = h5r -> read_intVector("directionBasis_loc", "node_locgho_solid");
    basis_r = h5r -> read_Vector_3_Vector("directionBasis_loc", "loc_basis_r");
    basis_c = h5r -> read_Vector_3_Vector("directionBasis_loc", "loc_basis_c");
    basis_l = h5r -> read_Vector_3_Vector("directionBasis_loc", "loc_basis_l");
  }

  delete h5r; H5Fclose( file_id );
}

Vector_3 Tissue_property::get_basis_r( const int &nn ) const
{
  const int index = node_locgho_solid[nn];
  ASSERT(index > -1, "ERROR: Tissue_property::get_basis_r the input node is a fluid node.");
  return basis_r[index];
}

Vector_3 Tissue_property::get_basis_c( const int &nn ) const
{
  const int index = node_locgho_solid[nn];
  ASSERT(index > -1, "ERROR: Tissue_property::get_basis_c the input node is a fluid node.");
  return basis_c[index];
}

Vector_3 Tissue_property::get_basis_l( const int &nn ) const
{
  const int index = node_locgho_solid[nn];
  ASSERT(index > -1, "ERROR: Tissue_property::get_basis_l the input node is a fluid node.");
  return basis_l[index];
}

void Tissue_property::print_info() const
{
  std::cout<<"\n -- cpu "<<cpu_rank<<" node info: \n";
  std::cout<<"\n nlocghonode_s: "<<nlocghonode_s<<std::endl;
  std::cout<<"\n node_locgho_solid: \n";
  VEC_T::print(node_locgho_solid);
  const int size = VEC_T::get_size( basis_r );
  std::cout<<"\n basis_r: \n";
  for(int ii=0; ii<size; ii++) basis_r[ii].print(); 
  std::cout<<"\n basis_c: \n";
  for(int ii=0; ii<size; ii++) basis_c[ii].print();
  std::cout<<"\n basis_l: \n";
  for(int ii=0; ii<size; ii++) basis_l[ii].print();
}

// EOF
