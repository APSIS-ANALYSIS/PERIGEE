#include "ElemBC_3D_tet.hpp"

ElemBC_3D_tet::ElemBC_3D_tet( const std::vector<std::string> &vtkfileList,
   const int &elemtype )
{
  num_ebc = static_cast<int>( vtkfileList.size() );

  num_node     = new int [num_ebc];
  num_cell     = new int [num_ebc];
  cell_nLocBas = new int [num_ebc];

  pt_xyz.resize(num_ebc);
  tri_ien.resize(num_ebc);
  global_node.resize(num_ebc);
  global_cell.resize(num_ebc);

  std::cout<<"===> ElemBC_3D_tet specified by \n";

  for(int ii=0; ii<num_ebc; ++ii)
  {
    std::cout<<"     ebc_id = "<<ii<<": "<<vtkfileList[ii]<<'\n';

    if(elemtype == 501)
    {
      cell_nLocBas[ii] = 3; // linear triangle
      TET_T::read_vtp_grid( vtkfileList[ii], num_node[ii], num_cell[ii],
          pt_xyz[ii], tri_ien[ii], global_node[ii], global_cell[ii] );
    }
    else if(elemtype == 502)
    {
      cell_nLocBas[ii] = 6; // quadratic triangle

      TET_T::read_vtu_grid( vtkfileList[ii], num_node[ii], num_cell[ii],
          pt_xyz[ii], tri_ien[ii], global_node[ii], global_cell[ii] );
    }
  }

  std::cout<<"     is generated. \n";
}


ElemBC_3D_tet::~ElemBC_3D_tet()
{
  delete [] num_node; num_node = nullptr;
  delete [] num_cell; num_cell = nullptr;
  delete [] cell_nLocBas; cell_nLocBas = nullptr;
}


void ElemBC_3D_tet::print_info() const
{
  std::cout<<"========================= \n";
  std::cout<<"ElemBC_3D_tet : ";
  std::cout<<" num_ebc = "<<num_ebc<<std::endl;
  for(int ii=0; ii<num_ebc; ++ii)
  {
    std::cout<<"== ebc_id = "<<ii<<'\n';
    std::cout<<" num_node = "<<num_node[ii]<<'\t';
    std::cout<<" num_cell = "<<num_cell[ii]<<'\t';
    std::cout<<" cell_nLocBas = "<<cell_nLocBas[ii]<<'\n';
    std::cout<<" global node : ";
    VEC_T::print(global_node[ii]);
    std::cout<<" global elem : ";
    VEC_T::print(global_cell[ii]);
    std::cout<<" surface IEN : ";
    VEC_T::print(tri_ien[ii]);
  }
  std::cout<<"========================= \n";
}


// EOF
