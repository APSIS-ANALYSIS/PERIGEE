#include "ElemBC_2D.hpp"

ElemBC_2D::ElemBC_2D( const std::vector<int> &edge_idx,
    const std::vector<int> &in_num_cell,
    const std::vector<int> &in_num_node,
    const std::vector<int> &in_nlocbas,
    const std::vector< std::vector<double> > &in_pt_coor,
    const std::vector< std::vector<int> > &in_ien_loc,
    const std::vector< std::vector<int> > &in_bcpt,
    const std::vector< std::vector<int> > &edge2elem,
    const IIEN * const &VIEN,
    const std::vector<double> &vctrlPts )
{
  num_ebc = static_cast<int>( edge_idx.size() );

  num_node = new int [num_ebc];
  num_cell = new int [num_ebc];
  cell_nLocBas = new int [num_ebc];

  pt_xyz.resize(num_ebc);
  cell_interior_pt.resize(num_ebc);
  tri_ien.resize(num_ebc);
  global_node.resize(num_ebc);
  global_cell.resize(num_ebc);

  const int max_idx = static_cast<int>( in_num_cell.size() );

  std::cout<<"===> ElemBC_2D specified by \n";
  for(int ii=0; ii<num_ebc; ++ii)
  {
    const int idx = edge_idx[ii];
    if(idx >= max_idx || idx < 0) SYS_T::print_fatal("Error: edge index out of range. \n");
  
    std::cout<<"     ebc_id = "<<ii<<" : edge "<<idx<<'\n';

    num_node[ii] = in_num_node[idx];
    num_cell[ii] = in_num_cell[idx];
    cell_nLocBas[ii] = in_nlocbas[idx];
  
    pt_xyz[ii] = in_pt_coor[idx];
    tri_ien[ii] = in_ien_loc[idx];
    global_node[ii] = in_bcpt[idx];
    global_cell[ii] = edge2elem[idx];
 
    // The idea of the algorithm is this. In Gmsh, the output of line, triangle,
    // and tet gives the corner points first, no matter what the degree is.
    // To identify the line in the triangle, we only need the first two points
    // for the line segment and the first three points in the triangle, and
    // the point not in the line segment is the interior point. It will help
    // define the outward normal vector for this line segment. Similarily, we
    // only need the first three point from a triangle and the first four
    // point in a tet, and the point not belonging to the triangle gives the 
    // volumetric interior point, it may help define the outward normal of the
    // surface triangle element.
    int linpt [2] = {0, 0};
    int tript [3] = {0, 0, 0}; 
    bool is0, is1, is2;
    int loc_int_idx; // interior point local(elemental) index
    int glo_int_idx; // interior point global index
    for(int jj=0; jj<num_cell[ii]; ++jj)
    {
      int loc_line_0 = tri_ien[ii][ cell_nLocBas[ii]*jj + 0 ];
      int loc_line_1 = tri_ien[ii][ cell_nLocBas[ii]*jj + 1 ];
      linpt[0] = global_node[ii][loc_line_0];
      linpt[1] = global_node[ii][loc_line_1];

      tript[0] = VIEN->get_IEN( global_cell[ii][jj], 0 );
      tript[1] = VIEN->get_IEN( global_cell[ii][jj], 1 );
      tript[2] = VIEN->get_IEN( global_cell[ii][jj], 2 );
    
      is0 = ( (tript[0] == linpt[0]) || (tript[0] == linpt[1]) );
      is1 = ( (tript[1] == linpt[0]) || (tript[1] == linpt[1]) );
      is2 = ( (tript[2] == linpt[0]) || (tript[2] == linpt[1]) );

      if(!is0 && is1 && is2) loc_int_idx = 0;
      else if(is0 && !is1 && is2) loc_int_idx = 1;
      else if(is0 && is1 && !is2) loc_int_idx = 2;
      else SYS_T::print_fatal("Error: edge is incompatible with triangle. \n");
    
      glo_int_idx = VIEN->get_IEN( global_cell[ii][jj], loc_int_idx );
      
      cell_interior_pt[ii].push_back( vctrlPts[3*glo_int_idx + 0] ); 
      cell_interior_pt[ii].push_back( vctrlPts[3*glo_int_idx + 1] ); 
      cell_interior_pt[ii].push_back( vctrlPts[3*glo_int_idx + 2] ); 
    }
  }
  std::cout<<"     is generated. \n";
}


ElemBC_2D::~ElemBC_2D()
{
  delete [] num_node; num_node = NULL;
  delete [] num_cell; num_cell = NULL;
  delete [] cell_nLocBas; cell_nLocBas = NULL;
}


void ElemBC_2D::print_info() const
{
  std::cout<<"========================= \n";
  std::cout<<"ElemBC_2D : ";
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
    std::cout<<" interior pt coordinates : ";
    VEC_T::print(cell_interior_pt[ii]);
  }
  std::cout<<"========================= \n";
}


// EOF
