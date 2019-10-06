#include "ElemBC_3D_tet4_outflow.hpp"

ElemBC_3D_tet4_outflow::ElemBC_3D_tet4_outflow( 
    const std::vector<std::string> &vtpfileList,
    const std::vector< std::vector<double> > &outlet_normal_vec )
: ElemBC_3D_tet4( vtpfileList )
{
  // Copy the outward normal vector in outNormal
  outNormal.resize( num_ebc );
  for(int ii=0; ii<num_ebc; ++ii) outNormal[ii]= outlet_normal_vec[ii];

  // Allocate the intNA array
  intNA.resize( num_ebc );
  for( int ii=0; ii<num_ebc; ++ii ) intNA[ii].resize( num_node[ii] );

  // Zero all the entries in intNA
  for( int ii=0; ii<num_ebc; ++ii )
    for(int jj=0; jj<num_node[ii]; ++jj )
      intNA[ii][jj] = 0.0;

  // Setup quadrature rule for surface triangles
  const int nqp_tri = 3;

  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );
  FEAElement * elems = new FEAElement_Triangle3_3D_der0( nqp_tri );

  double ectrl_x[3]; double ectrl_y[3]; double ectrl_z[3];
  int node_idx[3]; double R[3];
  double nx, ny, nz, surface_area;

  for( int fid=0; fid < num_ebc; ++fid )
  {
    for( int ee=0; ee<num_cell[fid]; ++ee )
    {
      node_idx[0] = get_ien(fid, ee, 0);
      node_idx[1] = get_ien(fid, ee, 1);
      node_idx[2] = get_ien(fid, ee, 2);

      ectrl_x[0] = get_pt_xyz(fid, node_idx[0], 0);
      ectrl_x[1] = get_pt_xyz(fid, node_idx[1], 0);
      ectrl_x[2] = get_pt_xyz(fid, node_idx[2], 0);

      ectrl_y[0] = get_pt_xyz(fid, node_idx[0], 1);
      ectrl_y[1] = get_pt_xyz(fid, node_idx[1], 1);
      ectrl_y[2] = get_pt_xyz(fid, node_idx[2], 1);

      ectrl_z[0] = get_pt_xyz(fid, node_idx[0], 2);
      ectrl_z[1] = get_pt_xyz(fid, node_idx[1], 2);
      ectrl_z[2] = get_pt_xyz(fid, node_idx[2], 2);

      elems -> buildBasis(quads, ectrl_x, ectrl_y, ectrl_z);

      for(int qua=0; qua<nqp_tri; ++qua)
      {
        elems -> get_R(qua, R);
        elems -> get_2d_normal_out(qua, nx, ny, nz, surface_area);
        
        const double gwts = surface_area * quads -> get_qw( qua );

        intNA[fid][node_idx[0]] += gwts * R[0];
        intNA[fid][node_idx[1]] += gwts * R[1];
        intNA[fid][node_idx[2]] += gwts * R[2];
      }
    }
  }
  // Clean the quadrature rule and the triangle element objects
  delete quads; delete elems;

  std::cout<<"===> ElemBC_3D_tet4_outflow intNA vector generated.\n";
}


ElemBC_3D_tet4_outflow::~ElemBC_3D_tet4_outflow()
{
  for(int ii=0; ii<num_ebc; ++ii) VEC_T::clean( intNA[ii] );

  VEC_T::clean( intNA );
}


void ElemBC_3D_tet4_outflow::print_info() const
{
  ElemBC_3D_tet4::print_info();

  for(int face=0; face<num_ebc; ++face)
  {
    std::cout<<"Surface id : "<<face<<std::endl;
    VEC_T::print( intNA[face] );
  }

}

// EOF
