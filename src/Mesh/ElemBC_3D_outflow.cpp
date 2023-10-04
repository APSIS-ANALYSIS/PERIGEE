#include "ElemBC_3D_outflow.hpp"

ElemBC_3D_outflow::ElemBC_3D_outflow(
    const std::vector<std::string> &vtkfileList,
    const std::vector< Vector_3 > &outlet_normal_vec,
    const int &elemtype )
: ElemBC_3D( vtkfileList, elemtype )
{
  SYS_T::print_fatal_if(outlet_normal_vec.size() != static_cast<unsigned int>( num_ebc ),
      "Error: ElemBC_3D_outflow constructor: the input normal vector length does not match the number of outlets.\n");

  outNormal.resize( num_ebc );
  for(int ii=0; ii<num_ebc; ++ii) outNormal[ii]= outlet_normal_vec[ii];

  intNA.resize( num_ebc );
  for( int ii=0; ii<num_ebc; ++ii ) intNA[ii].resize( num_node[ii] );

  for( int ii=0; ii<num_ebc; ++ii )
    for(int jj=0; jj<num_node[ii]; ++jj )
      intNA[ii][jj] = 0.0;

  if(elem_type == 501)
  {
    const int nqp_tri = 3;

    IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );
    FEAElement * elems = new FEAElement_Triangle3_3D_der0( nqp_tri );

    double ectrl_x[3]{}; double ectrl_y[3]{}; double ectrl_z[3]{};
    int node_idx[3]{}; double R[3]{};

    for( int fid=0; fid<num_ebc; ++fid )
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

          const double gwts = elems -> get_detJac( qua ) * quads -> get_qw( qua );

          intNA[fid][node_idx[0]] += gwts * R[0];
          intNA[fid][node_idx[1]] += gwts * R[1];
          intNA[fid][node_idx[2]] += gwts * R[2];
        }
      }
    }
    delete quads; delete elems;
  }
  else if(elem_type == 502)
  {
    const int nqp_tri = 6;

    IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );
    FEAElement * elems = new FEAElement_Triangle6_3D_der0( nqp_tri );
  
    double ectrl_x[6]{}; double ectrl_y[6]{}; double ectrl_z[6]{};
    int node_idx[6]{}; double R[6]{};
    
    for( int fid=0; fid<num_ebc; ++fid )
    {
      for( int ee=0; ee<num_cell[fid]; ++ee )
      {
        for(int ii=0; ii<6; ++ii)
        {
          node_idx[ii] = get_ien(fid, ee, ii);
          ectrl_x[ii] = get_pt_xyz(fid, node_idx[ii], 0);
          ectrl_y[ii] = get_pt_xyz(fid, node_idx[ii], 1);
          ectrl_z[ii] = get_pt_xyz(fid, node_idx[ii], 2);
        }

        elems -> buildBasis(quads, ectrl_x, ectrl_y, ectrl_z);

        for(int qua=0; qua<nqp_tri; ++qua)
        {
          elems -> get_R(qua, R);

          const double gwts = elems -> get_detJac(qua) * quads -> get_qw( qua );

          for(int ii=0; ii<6; ++ii) intNA[fid][node_idx[ii]] += gwts * R[ii];
        
        } // loop over quadrature points
      } // loop over elements
    } // loop over outlet surfaces
    delete quads; delete elems;
  }
  else if(elem_type == 601)
  {
    const int nqp_quad = 2;  // nqp_quad represents the number of quadrature points in one direction 

    IQuadPts * quads = new QuadPts_Gauss_Quad( nqp_quad );
    FEAElement * elems = new FEAElement_Quad4_3D_der0( nqp_quad * nqp_quad );

    double ectrl_x[4]{}; double ectrl_y[4]{}; double ectrl_z[4]{};
    int node_idx[4]{}; double R[4]{};

    for( int fid=0; fid<num_ebc; ++fid )
    {
      for( int ee=0; ee<num_cell[fid]; ++ee )
      {
        node_idx[0] = get_ien(fid, ee, 0);
        node_idx[1] = get_ien(fid, ee, 1);
        node_idx[2] = get_ien(fid, ee, 2);
        node_idx[3] = get_ien(fid, ee, 3);

        ectrl_x[0] = get_pt_xyz(fid, node_idx[0], 0);
        ectrl_x[1] = get_pt_xyz(fid, node_idx[1], 0);
        ectrl_x[2] = get_pt_xyz(fid, node_idx[2], 0);
        ectrl_x[3] = get_pt_xyz(fid, node_idx[3], 0);

        ectrl_y[0] = get_pt_xyz(fid, node_idx[0], 1);
        ectrl_y[1] = get_pt_xyz(fid, node_idx[1], 1);
        ectrl_y[2] = get_pt_xyz(fid, node_idx[2], 1);
        ectrl_y[3] = get_pt_xyz(fid, node_idx[3], 1);

        ectrl_z[0] = get_pt_xyz(fid, node_idx[0], 2);
        ectrl_z[1] = get_pt_xyz(fid, node_idx[1], 2);
        ectrl_z[2] = get_pt_xyz(fid, node_idx[2], 2);
        ectrl_z[3] = get_pt_xyz(fid, node_idx[3], 2);

        elems -> buildBasis(quads, ectrl_x, ectrl_y, ectrl_z);

        for(int qua=0; qua<nqp_quad * nqp_quad; ++qua)
        {
          elems -> get_R(qua, R);

          const double gwts = elems -> get_detJac( qua ) * quads -> get_qw( qua );

          intNA[fid][node_idx[0]] += gwts * R[0];
          intNA[fid][node_idx[1]] += gwts * R[1];
          intNA[fid][node_idx[2]] += gwts * R[2];
          intNA[fid][node_idx[3]] += gwts * R[3];
        }
      }
    }
    delete quads; delete elems;
  }
  else if(elem_type == 602)
  {
    const int nqp_quad = 4;

    IQuadPts * quads = new QuadPts_Gauss_Quad( nqp_quad );
    FEAElement * elems = new FEAElement_Quad9_3D_der0( nqp_quad * nqp_quad );
  
    double ectrl_x[9]{}; double ectrl_y[9]{}; double ectrl_z[9]{};
    int node_idx[9]{}; double R[9]{};  
    
    for( int fid=0; fid<num_ebc; ++fid )
    {
      for( int ee=0; ee<num_cell[fid]; ++ee )
      {
        for(int ii=0; ii<9; ++ii)
        {
          node_idx[ii] = get_ien(fid, ee, ii);
          ectrl_x[ii] = get_pt_xyz(fid, node_idx[ii], 0);
          ectrl_y[ii] = get_pt_xyz(fid, node_idx[ii], 1);
          ectrl_z[ii] = get_pt_xyz(fid, node_idx[ii], 2);
        }

        elems -> buildBasis(quads, ectrl_x, ectrl_y, ectrl_z);

        for(int qua=0; qua<nqp_quad * nqp_quad; ++qua)
        {
          elems -> get_R(qua, R);

          const double gwts = elems -> get_detJac(qua) * quads -> get_qw( qua );

          for(int ii=0; ii<9; ++ii) intNA[fid][node_idx[ii]] += gwts * R[ii];
        
        } // loop over quadrature points
      } // loop over elements
    } // loop over outlet surfaces
    delete quads; delete elems;
  }
  else SYS_T::print_exit("Error: ElemBC_3D_outflow constructor: unknown element type.\n");

  std::cout<<"     ElemBC_3D_outflow intNA vector generated.\n";
}


ElemBC_3D_outflow::~ElemBC_3D_outflow()
{
  for(int ii=0; ii<num_ebc; ++ii) VEC_T::clean( intNA[ii] );

  VEC_T::clean( intNA );
}


void ElemBC_3D_outflow::print_info() const
{
  ElemBC_3D::print_info();

  for(int face=0; face<num_ebc; ++face)
  {
    std::cout<<"Surface id : "<<face<<std::endl;
    VEC_T::print( intNA[face] );
  }
}

// EOF