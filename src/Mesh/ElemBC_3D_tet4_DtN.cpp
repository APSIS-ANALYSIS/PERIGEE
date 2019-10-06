#include "ElemBC_3D_tet4_DtN.hpp"

ElemBC_3D_tet4_DtN::ElemBC_3D_tet4_DtN( 
    const std::vector<std::string> &vtpfileList,
    const std::vector<int> &dtn_bc_list )
: ElemBC_3D_tet4( vtpfileList )
{
  num_dtn_surface = static_cast<int>( dtn_bc_list.size() );

  // Check the format of dtn_bc_list by basic logic
  SYS_T::print_fatal_if(num_dtn_surface > num_ebc,
      "Error: the number of DtN surface is larger than the number of ebc.\n" ); 

  for(int ii=0; ii<num_dtn_surface; ++ii)
    SYS_T::print_fatal_if( dtn_bc_list[ii] >= num_ebc,
      "Error: the dtn surface id is out of range. \n" );
   
  // Print the DtN surface on screen
  std::cout<<"     DtN surfaces: \n";
  for(int ii=0; ii<num_dtn_surface; ++ii)
    std::cout<<"     "<<dtn_bc_list[ii]<<'\t'<<vtpfileList[ dtn_bc_list[ii] ]<<'\n';
  
  // Allocate the intNA array
  intNA.resize( num_dtn_surface );
  for(int ii=0; ii<num_dtn_surface; ++ii) 
    intNA[ii].resize( 3*num_node[dtn_bc_list[ii]] );
  
  for( int ii=0; ii<num_dtn_surface; ++ii )
    for( int jj=0; jj<3*num_node[dtn_bc_list[ii]]; ++jj )
      intNA[ii][jj] = 0.0;

  // Quadrature rule for surface integration
  const int nqp_tri = 3;

  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );
  FEAElement * elems = new FEAElement_Triangle3_3D_der0( nqp_tri );
 
  double ectrl_x[3]; double ectrl_y[3]; double ectrl_z[3];
  int node_idx[3]; double R[3];
  double nx, ny, nz, surface_area;
  
  // elem_normal vector stores the normal vector and calculate their statistics
  std::vector<double> elem_nx, elem_ny, elem_nz;
  
  // Loop over DtN surface kk 
  for( int kk=0; kk<num_dtn_surface; ++kk )
  {
    const int face_id = dtn_bc_list[kk];

    elem_nx.clear(); elem_ny.clear(); elem_nz.clear();

    for( int ee=0; ee<num_cell[face_id]; ++ee )
    {
      node_idx[0] = get_ien(face_id, ee, 0);
      node_idx[1] = get_ien(face_id, ee, 1);
      node_idx[2] = get_ien(face_id, ee, 2);
      ectrl_x[0] = get_pt_xyz(face_id, node_idx[0], 0);
      ectrl_x[1] = get_pt_xyz(face_id, node_idx[1], 0);
      ectrl_x[2] = get_pt_xyz(face_id, node_idx[2], 0);
      ectrl_y[0] = get_pt_xyz(face_id, node_idx[0], 1);
      ectrl_y[1] = get_pt_xyz(face_id, node_idx[1], 1);
      ectrl_y[2] = get_pt_xyz(face_id, node_idx[2], 1);
      ectrl_z[0] = get_pt_xyz(face_id, node_idx[0], 2);
      ectrl_z[1] = get_pt_xyz(face_id, node_idx[1], 2);
      ectrl_z[2] = get_pt_xyz(face_id, node_idx[2], 2);

      elems -> buildBasis(quads, ectrl_x, ectrl_y, ectrl_z);
      
      for(int qua=0; qua<nqp_tri; ++qua)
      {
        elems -> get_R(qua, R);
        elems -> get_2d_normal_out(qua, nx, ny, nz, surface_area);
        
        elem_nx.push_back(nx); elem_ny.push_back(ny); elem_nz.push_back(nz);
         
        double gwts = surface_area * quads->get_qw(qua);

        intNA[kk][node_idx[0]*3+0] += gwts * R[0] * nx;
        intNA[kk][node_idx[0]*3+1] += gwts * R[0] * ny;
        intNA[kk][node_idx[0]*3+2] += gwts * R[0] * nz;

        intNA[kk][node_idx[1]*3+0] += gwts * R[1] * nx;
        intNA[kk][node_idx[1]*3+1] += gwts * R[1] * ny;
        intNA[kk][node_idx[1]*3+2] += gwts * R[1] * nz;

        intNA[kk][node_idx[2]*3+0] += gwts * R[2] * nx;
        intNA[kk][node_idx[2]*3+1] += gwts * R[2] * ny;
        intNA[kk][node_idx[2]*3+2] += gwts * R[2] * nz;
      }
    }
    std::cout<<vtpfileList[face_id]<<" normal vector: ";
    std::cout<<MATH_T::get_mean(elem_nx)<<'('<<MATH_T::get_std_dev(elem_nx)<<") ";
    std::cout<<MATH_T::get_mean(elem_ny)<<'('<<MATH_T::get_std_dev(elem_ny)<<") ";
    std::cout<<MATH_T::get_mean(elem_nz)<<'('<<MATH_T::get_std_dev(elem_nz)<<")\n";
  }
  
  delete quads; delete elems; // intNA has been assembled
  
  // write intNA into a h5 file
  hid_t cpbc_file = H5Fcreate("cpBC.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cpbch5w = new HDF5_Writer(cpbc_file);

  cpbch5w->write_intScalar("num_dtn_surface", num_dtn_surface);

  cpbch5w->write_intVector("dtn_bc_list", dtn_bc_list);

  const std::string intna_base("intNA_");
  const std::string nodem_base("node_mapping_");

  for(int kk=0; kk<num_dtn_surface; ++kk)
  {
    std::string intna_name( intna_base );
    intna_name.append(SYS_T::to_string(kk));
    cpbch5w->write_doubleVector( intna_name.c_str(), intNA[kk] );    
  
    std::string nodem( nodem_base );
    nodem.append(SYS_T::to_string(kk));
    cpbch5w->write_intVector( nodem.c_str(), global_node[ dtn_bc_list[kk] ] );
  }

  delete cpbch5w; H5Fclose(cpbc_file);
  // finish writing the coupled BC file
}


ElemBC_3D_tet4_DtN::~ElemBC_3D_tet4_DtN()
{
  for(int ii=0; ii<num_dtn_surface; ++ii) VEC_T::clean(intNA[ii]);

  VEC_T::clean(intNA);
}


void ElemBC_3D_tet4_DtN::print_info() const
{
  ElemBC_3D_tet4::print_info();

  std::cout<<"DtN info: \n";
  std::cout<<" num_dtn_surface = "<<num_dtn_surface<<std::endl;
}

// EOF
