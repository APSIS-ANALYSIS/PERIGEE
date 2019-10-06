#include "BoundaryCond_Test.hpp"

void TEST_T::BC_Periodic_coordinate_check( const std::vector<double> &cPts,
    const BoundaryCond * const &bc, const double &tol )
{
  const int num_per_nodes = bc->get_num_per_nodes();
  
  int master, slave;
  
  double mx, my, mz, mw;
  double sx, sy, sz, sw;
  
  for(int ii=0; ii<num_per_nodes; ++ii)
  {
    master = bc->get_per_master_nodes(ii);
    slave  = bc->get_per_slave_nodes(ii);

    mx = cPts[master*4 + 0];
    my = cPts[master*4 + 1];
    mz = cPts[master*4 + 2];
    mw = cPts[master*4 + 3];

    sx = cPts[slave*4  + 0];
    sy = cPts[slave*4  + 1];
    sz = cPts[slave*4  + 2];
    sw = cPts[slave*4  + 3];

    if(abs(mx - sx) > tol)
    {
      std::cerr<<"Warning: master "<<master<<" slave "<<slave<<" x coordinate ";
      std::cerr<<mx<<'\t'<<sx<<'\t'<<mx - sx<<'\n';
    }

    if(abs(my - sy) > tol)
    {
      std::cerr<<"Warning: master "<<master<<" slave "<<slave<<" y coordinate ";
      std::cerr<<my<<'\t'<<sy<<'\t'<<my - sy<<'\n';
    }

    if(abs(mz - sz) > tol)
    {
      std::cerr<<"Warning: master "<<master<<" slave "<<slave<<" z coordinate ";
      std::cerr<<mz<<'\t'<<sz<<'\t'<<mz - sz<<'\n';
    }
    
    // Not check weights. We only require physical coordinates to be identical
    // for master-slave pairs.
    if(abs(mw - sw) > 1.0e-13)
    {
      std::cerr<<"Error: master "<<master<<" slave "<<slave<<" w coordinate ";
      std::cerr<<mw<<'\t'<<sw<<'\n';
    }
  }

  std::cout<<"BC_Periodic_coordinate_check passed with tolerance = "<<tol<<std::endl;
}


void TEST_T::correct_masterslave_ctrlPts( std::vector<double> &cPts,
          const BoundaryCond * const &bc )
{
  const int num_per_nodes = bc->get_num_per_nodes();
  
  const double tol = 1.0e-16;

  int master, slave;
  
  double mx, my, mz;// mw;
  double sx, sy, sz;// sw;
  
  int counter = 0;

  for(int ii=0; ii<num_per_nodes; ++ii)
  {
    master = bc->get_per_master_nodes(ii);
    slave  = bc->get_per_slave_nodes(ii);

    mx = cPts[master*4 + 0];
    my = cPts[master*4 + 1];
    mz = cPts[master*4 + 2];
    //mw = cPts[master*4 + 3];

    sx = cPts[slave*4  + 0];
    sy = cPts[slave*4  + 1];
    sz = cPts[slave*4  + 2];
    //sw = cPts[slave*4  + 3];

    if(abs(mx - sx) > tol)
    {
      cPts[slave*4 + 0] = mx;
      counter ++;
    }

    if(abs(my - sy) > tol)
    {
      cPts[slave*4 + 1] = my;
      counter ++;
    }

    if(abs(mz - sz) > tol)
    {
      cPts[slave*4 + 2] = mz;
      counter ++;
    }
  }

  std::cout<<"With tolerance = "<<tol<<" totally "<<counter
    <<" coordinates have been corrected."<<std::endl;
}




void TEST_T::NBC_compare( BoundaryCond2D const * const &bc,
         INodalBC const * const &nbc )
{
  assert(bc->get_num_dir_nodes() == nbc->get_num_dir_nodes());
  assert(bc->get_num_per_nodes() == nbc->get_num_per_nodes());

  for(unsigned int ii=0; ii<nbc->get_num_dir_nodes(); ++ii)
    assert(bc->get_dir_nodes(ii) == int(nbc->get_dir_nodes(ii)));

  for(unsigned int ii=0; ii<nbc->get_num_per_nodes(); ++ii)
  {
    assert(bc->get_per_slave_nodes(ii) == int(nbc->get_per_slave_nodes(ii)));
    assert(bc->get_per_master_nodes(ii) == int(nbc->get_per_master_nodes(ii)));
  }

}


void TEST_T::EBC_node_compatibility_check( const ElemBC * const &bc,
          std::vector<double> &ctrlPts )
{
  const int nebc = bc->get_num_ebc();
  for(int ii=0; ii<nebc; ++ii)
  {
    const int nnode = bc->get_num_node(ii);
    for(int jj=0; jj<nnode; ++jj)
    {
      int nn = bc->get_global_node(ii,jj);
      double nx = bc->get_pt_xyz(ii, jj, 0);
      double ny = bc->get_pt_xyz(ii, jj, 1);
      double nz = bc->get_pt_xyz(ii, jj, 2);

      double vx = ctrlPts[nn*3];
      double vy = ctrlPts[nn*3+1];
      double vz = ctrlPts[nn*3+2];

      if(nx != vx)
      {
        std::cout<<"Error: ElemBC "<<ii<<" node "<<jj<<" x-coor does not match volume point. \n";
        std::cout<<nx<<'\t'<<vx<<std::endl;
        exit(EXIT_FAILURE);
      }
      
      if(ny != vy)
      {
        std::cout<<"Error: ElemBC "<<ii<<" node "<<jj<<" y-coor does not match volume point. \n";
        std::cout<<ny<<'\t'<<vy<<std::endl;
        exit(EXIT_FAILURE);
      }
      
      if(nz != vz)
      {
        std::cout<<"Error: ElemBC "<<ii<<" node "<<jj<<" z-coor does not match volume point. \n";
        std::cout<<nz<<'\t'<<vz<<std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  std::cout<<"EBC_node_compatibility_check passed! \n";
}


void TEST_T::EBC_cell_IEN_check( const ElemBC * const &bc,
    const IIEN * const &ien )
{
  const int nebc = bc->get_num_ebc();
  std::vector<int> tet_ien;
  std::vector<int>::iterator it;
  for(int ii=0; ii<nebc; ++ii)
  {
    const int ncell = bc->get_num_cell(ii);
    
    for(int jj=0; jj<ncell; ++jj)
    {
      tet_ien.clear();
      int eindex = bc->get_global_cell(ii, jj);
      tet_ien.push_back(ien->get_IEN(eindex, 0));
      tet_ien.push_back(ien->get_IEN(eindex, 1));
      tet_ien.push_back(ien->get_IEN(eindex, 2));
      tet_ien.push_back(ien->get_IEN(eindex, 3));

      int c0 = bc->get_ien(ii, jj, 0);
      int c1 = bc->get_ien(ii, jj, 1);
      int c2 = bc->get_ien(ii, jj, 2);

      c0 = bc->get_global_node(ii, c0);
      c1 = bc->get_global_node(ii, c1);
      c2 = bc->get_global_node(ii, c2);

      it = find(tet_ien.begin(), tet_ien.end(), c0);
      if(it == tet_ien.end())
      {
        std::cout<<"Error: ElemBC "<<ii<<" cell "<<jj<<" failed to find node 0. \n";
        exit(EXIT_FAILURE);
      }
  
      it = find(tet_ien.begin(), tet_ien.end(), c1);
      if(it == tet_ien.end())
      {
        std::cout<<"Error: ElemBC "<<ii<<" cell "<<jj<<" failed to find node 1. \n";
        exit(EXIT_FAILURE);
      }
  
      it = find(tet_ien.begin(), tet_ien.end(), c2);
      if(it == tet_ien.end())
      {
        std::cout<<"Error: ElemBC "<<ii<<" cell "<<jj<<" failed to find node 2. \n";
        exit(EXIT_FAILURE);
      }
    }
  }

  std::cout<<"EBC_cell_IEN_check passed! \n";
}

// EOF
