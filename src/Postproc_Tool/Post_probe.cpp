#include "Post_probe.hpp"

Post_probe::Post_probe( const std::string &part_bname_in,
    const std::string &sol_bname_in,
    const int &in_sdeg, const int &in_tdeg, const int &in_udeg,
    const int &nlocbas_in, const int &in_nfunc, const int &in_dof )
: partfile_bname( part_bname_in ), sol_bname( sol_bname_in ),
  s_degree(in_sdeg), t_degree(in_tdeg), u_degree(in_udeg), 
  nLocBas(nlocbas_in), nFunc( in_nfunc ), dof( in_dof )
{
  // ----- read epart.h5->part list 
  hid_t file_id = H5Fopen("epart.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  
  hid_t drank;
  hsize_t * ddims;
  HDF5_Reader * h5reader = new HDF5_Reader( file_id ); 
  h5reader->read_intArray( "/", "part", drank, ddims, elepart );
  
  H5Fclose( file_id ); 
  // ----- finish reading h5 file

  // format check
  if(drank > 1)
  {
    std::cerr<<"Error: The epart.h5 dataset part is read incorrectly. \n";
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  
  nelement = ddims[0];

  delete [] ddims; delete h5reader;
}


Post_probe::~Post_probe()
{
  delete [] elepart; elepart = NULL;
}


void Post_probe::print_info() const
{
  std::cout<<"Post_probe object: \n";
  std::cout<<" Partition files' base name is "<<partfile_bname<<std::endl;
  std::cout<<" Polynomial degrees are "<<s_degree<<'\t'<<t_degree<<'\t'<<u_degree<<std::endl;
  std::cout<<" nLocBas = "<<nLocBas<<std::endl;
  std::cout<<" Solution file base name is "<<sol_bname<<std::endl;
  std::cout<<" there are "<<dof<<" degrees of freedom. \n";
  std::cout<<" There are "<<nelement<<" elements totally. \n";
  //cout<<" Their partition is: \n";
  //for(int ii=0; ii<nelement; ++ii)
  //  cout<<ii<<'\t'<<elepart[ii]<<'\n';
}



int Post_probe::get_elem_cpu( const int &ee ) const
{
  if( ee >= nelement || ee < 0 )
  {
    std::cerr<<"Error: get_elem_locindex, the input element index "<<ee<<" is beyond scope. \n";
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  // based on the elepart array, give the cpu index of the element with global
  // index ee 
  return elepart[ee]; 
}


void Post_probe::readPETSc_vec( const int &sol_index,
    const std::vector<int> &gloIEN,
    double * const &veccopy ) const
{
  // Generate solution name
  int aux = 900000000 + sol_index;
  std::ostringstream temp;
  temp << aux;
  std::string sol_name(sol_bname);
  sol_name.append(temp.str());

  // read in binary file
  Vec sol_temp;
  VecCreate(PETSC_COMM_SELF, &sol_temp);
  VecSetType(sol_temp, VECSEQ);

  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_SELF, sol_name.c_str(), FILE_MODE_READ, &viewer);
  VecLoad(sol_temp, viewer);
  PetscViewerDestroy(&viewer);

  // Check the solution size
  int sol_length;
  VecGetSize(sol_temp, &sol_length);
  if( sol_length != nFunc * dof )
  {
    std::cerr<<"Error: The length of the binary vector does not match the partition file. \n";
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  // copy the correct solution entries
  double * avec;
  VecGetArray(sol_temp, &avec);

  for(unsigned int ii=0; ii<gloIEN.size(); ++ii)
  {
    for(int jj=0; jj<dof; ++jj)
      veccopy[ii*dof+jj] = avec[gloIEN[ii] * dof + jj];
  }
  
  VecRestoreArray(sol_temp, &avec); 
  VecDestroy(&sol_temp);
}



void Post_probe::readPETSc_full( const int &sol_index, std::vector<double> &out ) const
{
  // Generate solution name
  int aux = 900000000 + sol_index;
  std::ostringstream temp;
  temp << aux;
  std::string sol_name(sol_bname);
  sol_name.append(temp.str());

  // read in binary file
  Vec sol_temp;
  VecCreate(PETSC_COMM_SELF, &sol_temp);
  VecSetType(sol_temp, VECSEQ);

  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_SELF, sol_name.c_str(), FILE_MODE_READ, &viewer);
  VecLoad(sol_temp, viewer);
  PetscViewerDestroy(&viewer);

  // Check the solution size
  int sol_length;
  VecGetSize(sol_temp, &sol_length);
  if( sol_length != nFunc * dof )
  {
    std::cerr<<"Error: The length of the binary vector does not match the partition file. \n";
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
 
  // copy the solution
  double * avec;
  VecGetArray(sol_temp, &avec);

  out.clear(); out.resize(sol_length);
  for(int ii=0; ii<sol_length; ++ii)
    out[ii] = avec[ii]; 

  VecRestoreArray(sol_temp, &avec); 
  VecDestroy(&sol_temp);
}



void Post_probe::readPETSc_average( const int &sol_index_start, 
    const int &sol_index_step, const int &sol_index_end, 
    std::vector<double> &out ) const
{
  std::vector<double> temp;
  const int sol_size = nFunc * dof; 
  out.resize(sol_size);
  for(int jj=0; jj<sol_size; ++jj) out[jj] = 0.0;

  double counter = 0.0;
  for(int ii = sol_index_start; ii <= sol_index_end; ii += sol_index_step)
  {
    std::cout<<"Read solution with index = "<<ii<<std::endl;
    readPETSc_full(ii, temp);
    for(int jj=0; jj<sol_size; ++jj) out[jj] += temp[jj];
    counter += 1.0;
  }
  
  for(int jj=0; jj<sol_size; ++jj) out[jj] = out[jj] / counter;
}



void Post_probe::get_val_cood( const int sol_index,
    const int &ee, const double &xi_s, const double &xi_t,
    FEAElement * const &element,
    std::vector<double> &val, double &xx, double &yy ) const
{
  BernsteinBasis_Array * bs = new BernsteinBasis_Array(s_degree, xi_s);
  BernsteinBasis_Array * bt = new BernsteinBasis_Array(t_degree, xi_t);

  int ecpu = get_elem_cpu(ee); 
 
  HDF5_PartReader * h5r = new HDF5_PartReader( partfile_bname, ecpu );

  // get list of local elements and number of local elements
  std::vector<int> eloc;
  int nlocele;

  h5r->get_LE( eloc, nlocele );
  
  std::vector<int>::const_iterator it;
  
  it = find( eloc.begin(), eloc.end(), ee);
  if( it == eloc.end() )
  {
    std::cerr<<"Error: Failed to find "<<ee<<" element in part file of cpu "<<ecpu<<std::endl;
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  
  int lee = it - eloc.begin(); // lee gives the local index of this element in eloc

  // get the local_to_global mapping for this subdomain
  std::vector<int> l2g;
  int nlgnode;
  h5r->get_L2GN( nlgnode, l2g );

  IALocal_meshSize * locmSize = new ALocal_meshSize_2D_NURBS(h5r);
  FEANode * fNode = new FEANode(h5r);
  IAExtractor * fExt = new AExtractor_2D_NURBS_xy(h5r);

  // mesh size
  double ehx = locmSize->get_hx(lee);
  double ehy = locmSize->get_hy(lee);

  // extractor
  std::vector<double> ext_x, ext_y;
  fExt->get_EXT_x(lee, ext_x);
  fExt->get_EXT_y(lee, ext_y);

  // LIEN
  std::vector<int> locIEN;
  h5r->get_LIEN(lee, locIEN);
  
  if( (int) locIEN.size() != nLocBas )
  {
    std::cerr<<"Error: local IEN format wrong. \n";
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  // control pts
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];

  fNode->get_ctrlPts_xy( nLocBas, &locIEN[0], ectrl_x, ectrl_y );

  // build element basis function value at xi_s and xi_t point
  element->buildBasis( ehx, ehy, bs, bt, ectrl_x, ectrl_y, &ext_x[0], &ext_y[0] );

  // now map LIEN to its global indices
  for( int ii=0; ii<nLocBas; ++ii )
    locIEN[ii] = l2g[locIEN[ii]];

  // read solution vector
  double * local_sol = new double [dof * nLocBas];
  readPETSc_vec( sol_index, locIEN, local_sol );

  double * R = new double [nLocBas];

  element->get_R(0, R);

  val.clear();
  val.resize(dof); 

  for(int ii=0; ii<dof; ++ii)
    val[ii] = 0.0;

  xx = 0.0; yy = 0.0;
  for(int ii=0; ii<nLocBas; ++ii)
  {
    for(int jj=0; jj<dof; ++jj)
      val[jj] += R[ii] * local_sol[dof * ii + jj];
    
    xx  += R[ii] * ectrl_x[ii];
    yy  += R[ii] * ectrl_y[ii];
  }

  delete h5r;
  delete [] local_sol; delete [] R;
  delete [] ectrl_x; delete [] ectrl_y; 
  delete bs; delete bt;
  delete locmSize; delete fNode; delete fExt;
}



void Post_probe::get_lineAve_forNu( const int &sol_index, 
    const std::vector<double> &eta_list,
    const std::vector<int> &ey_list,
    const int &in_nqps, const int &nElem_x, const int &nElem_y,
    FEAElement * const &element,
    std::vector<double> &uzt, std::vector<double> &the, std::vector<double> &yyy ) const
{
  // Read solution vector
  std::vector<double> sol;
  readPETSc_full(sol_index, sol);

  // Check the eat_list and ey_list
  if(eta_list.size() != ey_list.size())
  {
    std::cerr<<"Error: size of eta and ey does not match. \n";
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  
  // Allocate containers
  double * R = new double [nLocBas];
  double * jac = new double [4];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];

  int len = (int) ey_list.size();
  
  uzt.resize(len);
  the.resize(len);
  yyy.resize(len);

  for(int mm = 0; mm<len; ++mm)
  {
    std::cout<<"index "<<mm<<std::endl;
    // We are integrating on the s-line; we sample Gauss points in the s-dir, fix
    // the eta coordinates
    IQuadPts * quad_s = new QuadPts_Gauss(in_nqps);
    std::vector<double> in_qp_t, in_qw_t;
    in_qp_t.push_back(eta_list[mm]);
    in_qw_t.push_back(1.0);
    IQuadPts * quad_t = new QuadPts_debug(1, in_qp_t, in_qw_t);

    AInt_Weight * Int_w = new AInt_Weight(quad_s, quad_t);

    BernsteinBasis_Array * bs = new BernsteinBasis_Array(s_degree, quad_s);
    BernsteinBasis_Array * bt = new BernsteinBasis_Array(t_degree, quad_t);

    double quad_const = 0.0;

    int ey = ey_list[mm]; 
    int offset = 0;
    double theta_ave = 0.0;
    double yy = 0.0, yye = 0.0;
    double uz_theta = 0.0;
    double uz_qua = 0.0, the_qua = 0.0;

    for(int ex = 0; ex < nElem_x; ++ex)
    {
      int ee = ex + ey * nElem_x;

      int ecpu = get_elem_cpu(ee);
      HDF5_PartReader * h5r = new HDF5_PartReader( partfile_bname, ecpu );

      std::vector<int> eloc;
      int nlocele;
      h5r->get_LE(eloc, nlocele);

      std::vector<int>::const_iterator it;
      it = find( eloc.begin(), eloc.end(), ee);
      if(it == eloc.end() )
      {
        std::cerr<<"Error: Failed to find "<<ee<<"th element in part file of cpu "<<ecpu<<std::endl;
        MPI_Abort(PETSC_COMM_WORLD, 1);
      }

      int lee = it - eloc.begin();

      std::vector<int> l2g; int nlgnode;
      h5r->get_L2GN( nlgnode, l2g );

      IALocal_meshSize * locmSize = new ALocal_meshSize_2D_NURBS(h5r);
      FEANode * fNode = new FEANode(h5r);
      IAExtractor * fExt = new AExtractor_2D_NURBS_xy(h5r);

      double ehx = locmSize->get_hx(lee);
      double ehy = locmSize->get_hy(lee);

      std::vector<double> ext_x, ext_y;
      fExt->get_EXT_x(lee, ext_x);
      fExt->get_EXT_y(lee, ext_y);

      std::vector<int> locIEN;
      h5r->get_LIEN(lee, locIEN);

      fNode->get_ctrlPts_xy( nLocBas, &locIEN[0], ectrl_x, ectrl_y );

      element->buildBasis( ehx, ehy, bs, bt, ectrl_x, ectrl_y, &ext_x[0], &ext_y[0] );

      for( int ii=0; ii<nLocBas; ++ii ) locIEN[ii] = l2g[locIEN[ii]];

      // Let me calcuate a base yy value
      if(ex == 0)
      {
        element->get_R(0, R);
        for(int ii=0; ii<nLocBas; ++ii) yy += R[ii] * ectrl_y[ii];
      }

      // Now calculate the surface averages
      for(int ii=0; ii<in_nqps; ++ii)
      {
        element->get_R(ii, R);
        element->get_Jacobian(ii, jac);
        quad_const = jac[0] * ehx * Int_w->get_weight(ii);
        yye = 0.0;          // y-coordinate at quad point ii of element ee
        uz_qua = 0.0;       // u_vertical / theta at qua point ii of element ee
        the_qua = 0.0;      // -1/theta at qua point ii of element ee
        for(int jj=0; jj<nLocBas; ++jj)
        {
          offset = locIEN[jj] * dof;
          yye     += ectrl_y[jj] * R[jj];
          uz_qua  += sol[offset + 2] * R[jj];
          the_qua += sol[offset + 3] * R[jj];
        }

        if( std::abs(yye - yy) > 1.0e-14 )
        {
          std::cerr<<"Error: The element "<<ee<<"'s line is not straight in y-direction. ";
          std::cerr<<yy-yye<<std::endl;
          MPI_Abort(PETSC_COMM_WORLD, 1);
        }

        theta_ave += ( -1.0 / the_qua ) * quad_const;
        uz_theta  += ( uz_qua / (the_qua * the_qua) ) * quad_const;    
      }

      delete h5r; delete locmSize; delete fNode; delete fExt; 
    }
    
    uzt[mm] = uz_theta;
    the[mm] = theta_ave;
    yyy[mm] = yy; 

    delete Int_w; delete quad_t; delete quad_s; delete bs; delete bt;
  }
  delete [] ectrl_x; delete [] ectrl_y; delete [] R; delete [] jac;
}



void Post_probe::get_lineAve_forNu( const int &sol_index_start, 
    const int &sol_index_step, const int &sol_index_end,
    const std::vector<double> &eta_list,
    const std::vector<int> &ey_list,
    const int &in_nqps, const int &nElem_x, const int &nElem_y,
    FEAElement * const &element,
    std::vector<double> &uzt, std::vector<double> &the, std::vector<double> &yyy ) const
{
  // Read solution vector
  std::cout<<"Read the solution vectors ... \n";
  std::vector<double> sol;
  readPETSc_average(sol_index_start, sol_index_step, sol_index_end, sol);

  // Check the eat_list and ey_list
  if(eta_list.size() != ey_list.size())
  {
    std::cerr<<"Error: size of eta and ey does not match. \n";
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  
  // Allocate containers
  double * R = new double [nLocBas];
  double * jac = new double [4];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];

  int len = (int) ey_list.size();
  
  uzt.resize(len);
  the.resize(len);
  yyy.resize(len);

  std::cout<<"Start line average ... \n";
  for(int mm = 0; mm<len; ++mm)
  {
    std::cout<<"index "<<mm<<std::endl;
    // We are integrating on the s-line; we sample Gauss points in the s-dir, fix
    // the eta coordinates
    IQuadPts * quad_s = new QuadPts_Gauss(in_nqps);
    std::vector<double> in_qp_t, in_qw_t;
    in_qp_t.push_back(eta_list[mm]);
    in_qw_t.push_back(1.0);
    IQuadPts * quad_t = new QuadPts_debug(1, in_qp_t, in_qw_t);

    AInt_Weight * Int_w = new AInt_Weight(quad_s, quad_t);

    BernsteinBasis_Array * bs = new BernsteinBasis_Array(s_degree, quad_s);
    BernsteinBasis_Array * bt = new BernsteinBasis_Array(t_degree, quad_t);

    double quad_const = 0.0;

    int ey = ey_list[mm]; 
    int offset = 0;
    double theta_ave = 0.0;
    double yy = 0.0, yye = 0.0;
    double uz_theta = 0.0;
    double uz_qua = 0.0, the_qua = 0.0;

    for(int ex = 0; ex < nElem_x; ++ex)
    {
      int ee = ex + ey * nElem_x;

      int ecpu = get_elem_cpu(ee);
      HDF5_PartReader * h5r = new HDF5_PartReader( partfile_bname, ecpu );

      std::vector<int> eloc;
      int nlocele;
      h5r->get_LE(eloc, nlocele);

      std::vector<int>::const_iterator it;
      it = find( eloc.begin(), eloc.end(), ee);
      if(it == eloc.end() )
      {
        std::cerr<<"Error: Failed to find "<<ee<<"th element in part file of cpu "<<ecpu<<std::endl;
        MPI_Abort(PETSC_COMM_WORLD, 1);
      }

      int lee = it - eloc.begin();

      std::vector<int> l2g; int nlgnode;
      h5r->get_L2GN( nlgnode, l2g );

      IALocal_meshSize * locmSize = new ALocal_meshSize_2D_NURBS(h5r);
      FEANode * fNode = new FEANode(h5r);
      IAExtractor * fExt = new AExtractor_2D_NURBS_xy(h5r);

      double ehx = locmSize->get_hx(lee);
      double ehy = locmSize->get_hy(lee);

      std::vector<double> ext_x, ext_y;
      fExt->get_EXT_x(lee, ext_x);
      fExt->get_EXT_y(lee, ext_y);

      std::vector<int> locIEN;
      h5r->get_LIEN(lee, locIEN);

      fNode->get_ctrlPts_xy( nLocBas, &locIEN[0], ectrl_x, ectrl_y );

      element->buildBasis( ehx, ehy, bs, bt, ectrl_x, ectrl_y, &ext_x[0], &ext_y[0] );

      for( int ii=0; ii<nLocBas; ++ii ) locIEN[ii] = l2g[locIEN[ii]];

      // Let me calcuate a base yy value
      if(ex == 0)
      {
        element->get_R(0, R);
        for(int ii=0; ii<nLocBas; ++ii) yy += R[ii] * ectrl_y[ii];
      }

      // Now calculate the surface averages
      for(int ii=0; ii<in_nqps; ++ii)
      {
        element->get_R(ii, R);
        element->get_Jacobian(ii, jac);
        quad_const = jac[0] * ehx * Int_w->get_weight(ii);
        yye = 0.0;          // y-coordinate at quad point ii of element ee
        uz_qua = 0.0;       // u_vertical / theta at qua point ii of element ee
        the_qua = 0.0;      // -1/theta at qua point ii of element ee
        for(int jj=0; jj<nLocBas; ++jj)
        {
          offset = locIEN[jj] * dof;
          yye     += ectrl_y[jj] * R[jj];
          uz_qua  += sol[offset + 2] * R[jj];
          the_qua += sol[offset + 3] * R[jj];
        }

        if( std::abs(yye - yy) > 1.0e-14 )
        {
          std::cerr<<"Error: The element "<<ee<<"'s line is not straight in y-direction. ";
          std::cerr<<yy-yye<<std::endl;
          MPI_Abort(PETSC_COMM_WORLD, 1);
        }

        theta_ave += ( -1.0 / the_qua ) * quad_const;
        uz_theta  += ( uz_qua / (the_qua * the_qua) ) * quad_const;    
      }

      delete h5r; delete locmSize; delete fNode; delete fExt; 
    }
    
    uzt[mm] = uz_theta;
    the[mm] = theta_ave;
    yyy[mm] = yy; 

    delete Int_w; delete quad_t; delete quad_s; delete bs; delete bt;
  }
  delete [] ectrl_x; delete [] ectrl_y; delete [] R; delete [] jac;
}



void Post_probe::get_lineAve_NSK_rho_theta( const int &sol_index, 
    const std::vector<double> &eta_list,
    const std::vector<int> &ey_list,
    const int &in_nqps, const int &nElem_x, const int &nElem_y,
    FEAElement * const &element,
    std::vector<double> &rho, std::vector<double> &the, 
    std::vector<double> &pre, std::vector<double> &yyy ) const
{
  const double fac8_27 = 8.0 / 27.0;
  
  // Read solution vector
  std::vector<double> sol;
  readPETSc_full(sol_index, sol);

  // Check the eat_list and ey_list
  if( eta_list.size() != ey_list.size() ) 
    SYS_T::print_fatal("Error: size of eta and ey does not match\n");
  
  // Allocate containers
  double * R = new double [nLocBas];
  double * Rx = new double [nLocBas];
  double * Ry = new double [nLocBas];
  double * jac = new double [4];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];

  int len = (int) ey_list.size();
  
  rho.resize(len);
  the.resize(len);
  pre.resize(len);
  yyy.resize(len);

  for(int mm = 0; mm<len; ++mm)
  {
    std::cout<<"index "<<mm<<std::endl;
    // We are integrating on the s-line; we sample Gauss points in the s-dir, fix
    // the eta coordinates
    IQuadPts * quad_s = new QuadPts_Gauss(in_nqps);
    std::vector<double> in_qp_t, in_qw_t;
    in_qp_t.push_back(eta_list[mm]);
    in_qw_t.push_back(1.0);
    IQuadPts * quad_t = new QuadPts_debug(1, in_qp_t, in_qw_t);

    AInt_Weight * Int_w = new AInt_Weight(quad_s, quad_t);

    BernsteinBasis_Array * bs = new BernsteinBasis_Array(s_degree, quad_s);
    BernsteinBasis_Array * bt = new BernsteinBasis_Array(t_degree, quad_t);

    double quad_const = 0.0;

    int ey = ey_list[mm]; 
    int offset = 0;
    double the_ave = 0.0;
    double rho_ave = 0.0;
    double pre_ave = 0.0;
    double yy = 0.0, yye = 0.0;
    double rho_qua = 0.0, the_qua = 0.0;

    for(int ex = 0; ex < nElem_x; ++ex)
    {
      int ee = ex + ey * nElem_x;

      int ecpu = get_elem_cpu(ee);
      HDF5_PartReader * h5r = new HDF5_PartReader( partfile_bname, ecpu );

      std::vector<int> eloc;
      int nlocele;
      h5r->get_LE(eloc, nlocele);

      std::vector<int>::const_iterator it;
      it = find( eloc.begin(), eloc.end(), ee);
      if( it == eloc.end() )
      {
        std::cerr<<"Error: Failed to find "<<ee<<"th element in part file of cpu "<<ecpu<<std::endl;
        MPI_Abort(PETSC_COMM_WORLD, 1);
      }

      int lee = it - eloc.begin();

      std::vector<int> l2g; int nlgnode;
      h5r->get_L2GN( nlgnode, l2g );

      IALocal_meshSize * locmSize = new ALocal_meshSize_2D_NURBS(h5r);
      FEANode * fNode = new FEANode(h5r);
      IAExtractor * fExt = new AExtractor_2D_NURBS_xy(h5r);

      double ehx = locmSize->get_hx(lee);
      double ehy = locmSize->get_hy(lee);

      std::vector<double> ext_x, ext_y;
      fExt->get_EXT_x(lee, ext_x);
      fExt->get_EXT_y(lee, ext_y);

      std::vector<int> locIEN;
      h5r->get_LIEN(lee, locIEN);

      fNode->get_ctrlPts_xy( nLocBas, &locIEN[0], ectrl_x, ectrl_y );

      element->buildBasis( ehx, ehy, bs, bt, ectrl_x, ectrl_y, &ext_x[0], &ext_y[0] );

      for( int ii=0; ii<nLocBas; ++ii ) locIEN[ii] = l2g[locIEN[ii]];

      // Let me calcuate a base yy value
      if(ex == 0)
      {
        element->get_R(0, R);
        for(int ii=0; ii<nLocBas; ++ii) yy += R[ii] * ectrl_y[ii];
      }

      // Now calculate the surface averages
      for(int ii=0; ii<in_nqps; ++ii)
      {
        element->get_R_gradR(ii, R, Rx, Ry);
        element->get_Jacobian(ii, jac);
        quad_const = jac[0] * ehx * Int_w->get_weight(ii);
        yye = 0.0;          // y-coordinate at quad point ii of element ee
        rho_qua = 0.0;      // rho at qua point ii of element ee 
        the_qua = 0.0;      // -1/theta at qua point ii of element ee
        for(int jj=0; jj<nLocBas; ++jj)
        {
          offset   = locIEN[jj] * dof;
          yye     += ectrl_y[jj] * R[jj];
          rho_qua += sol[offset + 0] * R[jj];
          the_qua += sol[offset + 3] * R[jj];
        }
        
        // check if the line is horizontal
        if( std::abs(yye - yy) > 1.0e-14 )
        {
          std::cerr<<"Error: The element "<<ee<<"'s line is not straight in y-direction. ";
          std::cerr<<yy-yye<<std::endl;
          MPI_Abort(PETSC_COMM_WORLD, 1);
        }

        rho_ave += rho_qua * quad_const;    
        the_ave += ( -1.0 / the_qua ) * quad_const;
        pre_ave += (fac8_27 * rho_qua / (the_qua * (rho_qua-1.0)) - rho_qua*rho_qua) * quad_const;
      }

      delete h5r; delete locmSize; delete fNode; delete fExt; 
    }
    
    rho[mm] = rho_ave;
    the[mm] = the_ave;
    pre[mm] = pre_ave;
    yyy[mm] = yy; 

    delete Int_w; delete quad_t; delete quad_s; delete bs; delete bt;
  }
  delete [] ectrl_x; delete [] ectrl_y; delete [] R; delete [] jac;
  delete [] Rx; delete [] Ry;
}

// EOF
