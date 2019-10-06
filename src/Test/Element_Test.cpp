#include "Element_Test.hpp"

void TEST_T::Element_JacxInvJac_check( const std::vector<FEAElement *> &earray,
    const ALocal_Elem * const &locelem,  const IALocal_meshSize * const &msize )
{
  double hx, hy, hz;
  double * id = new double [9];
  const double nElem = locelem->get_nlocalele();

  for(int ee=0; ee<nElem; ++ee)
  {
    hx = msize->get_hx(ee);
    hy = msize->get_hy(ee);
    hz = msize->get_hz(ee);

    if(hx*hy*hz > 0.0)
    {

      double * dx_ds = new double [9];
      double * ds_dx = new double [9];

      for(int ii=0; ii<earray[ee]->get_numQuapts(); ++ii)
      {
        earray[ee]->get_Jacobian(ii, dx_ds);
        earray[ee]->get_invJacobian(ii, ds_dx);

        id[0] = dx_ds[0] * hx * ds_dx[0] + dx_ds[1] * hy * ds_dx[3] + dx_ds[2] * hz * ds_dx[6];
        id[1] = dx_ds[0] * hx * ds_dx[1] + dx_ds[1] * hy * ds_dx[4] + dx_ds[2] * hz * ds_dx[7];
        id[2] = dx_ds[0] * hx * ds_dx[2] + dx_ds[1] * hy * ds_dx[5] + dx_ds[2] * hz * ds_dx[8];

        id[3] = dx_ds[3] * hx * ds_dx[0] + dx_ds[4] * hy * ds_dx[3] + dx_ds[5] * hz * ds_dx[6];
        id[4] = dx_ds[3] * hx * ds_dx[1] + dx_ds[4] * hy * ds_dx[4] + dx_ds[5] * hz * ds_dx[7];
        id[5] = dx_ds[3] * hx * ds_dx[2] + dx_ds[4] * hy * ds_dx[5] + dx_ds[5] * hz * ds_dx[8];

        id[6] = dx_ds[6] * hx * ds_dx[0] + dx_ds[7] * hy * ds_dx[3] + dx_ds[8] * hz * ds_dx[6];
        id[7] = dx_ds[6] * hx * ds_dx[1] + dx_ds[7] * hy * ds_dx[4] + dx_ds[8] * hz * ds_dx[7];
        id[8] = dx_ds[6] * hx * ds_dx[2] + dx_ds[7] * hy * ds_dx[5] + dx_ds[8] * hz * ds_dx[8];


        if(std::abs(id[0]-1.0) > 1.0e-15)
          std::cout<<ee<<'\t'<<ii<<'\t'<<id[0]-1.0<<std::endl;
        if(std::abs(id[1]-0.0) > 1.0e-15)
          std::cout<<ee<<'\t'<<ii<<'\t'<<id[1]<<std::endl;
        if(std::abs(id[2]-0.0) > 1.0e-15)
          std::cout<<ee<<'\t'<<ii<<'\t'<<id[2]<<std::endl;

        if(std::abs(id[4]-1.0) > 1.0e-15)
          std::cout<<ee<<'\t'<<ii<<'\t'<<id[0]-1.0<<std::endl;
        if(std::abs(id[3]-0.0) > 1.0e-15)
          std::cout<<ee<<'\t'<<ii<<'\t'<<id[3]<<std::endl;
        if(std::abs(id[5]-0.0) > 1.0e-15)
          std::cout<<ee<<'\t'<<ii<<'\t'<<id[5]<<std::endl;

        if(std::abs(id[8]- 1.0) > 1.0e-15)
          std::cout<<ee<<'\t'<<ii<<'\t'<<id[0]-1.0<<std::endl;
        if(std::abs(id[6]-0.0) > 1.0e-15)
          std::cout<<ee<<'\t'<<ii<<'\t'<<id[6]<<std::endl;
        if(std::abs(id[7]-0.0) > 1.0e-15)
          std::cout<<ee<<'\t'<<ii<<'\t'<<id[7]<<std::endl;
      }

      delete [] dx_ds; delete [] ds_dx;
    }
  }

  delete [] id;
}



void TEST_T::Element_JacxInvJac_check( FEAElement const * const &ele1 )
{
  const int ele_dim = ele1->get_elemDim();
  const int numQuapts = ele1->get_numQuapts();

  double id; 

  double * dx_ds = new double [ele_dim * ele_dim];

  double * ds_dx = new double [ele_dim * ele_dim];


  for(int qua=0; qua<numQuapts; ++qua)
  {
    ele1->get_Jacobian(qua, dx_ds);
    ele1->get_invJacobian(qua, ds_dx);
    for(int ii=0; ii<ele_dim; ++ii)
    {
      for(int jj=0; jj<ele_dim; ++jj)
      {
        id = 0.0;
        for(int kk=0; kk<ele_dim; ++kk)
        {
          id += dx_ds[ii*ele_dim+kk] * ds_dx[kk*ele_dim + jj];
        }
        if(ii == jj)
        {
          if( std::abs(id-1.0)>1.0e-14 ) std::cout<<id-1.0<<std::endl;
        }
        else
        {
          if(std::abs(id)>1.0e-14) std::cout<<id<<std::endl;
        }
      }
    }
  }

  delete [] dx_ds; delete [] ds_dx;
}


void TEST_T::Element_invJac_compare(
    const std::vector<FEAElement *> &earray1,
    const std::vector<FEAElement *> &earray2,
    const ALocal_Elem * const &locelem,
    const IALocal_meshSize * const &msize )
{
  double hx, hy, hz;
  const double nElem = locelem->get_nlocalele();

  for(int ee=0; ee<nElem; ++ee)
  {
    hx = msize->get_hx(ee);
    hy = msize->get_hy(ee);
    hz = msize->get_hz(ee);
    if( hx * hy * hz > 0.0 )
    {
      double * ds_dx = new double [9];
      double * dxi_dx = new double [9];
      for(int ii=0; ii<earray1[ee]->get_numQuapts(); ++ii)
      {
        earray1[ee]->get_invJacobian(ii, ds_dx);
        earray2[ee]->get_invJacobian(ii, dxi_dx);

        for(int jj=0; jj<9; ++jj)
        {
          if(dxi_dx[jj] != ds_dx[jj])
          {
            std::cout<<"Element "<<ee<<" qua pt "<<ii<<'\t'<<dxi_dx[jj]<<'\t'<<ds_dx[jj]<<'\n';
          }

        }
      }
      delete [] ds_dx; delete [] dxi_dx;
    }
  }
}



void TEST_T::Element_invJac_compare( FEAElement const * const &ele1,
    FEAElement const * const &ele2,
    const double &ehx, const double &ehy, const double &ehz,
    const double &tol )
{
  const int numQuapts = ele1->get_numQuapts();
  double * ds_dx_1 = new double [9];
  double * ds_dx_2 = new double [9];

  for(int qua = 0; qua<numQuapts; ++qua)
  {
    ele1->get_invJacobian(qua, ds_dx_1);
    ele2->get_invJacobian(qua, ds_dx_2);

    /*
       ds_dx_2[0] *= 1.0 / ehx;
       ds_dx_2[1] *= 1.0 / ehx;
       ds_dx_2[2] *= 1.0 / ehx;
       ds_dx_2[3] *= 1.0 / ehy;
       ds_dx_2[4] *= 1.0 / ehy;
       ds_dx_2[5] *= 1.0 / ehy;
       ds_dx_2[6] *= 1.0 / ehz;
       ds_dx_2[7] *= 1.0 / ehz;
       ds_dx_2[8] *= 1.0 / ehz;
       */
    for(int ii=0; ii<9; ++ii)
    {
      if( std::abs( ds_dx_1[ii] - ds_dx_2[ii] ) > tol )
        std::cout<<"Difference "<<std::setprecision(16)<<ds_dx_1[ii]<<'\t'<<ds_dx_2[ii]<<'\n';
    }
  }
  delete [] ds_dx_1; delete [] ds_dx_2;
}


void TEST_T::Element_Jac_compare(  FEAElement const * const &ele1,
    FEAElement const * const &ele2, const double &tol )
{
  const int numQuapts = ele1->get_numQuapts();
  double * dx_ds_1 = new double [9];
  double * dx_ds_2 = new double [9];

  for(int qua = 0; qua<numQuapts; ++qua)
  {
    ele1->get_Jacobian(qua, dx_ds_1);
    ele2->get_Jacobian(qua, dx_ds_2);
    for(int ii=0; ii<9; ++ii)
    {
      if( std::abs( dx_ds_1[ii] - dx_ds_2[ii] ) > tol )
        std::cout<<"Difference "<<std::setprecision(16)<<dx_ds_1[ii]<<'\t'<<dx_ds_2[ii]<<'\n';
    }
  }
  delete [] dx_ds_1; delete [] dx_ds_2;
}



void TEST_T::Element_R_dR_LapR_compare(
    const std::vector<FEAElement *> &earray1,
    const std::vector<FEAElement *> &earray2,
    const IAGlobal_Mesh_Info * const &gmiptr,
    const ALocal_Elem * const &locelem,
    const IALocal_meshSize * const &msize )
{
  double hx, hy, hz;
  const double nElem = locelem->get_nlocalele();
  const int nLocBas = gmiptr->get_nLocBas();

  double * R = new double [nLocBas];
  double * dR_dx = new double [nLocBas];
  double * dR_dy = new double [nLocBas];
  double * dR_dz = new double [nLocBas];
  double * d2R_dxx = new double [nLocBas];
  double * d2R_dyy = new double [nLocBas];
  double * d2R_dzz = new double [nLocBas];

  double * R2 = new double [nLocBas];
  double * dR2_dx = new double [nLocBas];
  double * dR2_dy = new double [nLocBas];
  double * dR2_dz = new double [nLocBas];
  double * d2R2_dxx = new double [nLocBas];
  double * d2R2_dyy = new double [nLocBas];
  double * d2R2_dzz = new double [nLocBas];

  for(int ee=0; ee<nElem; ++ee)
  {
    hx = msize->get_hx(ee);
    hy = msize->get_hy(ee);
    hz = msize->get_hz(ee);

    if( hx * hy * hz > 0.0 )
    {
      for(int ii=0; ii<earray1[ee]->get_numQuapts(); ++ii)
      {
        earray1[ee]->get_3D_R_gradR_LaplacianR(ii, R, dR_dx, dR_dy, dR_dz,
            d2R_dxx, d2R_dyy, d2R_dzz);
        earray2[ee]->get_3D_R_gradR_LaplacianR(ii, R2, dR2_dx, dR2_dy, dR2_dz,
            d2R2_dxx, d2R2_dyy, d2R2_dzz);

        for(int jj=0; jj<nLocBas; ++jj)
        {
          assert(R[jj] == R2[jj]);
          assert(dR_dx[jj] == dR2_dx[jj]);
          assert(dR_dy[jj] == dR2_dy[jj]);
          assert(dR_dz[jj] == dR2_dz[jj]);
          assert(d2R_dxx[jj] == d2R2_dxx[jj]);
          assert(d2R_dyy[jj] == d2R2_dyy[jj]);
          assert(d2R_dzz[jj] == d2R2_dzz[jj]);
        }
        assert(earray1[ee]->get_detJac(ii) == earray2[ee]->get_detJac(ii));
      }
    }
  }

  delete [] R; delete [] dR_dx; delete [] dR_dy; delete [] dR_dz;
  delete [] R2; delete [] dR2_dx; delete [] dR2_dy; delete [] dR2_dz;
  delete [] d2R_dxx; delete[] d2R_dyy; delete [] d2R_dzz;
  delete [] d2R2_dxx; delete[] d2R2_dyy; delete [] d2R2_dzz;
}


void TEST_T::Element_2D_R_dR_d2R_compare( FEAElement const * const &ele1,
    FEAElement const * const &ele2, const double &tol )
{
  const int nLocBas = ele1->get_nLocBas();
  const int numQuapts = ele1->get_numQuapts();

  assert(nLocBas == ele2->get_nLocBas());
  assert(numQuapts == ele2->get_numQuapts());

  double * R1 = new double [nLocBas];
  double * R2 = new double [nLocBas];

  double * dR1x = new double [nLocBas];
  double * dR2x = new double [nLocBas];

  double * dR1y = new double [nLocBas];
  double * dR2y = new double [nLocBas];

  double * dR1xx = new double [nLocBas];
  double * dR2xx = new double [nLocBas];

  double * dR1yy = new double [nLocBas];
  double * dR2yy = new double [nLocBas];

  double * dR1xy = new double [nLocBas];
  double * dR2xy = new double [nLocBas];

  for(int ii=0; ii<numQuapts; ++ii)
  {
    ele1->get_2D_R_dR_d2R(ii, R1, dR1x, dR1y, dR1xx, dR1yy, dR1xy);
    ele2->get_2D_R_dR_d2R(ii, R2, dR2x, dR2y, dR2xx, dR2yy, dR2xy);

    for(int jj=0; jj<nLocBas; ++jj)
    {
      if( !MATH_T::equals(R1[jj], R2[jj], tol) )
        std::cout<<std::setprecision(16)<<"Error R: "<<ii<<'\t'<<jj<<'\t'<<R1[jj]<<'\t'<<R2[jj]<<'\t'<<std::abs(R1[jj] - R2[jj])<<std::endl;
      if( !MATH_T::equals(dR1x[jj], dR2x[jj], tol) )
        std::cout<<std::setprecision(16)<<"Error dR_dx: "<<ii<<'\t'<<jj<<'\t'<<dR1x[jj]<<'\t'<<dR2x[jj]<<'\t'<<std::abs(dR1x[jj] - dR2x[jj])<<std::endl;
      if( !MATH_T::equals(dR1y[jj], dR2y[jj], tol) )
        std::cout<<std::setprecision(16)<<"Error dR_dy: "<<ii<<'\t'<<jj<<'\t'<<dR1y[jj]<<'\t'<<dR2y[jj]<<'\t'<<std::abs(dR1y[jj] - dR2y[jj])<<std::endl;
      if( !MATH_T::equals(dR1xx[jj], dR2xx[jj], tol) )
        std::cout<<std::setprecision(16)<<"Error d2R_dxx: "<<ii<<'\t'<<jj<<'\t'<<dR1xx[jj]<<'\t'<<dR2xx[jj]<<'\t'<<std::abs(dR1xx[jj] - dR2xx[jj])<<std::endl;
      if( !MATH_T::equals(dR1yy[jj], dR2yy[jj], tol) )
        std::cout<<std::setprecision(16)<<"Error d2R_dyy: "<<ii<<'\t'<<jj<<'\t'<<dR1yy[jj]<<'\t'<<dR2yy[jj]<<'\t'<<std::abs(dR1yy[jj] - dR2yy[jj])<<std::endl;
      if( !MATH_T::equals(dR1xy[jj], dR2xy[jj], tol) )
        std::cout<<std::setprecision(16)<<"Error d2R_dxy: "<<ii<<'\t'<<jj<<'\t'<<dR1xy[jj]<<'\t'<<dR2xy[jj]<<'\t'<<std::abs(dR1xy[jj] - dR2xy[jj])<<std::endl;
    }
    if( !MATH_T::equals(ele1->get_detJac(ii), ele2->get_detJac(ii), tol) )
      std::cout<<std::setprecision(16)<<"detJac: "<<ii<<'\t'<<ele1->get_detJac(ii) - ele2->get_detJac(ii)<<'\n';
  }


  delete [] R1; delete [] R2; 
  delete [] dR1x; delete [] dR2x;
  delete [] dR1y; delete [] dR2y; 
  delete [] dR1xx; delete [] dR2xx; 
  delete [] dR1yy; delete [] dR2yy;
  delete [] dR1xy; delete [] dR2xy;
}


void TEST_T::Element_R_compare( FEAElement const * const &ele1,
    FEAElement const * const &ele2, const double &tol )
{
  const int nLocBas = ele1->get_nLocBas();
  const int numQuapts = ele1->get_numQuapts();

  assert(nLocBas == ele2->get_nLocBas());
  assert(numQuapts == ele2->get_numQuapts());

  double * R1 = new double [nLocBas];
  double * R2 = new double [nLocBas];

  for(int ii=0; ii<numQuapts; ++ii)
  {
    ele1->get_R(ii, R1);
    ele2->get_R(ii, R2);

    for(int jj=0; jj<nLocBas; ++jj)
    {
      if(std::abs( R1[jj] -R2[jj] ) > tol )
        std::cout<<std::setprecision(16)<<"error R: "<<ii<<'\t'<<jj<<'\t'<<R1[jj]<<'\t'<<R2[jj]<<'\t'<<std::abs(R1[jj] - R2[jj])<<std::endl;
    }

    if(std::abs(ele1->get_detJac(ii) - ele2->get_detJac(ii)) > tol)
    {
      std::cout<<std::setprecision(16)<<"difference in detJac : "<<ii<<'\t'<<ele1->get_detJac(ii) 
        - ele2->get_detJac(ii)<<'\n';
    }
  }

  delete [] R1;
  delete [] R2;
}



void TEST_T::Element_grad2d_compare( FEAElement const * const &ele1,
    FEAElement const * const &ele2, const double &tol )
{
  const int nLocBas = ele1->get_nLocBas();
  const int numQuapts = ele1->get_numQuapts();

  assert(nLocBas == ele2->get_nLocBas());
  assert(numQuapts == ele2->get_numQuapts());

  double * R1 = new double [nLocBas];
  double * R2 = new double [nLocBas];

  double * R1x = new double [nLocBas];
  double * R2x = new double [nLocBas];

  double * R1y = new double [nLocBas];
  double * R2y = new double [nLocBas];

  for(int ii=0; ii<numQuapts; ++ii)
  {
    ele1->get_R_gradR(ii, R1, R1x, R1y);
    ele2->get_R_gradR(ii, R2, R2x, R2y);

    for(int jj=0; jj<nLocBas; ++jj)
    {
      if(std::abs( R1x[jj] -R2x[jj] ) > tol )
        std::cout<<std::setprecision(16)<<"error Rx: "<<ii<<'\t'<<jj<<'\t'<<R1x[jj]<<'\t'<<R2x[jj]<<'\t'<<std::abs(R1x[jj] - R2x[jj])<<std::endl;
      if(std::abs( R1y[jj] -R2y[jj] ) > tol )
        std::cout<<std::setprecision(16)<<"error Ry : "<<ii<<'\t'<<jj<<'\t'<<R1y[jj]<<'\t'<<R2y[jj]<<'\t'<<std::abs(R1y[jj] - R2y[jj])<<std::endl;
    }
  }

  delete [] R1;  delete [] R2;
  delete [] R1x;  delete [] R2x;
  delete [] R1y;  delete [] R2y;
}


void TEST_T::Element_grad3d_compare( FEAElement const * const &ele1,
    FEAElement const * const &ele2, const double &tol )
{
  const int nLocBas = ele1->get_nLocBas();
  const int numQuapts = ele1->get_numQuapts();

  assert(nLocBas == ele2->get_nLocBas());
  assert(numQuapts == ele2->get_numQuapts());

  double * R1 = new double [nLocBas];
  double * R2 = new double [nLocBas];

  double * R1x = new double [nLocBas];
  double * R2x = new double [nLocBas];

  double * R1y = new double [nLocBas];
  double * R2y = new double [nLocBas];

  double * R1z = new double [nLocBas];
  double * R2z = new double [nLocBas];

  for(int ii=0; ii<numQuapts; ++ii)
  {
    ele1->get_R(ii, R1);
    ele1->get_gradR(ii, R1x, R1y, R1z);
    ele2->get_R_gradR(ii, R2, R2x, R2y, R2z);

    for(int jj=0; jj<nLocBas; ++jj)
    {
      if(std::abs( R1x[jj] -R2x[jj] ) > tol )
        std::cout<<std::setprecision(16)<<"error : "<<ii<<'\t'<<jj<<'\t'<<R1x[jj]<<'\t'<<R2x[jj]<<std::endl;
      if(std::abs( R1y[jj] -R2y[jj] ) > tol )
        std::cout<<std::setprecision(16)<<"error : "<<ii<<'\t'<<jj<<'\t'<<R1y[jj]<<'\t'<<R2y[jj]<<std::endl;
      if(std::abs( R1z[jj] -R2z[jj] ) > tol )
        std::cout<<std::setprecision(16)<<"error : "<<ii<<'\t'<<jj<<'\t'<<R1z[jj]<<'\t'<<R2z[jj]<<std::endl;
    }
  }

  delete [] R1;  delete [] R2;  delete [] R1x;  delete [] R2x;
  delete [] R1y;  delete [] R2y; delete [] R1z;  delete [] R2z;
}


void TEST_T::Element_d2R_compare( FEAElement const * const &nurbs,
    FEAElement const * const &bs, const double &tol )
{
  const int nLocBas = nurbs->get_nLocBas();
  const int numQuapts = nurbs->get_numQuapts();

  assert(nLocBas == bs->get_nLocBas());
  assert(numQuapts == bs->get_numQuapts());

  double * R = new double [nLocBas];

  double * Rx = new double [nLocBas];
  double * Ry = new double [nLocBas];

  double * Rxx = new double [nLocBas];
  double * Ryy = new double [nLocBas];
  double * Rxy = new double [nLocBas];

  double * Nxx = new double [nLocBas];
  double * Nyy = new double [nLocBas];
  double * Nxy = new double [nLocBas];

  for(int ii=0; ii<numQuapts; ++ii)
  {
    bs->get_2D_R_dR_d2R(ii, R, Rx, Ry, Rxx, Ryy, Rxy);

    for(int jj=0; jj<nLocBas; ++jj)
    {
      if(std::abs( Rxx[jj] - Nxx[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error Rxx: "<<ii<<'\t'<<jj<<'\t'<<Rxx[jj]<<'\t'<<Nxx[jj]<<'\t'<<std::abs( Rxx[jj] - Nxx[jj] )<<std::endl;
      if(std::abs( Ryy[jj] - Nyy[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error Ryy: "<<ii<<'\t'<<jj<<'\t'<<Ryy[jj]<<'\t'<<Nyy[jj]<<'\t'<<std::abs( Ryy[jj] - Nyy[jj] )<<std::endl;
      if(std::abs( Rxy[jj] - Nxy[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error Rxy: "<<ii<<'\t'<<jj<<'\t'<<Rxy[jj]<<'\t'<<Nxy[jj]<<'\t'<<std::abs( Rxy[jj] - Nxy[jj] )<<std::endl;
    }
  }

  delete [] R; delete [] Rx; delete [] Ry; delete [] Rxx;
  delete [] Ryy; delete [] Rxy; delete [] Nxx; delete [] Nyy; delete [] Nxy;
}


void TEST_T::Element_d2R_3d_compare( FEAElement const * const &nurbs,
    FEAElement const * const &bs, const double &tol )
{
  const int nLocBas = nurbs->get_nLocBas();
  const int numQuapts = nurbs->get_numQuapts();

  assert(nLocBas == bs->get_nLocBas());
  assert(numQuapts == bs->get_numQuapts());

  double * R1 = new double [nLocBas];
  double * R1x = new double [nLocBas];
  double * R1y = new double [nLocBas];
  double * R1z = new double [nLocBas];
  double * R1xx = new double [nLocBas];
  double * R1yy = new double [nLocBas];
  double * R1zz = new double [nLocBas];
  double * R1xy = new double [nLocBas];
  double * R1xz = new double [nLocBas];
  double * R1yz = new double [nLocBas];

  double * R2 = new double [nLocBas];
  double * R2x = new double [nLocBas];
  double * R2y = new double [nLocBas];
  double * R2z = new double [nLocBas];
  double * R2xx = new double [nLocBas];
  double * R2yy = new double [nLocBas];
  double * R2zz = new double [nLocBas];
  double * R2xy = new double [nLocBas];
  double * R2xz = new double [nLocBas];
  double * R2yz = new double [nLocBas];

  for(int ii=0; ii<numQuapts; ++ii)
  {
    nurbs->get_3D_R_dR_d2R(ii, R1, R1x, R1y, R1z, R1xx, R1yy, R1zz, R1xy, R1xz, R1yz);
    bs->get_3D_R_dR_d2R(ii, R2, R2x, R2y, R2z, R2xx, R2yy, R2zz, R2xy, R2xz, R2yz);

    for(int jj=0; jj<nLocBas; ++jj)
    {
      if(std::abs( R1xx[jj] - R2xx[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error dxx : "<<ii<<'\t'<<jj<<'\t'<<R1xx[jj]<<'\t'<<R2xx[jj]<<std::endl;
      if(std::abs( R1yy[jj] - R2yy[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error dyy: "<<ii<<'\t'<<jj<<'\t'<<R1yy[jj]<<'\t'<<R2yy[jj]<<std::endl;
      if(std::abs( R1zz[jj] - R2zz[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error dzz: "<<ii<<'\t'<<jj<<'\t'<<R1zz[jj]<<'\t'<<R2zz[jj]<<std::endl;
      if(std::abs( R1xy[jj] - R2xy[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error dxy: "<<ii<<'\t'<<jj<<'\t'<<R1xy[jj]<<'\t'<<R2xy[jj]<<std::endl;
      if(std::abs( R1xz[jj] - R2xz[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error dxz: "<<ii<<'\t'<<jj<<'\t'<<R1xz[jj]<<'\t'<<R2xz[jj]<<std::endl;
      if(std::abs( R1yz[jj] - R2yz[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error dyz: "<<ii<<'\t'<<jj<<'\t'<<R1yz[jj]<<'\t'<<R2yz[jj]<<std::endl;
    }

    nurbs->get_3D_R_gradR_LaplacianR(ii, R1, R1x, R1y, R1z, R1xx, R1yy, R1zz);
    //nurbs->get_3D_R_gradR_LaplacianR(ii, R1, R1x, R1y, R1z, R1xx, R1yy, R1zz);
    bs->get_3D_R_gradR_LaplacianR(ii, R2, R2x, R2y, R2z, R2xx, R2yy, R2zz);

    for(int jj=0; jj<nLocBas; ++jj)
    {
      if(std::abs( R1[jj] - R2[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error R: "<<ii<<'\t'<<jj<<'\t'<<R1[jj]<<'\t'<<R2[jj]<<std::endl;
      if(std::abs( R1x[jj] - R2x[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error dx: "<<ii<<'\t'<<jj<<'\t'<<R1x[jj]<<'\t'<<R2x[jj]<<std::endl;
      if(std::abs( R1y[jj] - R2y[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error dy: "<<ii<<'\t'<<jj<<'\t'<<R1y[jj]<<'\t'<<R2y[jj]<<std::endl;
      if(std::abs( R1z[jj] - R2z[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error dz: "<<ii<<'\t'<<jj<<'\t'<<R1z[jj]<<'\t'<<R2z[jj]<<std::endl;
      if(std::abs( R1xx[jj] - R2xx[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error dxx: "<<ii<<'\t'<<jj<<'\t'<<R1xx[jj]<<'\t'<<R2xx[jj]<<std::endl;
      if(std::abs( R1yy[jj] - R2yy[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error dyy: "<<ii<<'\t'<<jj<<'\t'<<R1yy[jj]<<'\t'<<R2yy[jj]<<std::endl;
      if(std::abs( R1zz[jj] - R2zz[jj] ) > tol)
        std::cout<<std::setprecision(16)<<"Error dzz: "<<ii<<'\t'<<jj<<'\t'<<R1zz[jj]<<'\t'<<R2zz[jj]<<std::endl;
    }

  }

  delete [] R1; delete [] R1x; delete [] R1y; delete [] R1z; delete [] R1xx;
  delete [] R1yy; delete [] R1zz; delete [] R1xy; delete [] R1xz; delete [] R1yz;
  delete [] R2; delete [] R2x; delete [] R2y; delete [] R2z; delete [] R2xx;
  delete [] R2yy; delete [] R2zz; delete [] R2xy; delete [] R2xz; delete [] R2yz;
}



void TEST_T::Element_normal_vector_compare( FEAElement const * const &bs1,
    FEAElement const * const &bs2, const double &tol )
{
  const int nLocBas = bs1->get_nLocBas();
  const int numQuapts = bs1->get_numQuapts();

  assert(nLocBas == bs2->get_nLocBas());
  assert(numQuapts == bs2->get_numQuapts());

  double nx1, nx2, ny1, ny2, line1, line2; 
  for(int ii=0; ii<numQuapts; ++ii)
  {
    bs1->get_2d_normal_front(ii, nx1, ny1, line1);
    bs2->get_2d_normal_front(ii, nx2, ny2, line2);

    if(std::abs(nx1 - nx2) > tol)
      std::cout<<std::setprecision(16)<<"Error : front "<<nx1<<'\t'<<nx2<<std::endl;
    if(std::abs(ny1 - ny2) > tol)
      std::cout<<std::setprecision(16)<<"Error : front "<<ny1<<'\t'<<ny2<<std::endl;
    if(std::abs(line1 - line2) > tol)
      std::cout<<std::setprecision(16)<<"Error : front "<<line1<<'\t'<<line2<<std::endl;

    bs1->get_2d_normal_back(ii, nx1, ny1, line1);
    bs2->get_2d_normal_back(ii, nx2, ny2, line2);

    if(std::abs(nx1 - nx2) > tol)
      std::cout<<std::setprecision(16)<<"Error : back "<<nx1<<'\t'<<nx2<<std::endl;
    if(std::abs(ny1 - ny2) > tol)
      std::cout<<std::setprecision(16)<<"Error : back "<<ny1<<'\t'<<ny2<<std::endl;
    if(std::abs(line1 - line2) > tol)
      std::cout<<std::setprecision(16)<<"Error : back "<<line1<<'\t'<<line2<<std::endl;

    bs1->get_2d_normal_left(ii, nx1, ny1, line1);
    bs2->get_2d_normal_left(ii, nx2, ny2, line2);

    if(std::abs(nx1 - nx2) > tol)
      std::cout<<std::setprecision(16)<<"Error : left "<<nx1<<'\t'<<nx2<<std::endl;
    if(std::abs(ny1 - ny2) > tol)
      std::cout<<std::setprecision(16)<<"Error : left "<<ny1<<'\t'<<ny2<<std::endl;
    if(std::abs(line1 - line2) > tol)
      std::cout<<std::setprecision(16)<<"Error : left "<<line1<<'\t'<<line2<<std::endl;

    bs1->get_2d_normal_right(ii, nx1, ny1, line1);
    bs2->get_2d_normal_right(ii, nx2, ny2, line2);

    if(std::abs(nx1 - nx2) > tol)
      std::cout<<std::setprecision(16)<<"Error : right "<<nx1<<'\t'<<nx2<<std::endl;
    if(std::abs(ny1 - ny2) > tol)
      std::cout<<std::setprecision(16)<<"Error : right "<<ny1<<'\t'<<ny2<<std::endl;
    if(std::abs(line1 - line2) > tol)
      std::cout<<std::setprecision(16)<<"Error : right "<<line1<<'\t'<<line2<<std::endl;
  }
}


void TEST_T::Element_normal_vector_3d_compare( FEAElement const * const &bs1,
    FEAElement const * const &bs2, const double &tol )
{
  const int nLocBas = bs1->get_nLocBas();
  const int numQuapts = bs1->get_numQuapts();

  assert(nLocBas == bs2->get_nLocBas());
  assert(numQuapts == bs2->get_numQuapts());

  double nx1, nx2, ny1, ny2, nz1, nz2, line1, line2; 
  for(int ii=0; ii<numQuapts; ++ii)
  {
    bs1->get_3d_normal_bottom(ii, nx1, ny1, nz1, line1);
    bs2->get_3d_normal_bottom(ii, nx2, ny2, nz2, line2);

    if(std::abs(nx1 - nx2) > tol && nx1 > tol)
      std::cout<<std::setprecision(16)<<"Error : bottom "<<nx1<<'\t'<<nx2<<std::endl;
    if(std::abs(ny1 - ny2) > tol && ny1 > tol)
      std::cout<<std::setprecision(16)<<"Error : bottom "<<ny1<<'\t'<<ny2<<std::endl;
    if(std::abs(nz1 - nz2) > tol && nz1 > tol)
      std::cout<<std::setprecision(16)<<"Error : bottom "<<nz1<<'\t'<<nz2<<std::endl;
    if(std::abs(line1 - line2) > tol)
      std::cout<<std::setprecision(16)<<"Error : bottom "<<line1<<'\t'<<line2<<std::endl;

    bs1->get_3d_normal_top(ii, nx1, ny1, nz1, line1);
    bs2->get_3d_normal_top(ii, nx2, ny2, nz2, line2);

    if(std::abs(nx1 - nx2) > tol && nx1 > tol)
      std::cout<<std::setprecision(16)<<"Error : top "<<nx1<<'\t'<<nx2<<std::endl;
    if(std::abs(ny1 - ny2) > tol && ny1 > tol)
      std::cout<<std::setprecision(16)<<"Error : top "<<ny1<<'\t'<<ny2<<std::endl;
    if(std::abs(nz1 - nz2) > tol && nz1 > tol)
      std::cout<<std::setprecision(16)<<"Error : top "<<nz1<<'\t'<<nz2<<std::endl;
    if(std::abs(line1 - line2) > tol)
      std::cout<<std::setprecision(16)<<"Error : top "<<line1<<'\t'<<line2<<std::endl;

    bs1->get_3d_normal_left(ii, nx1, ny1, nz1, line1);
    bs2->get_3d_normal_left(ii, nx2, ny2, nz2, line2);

    if(std::abs(nx1 - nx2) > tol)
      std::cout<<std::setprecision(16)<<"Error : left "<<nx1<<'\t'<<nx2<<std::endl;
    if(std::abs(ny1 - ny2) > tol)
      std::cout<<std::setprecision(16)<<"Error : left "<<ny1<<'\t'<<ny2<<std::endl;
    if(std::abs(nz1 - nz2) > tol)
      std::cout<<std::setprecision(16)<<"Error : left "<<nz1<<'\t'<<nz2<<std::endl;
    if(std::abs(line1 - line2) > tol)
      std::cout<<std::setprecision(16)<<"Error : left "<<line1<<'\t'<<line2<<std::endl;

    bs1->get_3d_normal_right(ii, nx1, ny1, nz1, line1);
    bs2->get_3d_normal_right(ii, nx2, ny2, nz2, line2);

    if(std::abs(nx1 - nx2) > tol)
      std::cout<<std::setprecision(16)<<"Error : right "<<nx1<<'\t'<<nx2<<std::endl;
    if(std::abs(ny1 - ny2) > tol)
      std::cout<<std::setprecision(16)<<"Error : right "<<ny1<<'\t'<<ny2<<std::endl;
    if(std::abs(nz1 - nz2) > tol)
      std::cout<<std::setprecision(16)<<"Error : right "<<nz1<<'\t'<<nz2<<std::endl;
    if(std::abs(line1 - line2) > tol)
      std::cout<<std::setprecision(16)<<"Error : right "<<line1<<'\t'<<line2<<std::endl;

    bs1->get_3d_normal_front(ii, nx1, ny1, nz1, line1);
    bs2->get_3d_normal_front(ii, nx2, ny2, nz2, line2);

    if(std::abs(nx1 - nx2) > tol)
      std::cout<<std::setprecision(16)<<"Error : front "<<nx1<<'\t'<<nx2<<std::endl;
    if(std::abs(ny1 - ny2) > tol)
      std::cout<<std::setprecision(16)<<"Error : front "<<ny1<<'\t'<<ny2<<std::endl;
    if(std::abs(nz1 - nz2) > tol)
      std::cout<<std::setprecision(16)<<"Error : front "<<nz1<<'\t'<<nz2<<std::endl;
    if(std::abs(line1 - line2) > tol)
      std::cout<<std::setprecision(16)<<"Error : front "<<line1<<'\t'<<line2<<std::endl;

    bs1->get_3d_normal_back(ii, nx1, ny1, nz1, line1);
    bs2->get_3d_normal_back(ii, nx2, ny2, nz2, line2);

    if(std::abs(nx1 - nx2) > tol)
      std::cout<<std::setprecision(16)<<"Error : back "<<nx1<<'\t'<<nx2<<std::endl;
    if(std::abs(ny1 - ny2) > tol)
      std::cout<<std::setprecision(16)<<"Error : back "<<ny1<<'\t'<<ny2<<std::endl;
    if(std::abs(nz1 - nz2) > tol)
      std::cout<<std::setprecision(16)<<"Error : back "<<nz1<<'\t'<<nz2<<std::endl;
    if(std::abs(line1 - line2) > tol)
      std::cout<<std::setprecision(16)<<"Error : back "<<line1<<'\t'<<line2<<std::endl;
  }
}



// EOF
