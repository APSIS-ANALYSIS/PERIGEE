#include "Matrix_double_6by6_Array.hpp"

Matrix_double_6by6_Array::Matrix_double_6by6_Array(const double &aa, const double &bb,
    const double &cc, const double &dd, const double &ee, const double &ff,
    const double &gg, const double &hh, const double &ii)
{
  Mat[0] = aa*aa;     Mat[1] = dd*dd;     Mat[2] = gg*gg;
  Mat[3] = 2.0*aa*dd; Mat[4] = 2.0*aa*gg; Mat[5] = 2.0*dd*gg;

  Mat[6] = aa*bb;       Mat[7]  = dd*ee;       Mat[8]  = gg*hh;
  Mat[9] = aa*ee+bb*dd; Mat[10] = aa*hh+bb*gg; Mat[11] = dd*hh+ee*gg;
  
  Mat[12] = aa*cc;       Mat[13] = dd*ff;       Mat[14] = gg*ii;
  Mat[15] = aa*ff+cc*dd; Mat[16] = aa*ii+cc*gg; Mat[17] = dd*ii+ff*gg;
  
  Mat[18] = bb*bb;     Mat[19] = ee*ee;     Mat[20] = hh*hh;
  Mat[21] = 2.0*bb*ee; Mat[22] = 2.0*bb*hh; Mat[23] = 2.0*ee*hh;
  
  Mat[24] = bb*cc;       Mat[25] = ee*ff;       Mat[26] = hh*ii;
  Mat[27] = bb*ff+cc*ee; Mat[28] = bb*ii+cc*hh; Mat[29] = ee*ii+ff*hh;

  Mat[30] = cc*cc;     Mat[31] = ff*ff;     Mat[32] = ii*ii;
  Mat[33] = 2.0*cc*ff; Mat[34] = 2.0*cc*ii; Mat[35] = 2.0*ff*ii;

  pp[0] = 0; pp[1] = 1; pp[2] = 2; pp[3] = 3; pp[4] = 4; pp[5] = 5;
}

Matrix_double_6by6_Array::~Matrix_double_6by6_Array()
{}

void Matrix_double_6by6_Array::LU_fac()
{
  double max_value, temp;
  int max_index, int_temp;
  bool pivot_flag;
  for(int kk=0; kk<5; ++kk)
  {
    max_value = std::abs( Mat[6*kk+kk] );
    max_index = kk;
    pivot_flag = false;
    // find the column pivoting
    for(int ii=kk+1; ii<6; ++ii)
    {
      if( max_value < std::abs( Mat[6*ii+kk] ) )
      {
        max_value = std::abs( Mat[6*ii+kk] );
        max_index = ii;
        pivot_flag = true;
      }
    }
    
    if(pivot_flag)
    {
      int_temp = pp[kk];
      pp[kk] = pp[max_index];
      pp[max_index] = int_temp;
      
      for(int ii=0; ii<6; ++ii)
      {
        temp = Mat[6*kk+ii];
        Mat[6*kk+ii] = Mat[6*max_index+ii];
        Mat[6*max_index+ii] = temp;
      }
    } 
    
    for(int ii=kk+1; ii<6; ++ii)
    {
      Mat[6*ii+kk] = Mat[6*ii+kk] / Mat[6*kk+kk];
      for(int jj=kk+1; jj<6; ++jj)
        Mat[6*ii+jj] = Mat[6*ii+jj] - Mat[6*ii+kk] * Mat[6*kk+jj];
    }
  }
}

std::array<double, 6> Matrix_double_6by6_Array::LU_solve( const std::array<double, 6> &rhs ) const
{
  std::array<double, 6> xx {{ rhs[pp[0]], rhs[pp[1]], rhs[pp[2]], rhs[pp[3]], rhs[pp[4]], rhs[pp[5]] }};

  xx[1] = xx[1] - Mat[6]  * xx[0];
  xx[2] = xx[2] - Mat[12] * xx[0] - Mat[13] * xx[1];
  xx[3] = xx[3] - Mat[18] * xx[0] - Mat[19] * xx[1] - Mat[20] * xx[2];
  xx[4] = xx[4] - Mat[24] * xx[0] - Mat[25] * xx[1] - Mat[26] * xx[2] - Mat[27] * xx[3];
  xx[5] = xx[5] - Mat[30] * xx[0] - Mat[31] * xx[1] - Mat[32] * xx[2] - Mat[33] * xx[3] - Mat[34] * xx[4];

  xx[5] =  xx[5] / Mat[35];
  xx[4] = (xx[4] - Mat[29] * xx[5]) / Mat[28];
  xx[3] = (xx[3] - Mat[23] * xx[5]  - Mat[22] * xx[4]) / Mat[21];
  xx[2] = (xx[2] - Mat[17] * xx[5]  - Mat[16] * xx[4]  - Mat[15] * xx[3]) / Mat[14];
  xx[1] = (xx[1] - Mat[11] * xx[5]  - Mat[10] * xx[4]  - Mat[9]  * xx[3]  - Mat[8] * xx[2]) / Mat[7];
  xx[0] = (xx[0] - Mat[5]  * xx[5]  - Mat[4]  * xx[4]  - Mat[3]  * xx[3]  - Mat[2] * xx[2]  - Mat[1] * xx[1]) / Mat[0];

  return xx;
}

void Matrix_double_6by6_Array::print() const
{
  for(int ii=0; ii<6; ++ii)
  {
    for(int jj=0; jj<6; ++jj)
      std::cout<<Mat[6*ii+jj]<<'\t';
    std::cout<<'\n';
  }
  std::cout<<'\n';
  for(int ii=0; ii<6; ++ii)
    std::cout<<pp[ii]<<'\t';
  
  std::cout<<'\n'<<std::endl;
}

// EOF
