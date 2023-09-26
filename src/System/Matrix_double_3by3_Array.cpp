#include "Matrix_double_3by3_Array.hpp"

Matrix_double_3by3_Array::Matrix_double_3by3_Array()
{
  mat[0] = 1.0; mat[1] = 0.0; mat[2] = 0.0;
  mat[3] = 0.0; mat[4] = 1.0; mat[5] = 0.0;
  mat[6] = 0.0; mat[7] = 0.0; mat[8] = 1.0;

  pp[0] = 0; pp[1] = 1; pp[2] = 2;
}

Matrix_double_3by3_Array::Matrix_double_3by3_Array( 
    const double &a11, const double &a12, const double &a13,
    const double &a21, const double &a22, const double &a23,
    const double &a31, const double &a32, const double &a33 )
{
  mat[0] = a11;  mat[1] = a12;  mat[2] = a13;
  mat[3] = a21;  mat[4] = a22;  mat[5] = a23;
  mat[6] = a31;  mat[7] = a32;  mat[8] = a33;

  pp[0] = 0; pp[1] = 1; pp[2] = 2; 
}

Matrix_double_3by3_Array::~Matrix_double_3by3_Array()
{}

Matrix_double_3by3_Array& Matrix_double_3by3_Array::operator= (
    const Matrix_double_3by3_Array &input )
{
  if(this != &input)
  {
    for(int ii=0; ii<9; ++ii) mat[ii] = input.mat[ii];
    pp[0] = input.pp[0];
    pp[1] = input.pp[1];
    pp[2] = input.pp[2];
  }
  return *this;
}

void Matrix_double_3by3_Array::gen_id()
{
  mat[0] = 1.0;  mat[1] = 0.0;  mat[2] = 0.0;
  mat[3] = 0.0;  mat[4] = 1.0;  mat[5] = 0.0;
  mat[6] = 0.0;  mat[7] = 0.0;  mat[8] = 1.0;

  pp[0] = 0; pp[1] = 1; pp[2] = 2;
}

void Matrix_double_3by3_Array::gen_rand(const double &min, const double &max)
{
  std::random_device rd;
  std::mt19937_64 gen( rd() );
  std::uniform_real_distribution<double> dis(min, max);
  for(int ii=0; ii<9; ++ii) mat[ii] = dis(gen); 

  pp[0] = 0; pp[1] = 1; pp[2] = 2;
}

void Matrix_double_3by3_Array::gen_hilb()
{
  for(int ii=0; ii<3; ++ii)
    for(int jj=0; jj<3; ++jj)
      mat[ii*3+jj] = 1.0 / (ii + jj + 1.0);

  pp[0] = 0; pp[1] = 1; pp[2] = 2;
}

void Matrix_double_3by3_Array::LU_fac()
{
  double temp;
  int int_temp;
  // 1st row
  double max_value = std::abs(mat[0]);
  int max_index = 0;
  bool pivot_flag = false;

  // Pivoting the first column
  if( max_value < std::abs(mat[3]) )
  {
    max_value = std::abs(mat[3]);
    max_index = 1;
    pivot_flag = true;
  }

  if( max_value < std::abs(mat[6]) )
  {
    max_value = std::abs(mat[6]);
    max_index = 2;
    pivot_flag = true;
  }

  if( pivot_flag )
  {
    int_temp = pp[0];
    pp[0] = pp[max_index];
    pp[max_index] = int_temp;

    temp = mat[0];
    mat[0] = mat[max_index*3];
    mat[max_index*3] = temp;

    temp = mat[1];
    mat[1] = mat[max_index*3+1];
    mat[max_index*3+1] = temp;

    temp = mat[2];
    mat[2] = mat[max_index*3+2];
    mat[max_index*3+2] = temp;
  }

  mat[3] = mat[3] / mat[0];
  mat[4] = mat[4] - mat[3] * mat[1];
  mat[5] = mat[5] - mat[3] * mat[2];
  mat[6] = mat[6] / mat[0];
  mat[7] = mat[7] - mat[6] * mat[1];
  mat[8] = mat[8] - mat[6] * mat[2];


  // Pivoting the second column
  max_value = std::abs(mat[4]);
  max_index = 1;
  pivot_flag = false;

  if(max_value< std::abs(mat[7]))
  {
    max_value = std::abs(mat[7]);
    max_index = 2;
    pivot_flag = true;
  }

  if(pivot_flag)
  {
    int_temp = pp[1];
    pp[1] = pp[2];
    pp[2] = int_temp;

    temp = mat[3];
    mat[3] = mat[6];
    mat[6] = temp;

    temp = mat[4];
    mat[4] = mat[7];
    mat[7] = temp;

    temp = mat[5];
    mat[5] = mat[8];
    mat[8] = temp;
  }


  mat[7] = mat[7] / mat[4];
  mat[8] = mat[8] - mat[7] * mat[5];
}

Vector_3 Matrix_double_3by3_Array::LU_solve( const Vector_3 &bb ) const
{
  Vector_3 xx( bb(pp[0]), bb(pp[1]), bb(pp[2]) );

  xx(1) = xx(1) - mat[3] * xx(0);
  xx(2) = xx(2) - mat[6] * xx(0) - mat[7] * xx(1);

  xx(2) =  xx(2) / mat[8];
  xx(1) = (xx(1) - mat[5] * xx(2)) / mat[4];
  xx(0) = (xx(0) - mat[2] * xx(2)  - mat[1] * xx(1)) / mat[0];

  return xx;
}

std::array<double, 3> Matrix_double_3by3_Array::LU_solve( const std::array<double, 3> &bb ) const
{
  std::array<double, 3> xx {{ bb[pp[0]], bb[pp[1]], bb[pp[2]] }};
  
  xx[1] = xx[1] - mat[3] * xx[0];
  xx[2] = xx[2] - mat[6] * xx[0] - mat[7] * xx[1];

  xx[2] =  xx[2] / mat[8];
  xx[1] = (xx[1] - mat[5] * xx[2]) / mat[4];
  xx[0] = (xx[0] - mat[2] * xx[2]  - mat[1] * xx[1]) / mat[0];

  return xx;
}

void Matrix_double_3by3_Array::LU_solve(const double &b1, const double &b2, const double &b3,
            double &x1, double &x2, double &x3) const
{
  const double bb[3] = {b1, b2, b3};
  
  x1 = bb[pp[0]];
  x2 = bb[pp[1]];
  x3 = bb[pp[2]];
  
  x2 = x2 - mat[3] * x1;
  x3 = x3 - mat[6] * x1 - mat[7] * x2;

  x3 = x3 / mat[8];
  x2 = (x2 - mat[5] * x3) / mat[4];
  x1 = (x1 - mat[2] * x3 - mat[1] * x2) / mat[0];
}

void Matrix_double_3by3_Array::transpose()
{
  double temp;
  temp = mat[1]; mat[1] = mat[3]; mat[3] = temp;
  temp = mat[2]; mat[2] = mat[6]; mat[6] = temp;
  temp = mat[5]; mat[5] = mat[7]; mat[7] = temp; 
}

void Matrix_double_3by3_Array::inverse()
{
  const double invdetA = 1.0 / det();

  double temp[9];
  
  temp[0] = invdetA * (mat[4] * mat[8] - mat[5] * mat[7]);
  temp[1] = invdetA * (mat[2] * mat[7] - mat[1] * mat[8]);
  temp[2] = invdetA * (mat[1] * mat[5] - mat[2] * mat[4]);
  temp[3] = invdetA * (mat[5] * mat[6] - mat[3] * mat[8]);
  temp[4] = invdetA * (mat[0] * mat[8] - mat[2] * mat[6]);
  temp[5] = invdetA * (mat[2] * mat[3] - mat[0] * mat[5]);
  temp[6] = invdetA * (mat[3] * mat[7] - mat[4] * mat[6]);
  temp[7] = invdetA * (mat[1] * mat[6] - mat[0] * mat[7]);
  temp[8] = invdetA * (mat[0] * mat[4] - mat[1] * mat[3]);

  for(int ii=0; ii<9; ++ii) mat[ii] = temp[ii];
}

double Matrix_double_3by3_Array::det() const
{
  return mat[0] * mat[4] * mat[8] + mat[1] * mat[5] * mat[6] 
    + mat[2] * mat[3] * mat[7] - mat[2] * mat[4] * mat[6] 
    - mat[0] * mat[5] * mat[7] - mat[1] * mat[3] * mat[8];
}

void Matrix_double_3by3_Array::VecMult(const double * const &xx, double * const &yy) const
{
  yy[0] = mat[0] * xx[0] + mat[1] * xx[1] + mat[2] * xx[2];
  yy[1] = mat[3] * xx[0] + mat[4] * xx[1] + mat[5] * xx[2];
  yy[2] = mat[6] * xx[0] + mat[7] * xx[1] + mat[8] * xx[2];
}

void Matrix_double_3by3_Array::MatMult( const Matrix_double_3by3_Array &mleft,
    const Matrix_double_3by3_Array &mright )
{
  mat[0] = mleft(0) * mright(0) + mleft(1) * mright(3) + mleft(2) * mright(6);
  mat[1] = mleft(0) * mright(1) + mleft(1) * mright(4) + mleft(2) * mright(7);
  mat[2] = mleft(0) * mright(2) + mleft(1) * mright(5) + mleft(2) * mright(8);
  
  mat[3] = mleft(3) * mright(0) + mleft(4) * mright(3) + mleft(5) * mright(6);
  mat[4] = mleft(3) * mright(1) + mleft(4) * mright(4) + mleft(5) * mright(7);
  mat[5] = mleft(3) * mright(2) + mleft(4) * mright(5) + mleft(5) * mright(8);
  
  mat[6] = mleft(6) * mright(0) + mleft(7) * mright(3) + mleft(8) * mright(6);
  mat[7] = mleft(6) * mright(1) + mleft(7) * mright(4) + mleft(8) * mright(7);
  mat[8] = mleft(6) * mright(2) + mleft(7) * mright(5) + mleft(8) * mright(8);
}

void Matrix_double_3by3_Array::print() const
{
  std::cout<<std::setprecision(9)<<mat[0]<<'\t'<<mat[1]<<'\t'<<mat[2]<<std::endl;
  std::cout<<std::setprecision(9)<<mat[3]<<'\t'<<mat[4]<<'\t'<<mat[5]<<std::endl;
  std::cout<<std::setprecision(9)<<mat[6]<<'\t'<<mat[7]<<'\t'<<mat[8]<<std::endl;
}

void Matrix_double_3by3_Array::print_full() const
{
  std::cout<<"Matrix: \n";
  print();
  std::cout<<"pivoting flag: \t";
  std::cout<<std::setprecision(9)<<pp[0]<<'\t'<<pp[1]<<'\t'<<pp[2]<<'\n';
}

// EOF
