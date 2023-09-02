#include "Matrix_double_3by3_Array.hpp"

Matrix_double_3by3_Array::Matrix_double_3by3_Array()
{
  mat[0] = 1.0; mat[1] = 0.0; mat[2] = 0.0;
  mat[3] = 0.0; mat[4] = 1.0; mat[5] = 0.0;
  mat[6] = 0.0; mat[7] = 0.0; mat[8] = 1.0;

  p[0] = 0; p[1] = 1; p[2] = 2;
 
  invm0 = 1.0; invm1 = 1.0; invm2 = 1.0;
}

Matrix_double_3by3_Array::Matrix_double_3by3_Array( const double * const &in_mat )
{
  for(int ii=0; ii<9; ++ii) mat[ii] = in_mat[ii];

  p[0] = 0; p[1] = 1; p[2] = 2;
  
  invm0 = 1.0; invm1 = 1.0; invm2 = 1.0;
}

Matrix_double_3by3_Array::Matrix_double_3by3_Array( 
    const double &a11, const double &a12, const double &a13,
    const double &a21, const double &a22, const double &a23,
    const double &a31, const double &a32, const double &a33 )
{
  mat[0] = a11;  mat[1] = a12;  mat[2] = a13;
  mat[3] = a21;  mat[4] = a22;  mat[5] = a23;
  mat[6] = a31;  mat[7] = a32;  mat[8] = a33;

  p[0] = 0; p[1] = 1; p[2] = 2;
  
  invm0 = 1.0; invm1 = 1.0; invm2 = 1.0;
}

Matrix_double_3by3_Array::Matrix_double_3by3_Array( 
    const Matrix_double_3by3_Array &other )
{
  for(int ii=0; ii<9; ++ii) mat[ii] = other.mat[ii];
  
  p[0] = other.p[0]; 
  p[1] = other.p[1]; 
  p[2] = other.p[2];
  
  invm0 = other.invm0;
  invm1 = other.invm1;
  invm2 = other.invm2; 
}

Matrix_double_3by3_Array::~Matrix_double_3by3_Array()
{}

Matrix_double_3by3_Array& Matrix_double_3by3_Array::operator= (
    const Matrix_double_3by3_Array &input )
{
  if(this != &input)
  {
    for(int ii=0; ii<9; ++ii) mat[ii] = input.mat[ii];
    p[0] = input.p[0];
    p[1] = input.p[1];
    p[2] = input.p[2];
    invm0 = input.invm0;
    invm1 = input.invm1;
    invm2 = input.invm2;
  }
  return *this;
}

void Matrix_double_3by3_Array::gen_id()
{
  mat[0] = 1.0;  mat[1] = 0.0;  mat[2] = 0.0;
  mat[3] = 0.0;  mat[4] = 1.0;  mat[5] = 0.0;
  mat[6] = 0.0;  mat[7] = 0.0;  mat[8] = 1.0;

  p[0] = 0; p[1] = 1; p[2] = 2;
  
  invm0 = 1.0; invm1 = 1.0; invm2 = 1.0;
}

void Matrix_double_3by3_Array::gen_rand()
{
  std::random_device rd;
  std::mt19937_64 gen( rd() );
  std::uniform_real_distribution<double> dis(min, max);
  for(int ii=0; ii<9; ++ii) mat[ii] = dis(gen); 

  p[0] = 0; p[1] = 1; p[2] = 2;
  
  invm0 = 1.0; invm1 = 1.0; invm2 = 1.0;
}

void Matrix_double_3by3_Array::gen_hilb()
{
  for(int ii=0; ii<3; ++ii)
    for(int jj=0; jj<3; ++jj)
      mat[ii*3+jj] = 1.0 / (ii + jj + 1.0);

  p[0] = 0; p[1] = 1; p[2] = 2;
  
  invm0 = 1.0; invm1 = 1.0; invm2 = 1.0;
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
    int_temp = p[0];
    p[0] = p[max_index];
    p[max_index] = int_temp;

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

  double invA0 = 1.0 / mat[0];
  mat[3] = mat[3] * invA0;
  mat[4] = mat[4] - mat[3] * mat[1];
  mat[5] = mat[5] - mat[3] * mat[2];
  mat[6] = mat[6] * invA0;
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
    int_temp = p[1];
    p[1] = p[2];
    p[2] = int_temp;

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

  double invA1 = 1.0 / mat[4];
  mat[7] = mat[7] * invA1;
  mat[8] = mat[8] - mat[7] * mat[5];

  invm0 = 1.0 / mat[0]; invm1 = 1.0 / mat[4]; invm2 = 1.0 / mat[8];
}

Vector_3 Matrix_double_3by3_Array::LU_solve( const Vector_3 &b ) const
{
  Vector_3 x( b(p[0]), b(p[1]), b(p[2]) );

  x(1) = x(1) - mat[3] * x(0);
  x(2) = x(2) - mat[6] * x(0) - mat[7] * x(1);

  x(2) = x(2) * invm2;
  x(1) = (x(1) - mat[5] * x(2)) * invm1;
  x(0) = (x(0) - mat[2] * x(2) - mat[1] * x(1)) * invm0;

  return x;
}

std::array<double, 3> Matrix_double_3by3_Array::LU_solve( const std::array<double, 3> &b ) const
{
  std::array<double, 3> x {{ b[p[0]], b[p[1]], b[p[2]] }};
  
  x[1] = x[1] - mat[3] * x[0];
  x[2] = x[2] - mat[6] * x[0] - mat[7] * x[1];

  x[2] = x[2] * invm2;
  x[1] = (x[1] - mat[5] * x[2]) * invm1;
  x[0] = (x[0] - mat[2] * x[2] - mat[1] * x[1]) * invm0;

  return x;
}

void Matrix_double_3by3_Array::LU_solve(const double &b1, const double &b2, const double &b3,
            double &x1, double &x2, double &x3) const
{
  const double b[3] = {b1, b2, b3};
  
  x1 = b[p[0]];
  x2 = b[p[1]];
  x3 = b[p[2]];
  
  x2 = x2 - mat[3] * x1;
  x3 = x3 - mat[6] * x1 - mat[7] * x2;

  x3 = x3 * invm2;
  x2 = (x2 - mat[5] * x3) * invm1;
  x1 = (x1 - mat[2] * x3 - mat[1] * x2) * invm0;
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

void Matrix_double_3by3_Array::VecMult(const double * const &x, double * const &y) const
{
  y[0] = mat[0] * x[0] + mat[1] * x[1] + mat[2] * x[2];
  y[1] = mat[3] * x[0] + mat[4] * x[1] + mat[5] * x[2];
  y[2] = mat[6] * x[0] + mat[7] * x[1] + mat[8] * x[2];
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
  std::cout<<std::setprecision(9)<<p[0]<<'\t'<<p[1]<<'\t'<<p[2]<<'\n';
  std::cout<<"invm0 = "<<invm0<<'\t';
  std::cout<<"invm1 = "<<invm1<<'\t';
  std::cout<<"invm2 = "<<invm2<<'\n';
}

// EOF
