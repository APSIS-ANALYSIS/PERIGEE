#include "Matrix_double_6by6_Array.hpp"

Matrix_double_6by6_Array::Matrix_double_6by6_Array(const double * const &in_array)
{
  int counter = 0;
  for(int ii=0; ii<6; ++ii)
  {
    for(int jj=0; jj<6; ++jj)
    {
      Mat[ii][jj] = in_array[counter];
      counter += 1;
    }
    p[ii] = ii;
  }
}

Matrix_double_6by6_Array::Matrix_double_6by6_Array(const double &a, const double &b,
    const double &c, const double &d, const double &e, const double &f,
    const double &g, const double &h, const double &i)
{
  Mat[0][0] = a*a;   Mat[0][1] = d*d;   Mat[0][2] = g*g;
  Mat[0][3] = 2.0*a*d; Mat[0][4] = 2.0*a*g; Mat[0][5] = 2.0*d*g;

  Mat[1][0] = a*b; Mat[1][1] = d*e; Mat[1][2] = g*h;
  Mat[1][3] = a*e+b*d; Mat[1][4] = a*h+b*g; Mat[1][5] = d*h+e*g;
  
  Mat[2][0] = a*c; Mat[2][1] = d*f; Mat[2][2] = g*i;
  Mat[2][3] = a*f+c*d; Mat[2][4] = a*i+c*g; Mat[2][5] = d*i+f*g;
  
  Mat[3][0] = b*b;   Mat[3][1] = e*e;   Mat[3][2] = h*h;
  Mat[3][3] = 2.0*b*e; Mat[3][4] = 2.0*b*h; Mat[3][5] = 2.0*e*h;
  
  Mat[4][0] = b*c; Mat[4][1] = e*f; Mat[4][2] = h*i;
  Mat[4][3] = b*f+c*e; Mat[4][4] = b*i+c*h; Mat[4][5] = e*i+f*h;

  Mat[5][0] = c*c; Mat[5][1] = f*f; Mat[5][2] = i*i;
  Mat[5][3] = 2.0*c*f; Mat[5][4] = 2.0*c*i; Mat[5][5] = 2.0*f*i;

  p[0] = 0; p[1] = 1; p[2] = 2; p[3] = 3; p[4] = 4; p[5] = 5;
}

Matrix_double_6by6_Array::~Matrix_double_6by6_Array()
{}

void Matrix_double_6by6_Array::LU_fac()
{
  double max_value, temp, invAkk;
  int max_index, int_temp;
  bool pivot_flag;
  for(int kk=0; kk<5; ++kk)
  {
    max_value = std::abs(Mat[kk][kk]);
    max_index = kk;
    pivot_flag = false;
    // find the column pivoting
    for(int ii=kk+1; ii<6; ++ii)
    {
      if(max_value < std::abs(Mat[ii][kk]))
      {
        max_value = std::abs(Mat[ii][kk]);
        max_index = ii;
        pivot_flag = true;
      }
    }
    
    if(pivot_flag)
    {
      int_temp = p[kk];
      p[kk] = p[max_index];
      p[max_index] = int_temp;
      
      for(int ii=0; ii<6; ++ii)
      {
        temp = Mat[kk][ii];
        Mat[kk][ii] = Mat[max_index][ii];
        Mat[max_index][ii] = temp;
      }
    }

    invAkk = 1.0 / Mat[kk][kk];
    
    for(int ii=kk+1; ii<6; ++ii)
    {
      Mat[ii][kk] = Mat[ii][kk] * invAkk;
      for(int jj=kk+1; jj<6; ++jj)
        Mat[ii][jj] = Mat[ii][jj] - Mat[ii][kk] * Mat[kk][jj];
    }
  }
  invm0 = 1.0 / Mat[0][0]; invm1 = 1.0 / Mat[1][1]; invm2 = 1.0 / Mat[2][2];
  invm3 = 1.0 / Mat[3][3]; invm4 = 1.0 / Mat[4][4]; invm5 = 1.0 / Mat[5][5];
}

void Matrix_double_6by6_Array::LU_solve(const double * const &b, double * const &x) const
{
  for(int ii=0; ii<6; ++ii)
    x[ii] = b[p[ii]];

  //x[0] = x[0];
  x[1] = x[1] - Mat[1][0] * x[0];
  x[2] = x[2] - Mat[2][0] * x[0] - Mat[2][1] * x[1];
  x[3] = x[3] - Mat[3][0] * x[0] - Mat[3][1] * x[1] - Mat[3][2] * x[2];
  x[4] = x[4] - Mat[4][0] * x[0] - Mat[4][1] * x[1] - Mat[4][2] * x[2] - Mat[4][3] * x[3];
  x[5] = x[5] - Mat[5][0] * x[0] - Mat[5][1] * x[1] - Mat[5][2] * x[2] - Mat[5][3] * x[3] - Mat[5][4]* x[4];

  x[5] = x[5] * invm5;
  x[4] = (x[4] - Mat[4][5]*x[5]) * invm4;
  x[3] = (x[3] - Mat[3][5]*x[5] - Mat[3][4] * x[4]) * invm3;
  x[2] = (x[2] - Mat[2][5]*x[5] - Mat[2][4]*x[4] - Mat[2][3] * x[3] ) * invm2;
  x[1] = (x[1] - Mat[1][5]*x[5] - Mat[1][4]*x[4] - Mat[1][3]*x[3] - Mat[1][2]*x[2]) * invm1;
  x[0] = (x[0] - Mat[0][5]*x[5] - Mat[0][4]*x[4] - Mat[0][3]*x[3] - Mat[0][2]*x[2] -Mat[0][1]*x[1]) * invm0;
}


void Matrix_double_6by6_Array::print() const
{
  for(int ii=0; ii<6; ++ii)
  {
    for(int jj=0; jj<6; ++jj)
      std::cout<<Mat[ii][jj]<<'\t';
    std::cout<<'\n';
  }
  std::cout<<'\n';
  for(int ii=0; ii<6; ++ii)
    std::cout<<p[ii]<<'\t';
  
  std::cout<<'\n'<<std::endl;
}

// EOF
