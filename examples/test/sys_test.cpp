#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"
#include "IEN_FEM.hpp"
#include "NodalBC.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D.hpp"
#include "ElemBC_3D_outflow.hpp"

#include <engine.h>
#include <matrix.h>

int main(int argc, char *argv[])
{
  Engine *ep;
  ep = engOpen(nullptr);

  Tensor2_3D ten; ten.gen_rand();
  Tensor2_3D exp1; exp1 = Ten2::exp(ten);
  const int arraySize = 9;
  double randomNumbers[arraySize] = { ten(0), ten(1), ten(2),
                                      ten(3), ten(4), ten(5),
                                      ten(6), ten(7), ten(8) };

  mxArray *matlabArray = mxCreateDoubleMatrix(1, arraySize, mxREAL);
  double *matlabArrayData = mxGetPr(matlabArray);

  for (int ii=0; ii<arraySize; ii++) matlabArrayData[ii] = randomNumbers[ii];

  engPutVariable(ep, "inputArray", matlabArray);
  const char *matlabFunction = "result = expm(inputArray);";
  engEvalString(ep, matlabFunction);

  mxArray *result = engGetVariable(ep, "result");
  double *resultData = mxGetPr(result);

  Tensor2_3D exp2 = Ten2::gen_zero();

  for(mwIndex ii=0; ii<9; ii++) exp2(ii) = resultData[ii];

  mxDestroyArray(result);
  mxDestroyArray(matlabArray);

  engClose(ep);

  if(exp1.is_identical(exp2)) std::cout<<"SUCESS"<<std::endl;
  else std::cout<<"WRONG!"<<std::endl;

  engClose(ep);
  return EXIT_SUCCESS;
}

// EOF
