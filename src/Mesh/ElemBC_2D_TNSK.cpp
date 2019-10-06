#include "ElemBC_2D_TNSK.hpp"


ElemBC_2D_TNSK::ElemBC_2D_TNSK(const int &nElem_x, const int &nElem_y,
            const int &bc_type)
{
  clock_t log_time = clock();

  switch(bc_type)
  {
    case 1:
      BC_type_1();
      break;
    case 2:
      BC_type_2(nElem_x, nElem_y);
      break;
    default:
      std::cerr<<"Error: ElemBC_2D_TNSK with type "<<bc_type<<" is not implemented. \n";
      exit(EXIT_FAILURE);
  }

  log_time = clock() - log_time;
  std::cout<<"=== Elem BC for 2D TNSK with type "<<bc_type<<" generated. ";
  std::cout<<"Time taken: "<<((float) log_time)/CLOCKS_PER_SEC<<" seconds. \n";
}



ElemBC_2D_TNSK::~ElemBC_2D_TNSK()
{}



void ElemBC_2D_TNSK::BC_type_1()
{
  std::cout<<"-----> Do-nothing on Element boundary faces. \n";
}



void ElemBC_2D_TNSK::BC_type_2(const int &nElem_x, const int &nElem_y)
{
  Generate_BCElem_2D(nElem_x, nElem_y, front_elem, back_elem, left_elem, right_elem);
  std::cout<<"-----> Element boundary faces are specified for front, back, left, and right. \n";
}


// EOF
