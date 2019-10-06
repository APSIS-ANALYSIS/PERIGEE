#include "ElemBC_3D_TNSK.hpp"


ElemBC_3D_TNSK::ElemBC_3D_TNSK( const int &nElem_x, const int &nElem_y, const int &nElem_z,
            const int &bc_type )
{
  clock_t log_time = clock();
  
  left_elem.clear(); right_elem.clear();
  top_elem.clear(); bottom_elem.clear();
  front_elem.clear(); back_elem.clear();

  switch(bc_type)
  {
    case 1:
      BC_type_1();
      break;
    case 2:
      BC_type_2(nElem_x, nElem_y, nElem_z);
      break;
    default:
      std::cerr<<"Error: ElemBC_3D_TNSK with type "<<bc_type<<" is not implemented. \n";
      exit(EXIT_FAILURE);
  }

  log_time = clock() - log_time;
  std::cout<<"=== Elem BC for 3D TNSK with type "<<bc_type<<" generated. ";
  std::cout<<"Time taken: "<<((float) log_time)/CLOCKS_PER_SEC<<" seconds. \n";
}



ElemBC_3D_TNSK::~ElemBC_3D_TNSK()
{}



void ElemBC_3D_TNSK::BC_type_1()
{
  std::cout<<"-----> Do-nothing on Element boundary faces. \n";
}



void ElemBC_3D_TNSK::BC_type_2(const int &nElem_x, const int &nElem_y, const int &nElem_z)
{
  Generate_BCElem_3D_A( nElem_x, nElem_y, nElem_z,
      front_elem, back_elem, left_elem, right_elem, top_elem, bottom_elem );

  std::cout<<"-----> Element boundary faces are specified for left, right, bottom top, front, back. \n";
}


// EOF
