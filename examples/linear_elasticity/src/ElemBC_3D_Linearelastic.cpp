#include "ElemBC_3D_Linearelastic.hpp"

ElemBC_3D_Linearelastic::ElemBC_3D_Linearelastic( const int &nElem_x, 
    const int &nElem_y, const int &nElem_z, const int &bc_type )
{
  clock_t log_time = clock();

  left_elem.clear(); right_elem.clear();
  top_elem.clear(); bottom_elem.clear();
  front_elem.clear(); back_elem.clear();

  switch(bc_type)
  {
    case 0:
      BC_test( nElem_x, nElem_y, nElem_z );
      break;
    case 1:
      BC_type_1();
      break;
    case 2:
      BC_type_2(nElem_x, nElem_y, nElem_z);
      break;
    case 3:
      BC_type_3(nElem_x, nElem_y, nElem_z);
      break;
    case 4:
      BC_type_4(nElem_x, nElem_y, nElem_z);
      break;
    case 5:
      BC_type_5(nElem_x, nElem_y, nElem_z);
      break;
    case 6:
      BC_type_6(nElem_x, nElem_y, nElem_z);
      break;
    case 7:
      BC_type_7(nElem_x, nElem_y, nElem_z);
      break;
    default:
      std::cerr<<"Error: ElemBC_3D_Linearelastic with type "
        <<bc_type<<" is not implemented. \n";
      exit(EXIT_FAILURE);
  }
  
  log_time = clock() - log_time;
  std::cout<<"=== Elem BC for 3D Hyperelastic with type "<<bc_type<<" generated. ";
  std::cout<<"Time taken: "<<((float) log_time)/CLOCKS_PER_SEC<<" seconds. \n";
}



ElemBC_3D_Linearelastic::~ElemBC_3D_Linearelastic()
{}



void ElemBC_3D_Linearelastic::BC_test( const int &nElem_x, 
    const int &nElem_y, const int &nElem_z)
{
  Generate_BCElem_3D_A( nElem_x, nElem_y, nElem_z,
      front_elem, back_elem, left_elem, right_elem, top_elem, bottom_elem );

  //left_elem.clear(); 
  right_elem.clear();
  top_elem.clear();  
  bottom_elem.clear();
  front_elem.clear(); 
  back_elem.clear();
  
  std::cout<<"-----> ElemBC_3D_Linearelastic: DEBUGGING MODE. \n";
}


void ElemBC_3D_Linearelastic::BC_type_1()
{
  std::cout<<"-----> No boundary surface intergal. \n";
}



void ElemBC_3D_Linearelastic::BC_type_2( const int &nElem_x, 
    const int &nElem_y, const int &nElem_z)
{
  Generate_BCElem_3D_A( nElem_x, nElem_y, nElem_z,
      front_elem, back_elem, left_elem, right_elem, top_elem, bottom_elem );

  std::cout<<"-----> All six boundary faces require surface intergal. \n";
}



void ElemBC_3D_Linearelastic::BC_type_3( const int &nElem_x, 
    const int &nElem_y, const int &nElem_z)
{
  Generate_BCElem_3D_A( nElem_x, nElem_y, nElem_z,
      front_elem, back_elem, left_elem, right_elem, top_elem, bottom_elem );

  left_elem.clear(); right_elem.clear();
  bottom_elem.clear();
  front_elem.clear(); back_elem.clear();
  
  std::cout<<"-----> TOP boundary faces require surface intergal. \n";
}



void ElemBC_3D_Linearelastic::BC_type_4( const int &nElem_x, 
    const int &nElem_y, const int &nElem_z)
{
  Generate_BCElem_3D_A( nElem_x, nElem_y, nElem_z,
      front_elem, back_elem, left_elem, right_elem, top_elem, bottom_elem );

  left_elem.clear(); right_elem.clear();
  top_elem.clear();
  front_elem.clear(); back_elem.clear();
  
  std::cout<<"-----> BOTTOM boundary faces require surface intergal. \n";
}


void ElemBC_3D_Linearelastic::BC_type_5( const int &nElem_x, 
    const int &nElem_y, const int &nElem_z)
{
  Generate_BCElem_3D_A( nElem_x, nElem_y, nElem_z,
      front_elem, back_elem, left_elem, right_elem, top_elem, bottom_elem );

  back_elem.clear(); 
  
  std::cout<<"-----> TOP & BOTTOM, FRONT, LEFT & RIGHT boundary faces require surface intergal. \n";
}


void ElemBC_3D_Linearelastic::BC_type_6( const int &nElem_x, 
    const int &nElem_y, const int &nElem_z)
{
  Generate_BCElem_3D_A( nElem_x, nElem_y, nElem_z,
      front_elem, back_elem, left_elem, right_elem, top_elem, bottom_elem );

  back_elem.clear(); 
  front_elem.clear();

  std::cout<<"-----> TOP & BOTTOM, LEFT & RIGHT boundary faces require surface intergal. \n";
}


void ElemBC_3D_Linearelastic::BC_type_7( const int &nElem_x, 
    const int &nElem_y, const int &nElem_z)
{
  Generate_BCElem_3D_A( nElem_x, nElem_y, nElem_z,
      front_elem, back_elem, left_elem, right_elem, top_elem, bottom_elem );

  bottom_elem.clear();

  std::cout<<"-----> TOP, LEFT & RIGHT, FRONT & BACK boundary faces require surface intergal. \n";
}


// EOF
