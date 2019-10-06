#include "NodalBC_2D_TNSK.hpp"


NodalBC_2D_TNSK::NodalBC_2D_TNSK(const int &nFunc, const int &nFunc_x, 
    const int &nFunc_y, const int &bc_type)
{
  clock_t log_time = clock();
  
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_dir_nodes = 0;
  num_per_nodes = 0;

  switch( bc_type )
  {
    case 1:
      BC_type_1(nFunc_x, nFunc_y);
      break;
    case 2:
      BC_type_2(nFunc_x, nFunc_y);
      break;
    case 3:
      BC_type_3();
      break;
    case 4:
      BC_type_4(nFunc_x, nFunc_y);
      break;
    case 5:
      BC_type_5(nFunc_x, nFunc_y);
      break;
    case 6:
      BC_type_6(nFunc_x, nFunc_y);
      break;
    case 7:
      BC_type_7(nFunc_x, nFunc_y);
      break;
    case 8:
      BC_type_8(nFunc_x, nFunc_y);
      break;
    case 9:
      BC_type_9(nFunc_x, nFunc_y);
      break;
    default:
      std::cerr<<"Error : NodalBC_2D_TNSK with bc_type = "<<bc_type<<" is not implemented. \n";
      exit(EXIT_FAILURE);
  }

  Create_ID(nFunc);

  log_time = clock() - log_time;
  std::cout<<"=== Nodal BC for 2D TNSK, type "<<bc_type<<" generated. ";
  std::cout<<"Time taken: "<<((float) log_time)/CLOCKS_PER_SEC<<" seconds"<<std::endl;
}



NodalBC_2D_TNSK::~NodalBC_2D_TNSK()
{}



void NodalBC_2D_TNSK::BC_type_1( const int &nFunc_x, const int &nFunc_y )
{
  std::vector<unsigned int> lef, rig, top, bot, cor;
  Generate_BCNode_2D_A(nFunc_x, nFunc_y, lef, rig, top, bot, cor);
  
  dir_nodes.insert(dir_nodes.end(), cor.begin(), cor.end());
  dir_nodes.insert(dir_nodes.end(), bot.begin(), bot.end());
  dir_nodes.insert(dir_nodes.end(), lef.begin(), lef.end());
  dir_nodes.insert(dir_nodes.end(), rig.begin(), rig.end());
  dir_nodes.insert(dir_nodes.end(), top.begin(), top.end());

  num_dir_nodes = dir_nodes.size();

  num_per_nodes = 0;

  std::cout<<"-----> C0 Dirichlet BC. \n";
}


void NodalBC_2D_TNSK::BC_type_2( const int &nFunc_x, const int &nFunc_y )
{
  std::vector<unsigned int> lef1, lef2, rig1, rig2, top1, top2, bot1, bot2;
  std::vector<unsigned int> cor1, cor2, cor3, cor4;

  Generate_BCNode_2D_B(nFunc_x, nFunc_y, lef1, rig1, top1, bot1, lef2,
      rig2, top2, bot2, cor1, cor2, cor3, cor4);

  per_slave_nodes.insert(per_slave_nodes.end(), lef1.begin(), lef1.end());
  per_master_nodes.insert(per_master_nodes.end(), lef2.begin(), lef2.end());
  per_slave_nodes.insert(per_slave_nodes.end(), rig1.begin(), rig1.end());
  per_master_nodes.insert(per_master_nodes.end(), rig2.begin(), rig2.end());
  per_slave_nodes.insert(per_slave_nodes.end(), bot1.begin(), bot1.end());
  per_master_nodes.insert(per_master_nodes.end(), bot2.begin(), bot2.end());
  per_slave_nodes.insert(per_slave_nodes.end(), top1.begin(), top1.end());
  per_master_nodes.insert(per_master_nodes.end(), top2.begin(), top2.end());

  per_master_nodes.push_back(cor4[0]);
  per_master_nodes.push_back(cor4[0]);
  per_master_nodes.push_back(cor4[0]);
  per_slave_nodes.push_back(cor1[0]);
  per_slave_nodes.push_back(cor2[0]);
  per_slave_nodes.push_back(cor3[0]);

  per_master_nodes.push_back(cor4[1]);
  per_master_nodes.push_back(cor4[1]);
  per_master_nodes.push_back(cor4[1]);
  per_slave_nodes.push_back(cor1[1]);
  per_slave_nodes.push_back(cor2[1]);
  per_slave_nodes.push_back(cor3[1]);

  per_master_nodes.push_back(cor4[2]);
  per_master_nodes.push_back(cor4[2]);
  per_master_nodes.push_back(cor4[2]);
  per_slave_nodes.push_back(cor1[2]);
  per_slave_nodes.push_back(cor2[2]);
  per_slave_nodes.push_back(cor3[2]);

  per_master_nodes.push_back(cor4[3]);
  per_master_nodes.push_back(cor4[3]);
  per_master_nodes.push_back(cor4[3]);
  per_slave_nodes.push_back(cor1[3]);
  per_slave_nodes.push_back(cor2[3]);
  per_slave_nodes.push_back(cor3[3]);

  num_per_nodes = per_slave_nodes.size();

  std::cout<<"-----> C1 grad c dot n = 0 \n";
}


void NodalBC_2D_TNSK::BC_type_3()
{
  std::cout<<"-----> Do-nothing on nodal BC. \n";
}



void NodalBC_2D_TNSK::BC_type_4( const int &nFunc_x, const int &nFunc_y )
{
  std::vector<unsigned int> lef, rig, top, bot, cor;

  Generate_BCNode_2D_A(nFunc_x, nFunc_y, lef, rig, top, bot, cor);

  dir_nodes.insert(dir_nodes.end(), cor.begin(), cor.end());
  dir_nodes.insert(dir_nodes.end(), bot.begin(), bot.end());
  dir_nodes.insert(dir_nodes.end(), top.begin(), top.end());

  num_dir_nodes = dir_nodes.size();

  num_per_nodes = 0;

  std::cout<<"-----> C0 Dirichlet BC on Top, Bottom, Corner; ";
  std::cout<<"No BC on Left / Right. \n";
}


void NodalBC_2D_TNSK::BC_type_5( const int &nFunc_x, const int &nFunc_y )
{
  std::vector<unsigned int> lef, rig, top, bot, cor;

  Generate_BCNode_2D_A(nFunc_x, nFunc_y, lef, rig, top, bot, cor);

  dir_nodes.insert(dir_nodes.end(), cor.begin(), cor.end());
  dir_nodes.insert(dir_nodes.end(), lef.begin(), lef.end());
  dir_nodes.insert(dir_nodes.end(), rig.begin(), rig.end());

  num_dir_nodes = dir_nodes.size();

  num_per_nodes = 0;

  std::cout<<"-----> C0 Dirichlet BC on Left, Right, and Corner; ";
  std::cout<<"No BC on Top / Bottom.\n";
}


void NodalBC_2D_TNSK::BC_type_6( const int &nFunc_x, const int &nFunc_y )
{
  std::vector<unsigned int> lef, rig, top, bot, cor;
  Generate_BCNode_2D_A(nFunc_x, nFunc_y, lef, rig, top, bot, cor);

  per_master_nodes.insert(per_master_nodes.end(), bot.begin(), bot.end());
  per_slave_nodes.insert(per_slave_nodes.end(), top.begin(), top.end());

  per_master_nodes.insert(per_master_nodes.end(), lef.begin(), lef.end());
  per_slave_nodes.insert(per_slave_nodes.end(), rig.begin(), rig.end());

  per_master_nodes.push_back(cor[0]);
  per_master_nodes.push_back(cor[0]);
  per_master_nodes.push_back(cor[0]);

  per_slave_nodes.push_back(cor[1]);
  per_slave_nodes.push_back(cor[2]);
  per_slave_nodes.push_back(cor[3]);

  num_per_nodes = per_slave_nodes.size();

  std::cout<<"-----> C0 periodic strong nodal imposition. \n";
}



void NodalBC_2D_TNSK::BC_type_7( const int &nFunc_x, const int &nFunc_y )
{
  std::vector<unsigned int> lef1, lef2, rig1, rig2, top1, top2, bot1, bot2;
  std::vector<unsigned int> cor1, cor2, cor3, cor4;

  Generate_BCNode_2D_B(nFunc_x, nFunc_y, lef1, rig1, top1, bot1, lef2,
      rig2, top2, bot2, cor1, cor2, cor3, cor4);

  VEC_T::insert_end( per_master_nodes, lef1 );
  VEC_T::insert_end( per_slave_nodes, rig2 );

  VEC_T::insert_end( per_master_nodes, lef2 );
  VEC_T::insert_end( per_slave_nodes, rig1 );

  VEC_T::insert_end( per_master_nodes, bot2 );
  VEC_T::insert_end( per_slave_nodes, bot1 );
  
  VEC_T::insert_end( per_master_nodes, top2 );
  VEC_T::insert_end( per_slave_nodes, top1 );

  per_master_nodes.push_back(cor4[0]); 
  per_master_nodes.push_back(cor4[0]); 
  per_master_nodes.push_back(cor4[0]); 
  per_slave_nodes.push_back(cor2[0]);
  per_slave_nodes.push_back(cor3[1]);
  per_slave_nodes.push_back(cor1[1]);

  
  per_master_nodes.push_back(cor3[0]);
  per_master_nodes.push_back(cor3[0]);
  per_master_nodes.push_back(cor3[0]);
  per_slave_nodes.push_back(cor1[0]);
  per_slave_nodes.push_back(cor4[1]);
  per_slave_nodes.push_back(cor2[1]);

  
  per_master_nodes.push_back(cor4[2]);
  per_master_nodes.push_back(cor4[2]);
  per_master_nodes.push_back(cor4[2]);
  per_slave_nodes.push_back(cor2[2]);
  per_slave_nodes.push_back(cor3[3]);
  per_slave_nodes.push_back(cor1[3]);


  per_master_nodes.push_back(cor3[2]);
  per_master_nodes.push_back(cor3[2]);
  per_master_nodes.push_back(cor3[2]);
  per_slave_nodes.push_back(cor1[2]);
  per_slave_nodes.push_back(cor4[3]);
  per_slave_nodes.push_back(cor2[3]);

  dir_nodes.clear();
  num_per_nodes = per_slave_nodes.size();
  num_dir_nodes = 0;

  std::cout<<"-----> C1 periodic in s-direction (left - right face) ";
  std::cout<<"grad c dot n = 0 on top and bottom. \n";
  std::cout<<"       Note: This BC only works for periodic splines. See JCP 242:321-350. \n";
}



void NodalBC_2D_TNSK::BC_type_8( const int &nFunc_x, const int &nFunc_y )
{
  std::vector<unsigned int> lef1, lef2, rig1, rig2, top1, top2, bot1, bot2;
  std::vector<unsigned int> cor1, cor2, cor3, cor4;

  Generate_BCNode_2D_B(nFunc_x, nFunc_y, lef1, rig1, top1, bot1, lef2,
      rig2, top2, bot2, cor1, cor2, cor3, cor4);

  VEC_T::insert_end( per_master_nodes, lef1 );
  VEC_T::insert_end( per_slave_nodes, rig2 );

  VEC_T::insert_end( per_master_nodes, lef2 );
  VEC_T::insert_end( per_slave_nodes, rig1 );

  per_master_nodes.push_back(cor4[0]);
  per_slave_nodes.push_back(cor3[1]);

  per_master_nodes.push_back(cor3[0]);
  per_slave_nodes.push_back(cor4[1]);

  per_master_nodes.push_back(cor4[2]);
  per_slave_nodes.push_back(cor3[3]);
  
  per_master_nodes.push_back(cor3[2]);
  per_slave_nodes.push_back(cor4[3]);

  VEC_T::insert_end( dir_nodes, bot1 );
  VEC_T::insert_end( dir_nodes, top1 );

  dir_nodes.push_back(cor1[0]);
  dir_nodes.push_back(cor2[0]);
  dir_nodes.push_back(cor2[1]);
  dir_nodes.push_back(cor1[1]);
  
  dir_nodes.push_back(cor1[2]);
  dir_nodes.push_back(cor2[2]);
  dir_nodes.push_back(cor2[3]);
  dir_nodes.push_back(cor1[3]);

  num_per_nodes = per_slave_nodes.size();
  num_dir_nodes = dir_nodes.size();

  std::cout<<"-----> C1 periodic in s-direction (left - right face) ";
  std::cout<<"Dirichlet on top and bottom. \n";
  std::cout<<"       Note: This BC only works for periodic splines. See JCP 242: 321-350. \n";
}



void NodalBC_2D_TNSK::BC_type_9( const int &nFunc_x, const int &nFunc_y )
{
  std::vector<unsigned int> lef1, lef2, rig1, rig2, top1, top2, bot1, bot2;
  std::vector<unsigned int> cor1, cor2, cor3, cor4;

  Generate_BCNode_2D_B(nFunc_x, nFunc_y, lef1, rig1, top1, bot1, lef2,
      rig2, top2, bot2, cor1, cor2, cor3, cor4);


  VEC_T::insert_end( per_master_nodes, lef1 );
  VEC_T::insert_end( per_slave_nodes, rig2 );
  
  VEC_T::insert_end( per_master_nodes, lef2 );
  VEC_T::insert_end( per_slave_nodes, rig1 );

  per_master_nodes.push_back(cor4[0]); // 6
  per_slave_nodes.push_back(cor3[1]);  // 9

  per_master_nodes.push_back(cor3[0]); // 5
  per_slave_nodes.push_back(cor4[1]);  // 8

  per_master_nodes.push_back(cor4[2]); // 16
  per_slave_nodes.push_back(cor3[3]);  // 19
  
  per_master_nodes.push_back(cor3[2]); // 15
  per_slave_nodes.push_back(cor4[3]);  // 18

  per_master_nodes.push_back(cor1[0]); // 0
  per_slave_nodes.push_back(cor2[1]);  // 3

  per_master_nodes.push_back(cor2[0]); // 1
  per_slave_nodes.push_back(cor1[1]);  // 4

  per_master_nodes.push_back(cor2[2]); // 21
  per_slave_nodes.push_back(cor1[3]);  // 24

  per_master_nodes.push_back(cor1[2]); // 20
  per_slave_nodes.push_back(cor2[3]);  // 23


  num_per_nodes = per_slave_nodes.size();
  num_dir_nodes = 0;

  std::cout<<"-----> C1 periodic in s-direction (left-right face). ";
  std::cout<<"Nothing on top and bottom. \n";
  std::cout<<"       Note: This BC only works for periodic splines. See JCP 242: 321-350. \n";
}










// EOF
