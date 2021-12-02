#include "NodalBC_3D_FSI.hpp"

NodalBC_3D_FSI::NodalBC_3D_FSI( const std::string &fluid_file,
    const std::string &solid_file,
    const std::string &fluid_wall_file,
    const std::string &solid_wall_file,
    const std::vector<std::string> &fluid_inlet_files,
    const std::vector<std::string> &fluid_outlet_files,
    const std::vector<std::string> &solid_inlet_files,
    const std::vector<std::string> &solid_outlet_files,
    const int &nFunc,
    const int &comp,
    const int &ringBC_type,
    const int &fsiBC_type )
{
  dir_nodes.clear();

  switch( fsiBC_type )
  {
    // ====== deformable FSI ======
    case 0:
      {
        if(ringBC_type == 0)
        {}
        else if(ringBC_type == 1)
        {}
        else SYS_T::print_fatal( "Error: NodalBC_3D_FSI has no such type of essential bc for ring nodes.\n" );

        break;
      }
    
      // ====== Rigid wall ======
    case 1:
      {
        if( comp == 0 ) dir_nodes.clear();
        else if( comp == 1 || comp == 2 || comp == 3 )
        {
          dir_nodes = TET_T::read_int_PointData( solid_file, "GlobalNodeID" );

          VEC_T::insert_end( dir_nodes, get_vtp_nodal_id( fluid_inlet_files ) );
        
          VEC_T::sort_unique_resize( dir_nodes );
        }
        else SYS_T::print_fatal( "Error: NodalBC_3D_FSI has no such type of component index.\n" );

        break;
      }

      // ====== Wall prestress solver =====
    case 2:
      {
        if(ringBC_type == 0)
        {
          if( comp == 0 )
            dir_nodes = TET_T::read_int_PointData( fluid_file, "GlobalNodeID" );
          else if( comp == 1 || comp == 2 || comp == 3 )
          {
          }
          else SYS_T::print_fatal( "Error: NodalBC_3D_FSI has no such type of component index.\n" );

        }
        else if(ringBC_type == 1)
        {
        
        }
        else SYS_T::print_fatal( "Error: NodalBC_3D_FSI has no such type of essential bc for ring nodes.\n" );

        break;
      }

    default:
      SYS_T::print_fatal( "NodalBC_3D_FSI Error: No such type of FSI BC.\n" );
      break;
  }

  // count the number of dirichlet nodes
  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  // Print info
}

std::vector<int> NodalBC_3D_FSI::get_vtp_nodal_id( const std::vector<std::string> &vtpfileList ) const
{
  std::vector<int> output; output.clear();

  const unsigned int num_file = vtpfileList.size();

  for(unsigned int ii=0; ii<num_file; ++ii)
  {
    SYS_T::file_check( vtpfileList[ii] );

    VEC_T::insert_end( output, TET_T::read_int_PointData( vtpfileList[ii], "GlobalNodeID" ) );
  }

  return output;
}

]// EOF
