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
  SYS_T::print_fatal_if( comp != 0 && comp != 1 && comp != 2, "Error: NodalBC_3D_FSI comp argument should be 0, 1, or 2. \n");
  SYS_T::print_fatal_if( ringBC_type != 0 && ringBC_type != 1, "Error: NodalBC_3D_FSI has no such type of ringBC.\n");

  dir_nodes.clear();

  switch( fsiBC_type )
  {
    // ====== deformable FSI ======
    case 0:
      {
        if(ringBC_type == 0)
        {
          // VEC_T::insert_end( dir_nodes, get_vtk_nodal_id( solid_inlet_files ) );
          // VEC_T::insert_end( dir_nodes, get_vtk_nodal_id( solid_outlet_files ) );

          // For coronary artery benchmark:
          if(comp == 0) // fix x-direction disp of solid outlet
            VEC_T::insert_end( dir_nodes, get_vtk_nodal_id( solid_outlet_files ) );
          else if(comp == 1)  // fix y-direction disp of redlines
          {
            VEC_T::insert_end( dir_nodes, CA_benchmark_BC( solid_inlet_files ) );
            VEC_T::insert_end( dir_nodes, CA_benchmark_BC( solid_outlet_files ) );
          }
          else if(comp == 2) // fix z-direction disp of solid inlet
          {
            dir_nodes = get_vtk_nodal_id( fluid_inlet_files );
            VEC_T::insert_end( dir_nodes, get_vtk_nodal_id( solid_inlet_files ) );
          }
          else
            ;
          VEC_T::sort_unique_resize( dir_nodes );

          std::cout<<"===> NodalBC_3D_FSI for deformable wall (fsiBC_type = 0) with cap surface \n";
          std::cout<<"     fully clamped (ringBC_type = 0) is generated for displacement component ["<<comp<<"]. \n";
        }
        else
        {
          // TO BE FILLED
        }
        break;
      }

      // ====== Rigid wall ======
    case 1:
      {
        dir_nodes = VEC_T::cast_to_unsigned_int( VTK_T::read_int_PointData( solid_file, "GlobalNodeID" ) );

        VEC_T::insert_end( dir_nodes, get_vtk_nodal_id( fluid_inlet_files ) );

        VEC_T::sort_unique_resize( dir_nodes );

        std::cout<<"===> NodalBC_3D_FSI for rigid wall (fsiBC_type = 1) is generated for displacement component ["<<comp<<"]. \n";
        break;
      }

      // ====== Wall prestress solver =====
    case 2:
      {
        if(ringBC_type == 0)
        {
          const std::vector<int> f_node = VTK_T::read_int_PointData( fluid_file, "GlobalNodeID" );
          const std::vector<int> fwall_node = VTK_T::read_int_PointData( fluid_wall_file, "GlobalNodeID" );

          dir_nodes = VEC_T::cast_to_unsigned_int( VEC_T::set_diff( f_node, fwall_node ) );

          VEC_T::insert_end( dir_nodes, get_vtk_nodal_id( solid_inlet_files ) );
          VEC_T::insert_end( dir_nodes, get_vtk_nodal_id( solid_outlet_files ) );
          VEC_T::sort_unique_resize( dir_nodes );

          std::cout<<"===> NodalBC_3D_FSI for wall prestressing (fsiBC_type = 2) with cap surface \n";
          std::cout<<"     fully clamped (ringBC_type = 0) is generated for displacement component ["<<comp<<"]. \n";
        }
        else
        {
          // TO BE FILLED
        }
        break;
      }

    default:
      SYS_T::print_fatal( "NodalBC_3D_FSI Error: No such type of FSI BC.\n" );
      break;
  }

  // count the number of dirichlet nodes
  num_dir_nodes = dir_nodes.size();

  // generate the ID array
  Create_ID( nFunc );
}

NodalBC_3D_FSI::NodalBC_3D_FSI( const std::string &fluid_file,
    const int &nFunc,
    const int &fsiBC_type )
{
  dir_nodes.clear();

  switch( fsiBC_type )
  {
    // ====== deformable FSI ======
    case 0:
      dir_nodes.clear();
      std::cout<<"===> NodalBC_3D_FSI for deformable wall (fsiBC_type = 0) is generated for pressure.\n";
      break;
      // ====== Rigid wall ======
    case 1:
      dir_nodes.clear();
      std::cout<<"===> NodalBC_3D_FSI for rigid wall (fsiBC_type = 1) is generated for pressure.\n";
      break;
      // ====== Wall prestress solver =====
    case 2:
      dir_nodes = VEC_T::cast_to_unsigned_int( VTK_T::read_int_PointData( fluid_file, "GlobalNodeID" ) );
      std::cout<<"===> NodalBC_3D_FSI for wall prestressing (fsiBC_type = 2) is generated for pressure \n";
      break;
    default:
      SYS_T::print_fatal( "NodalBC_3D_FSI Error: No such type of FSI BC.\n" );
      break;
  }

  // count the number of dirichlet nodes
  num_dir_nodes = dir_nodes.size();

  // generate the ID array
  Create_ID( nFunc );
}

std::vector<unsigned int> NodalBC_3D_FSI::get_vtk_nodal_id( const std::vector<std::string> &vtkfileList ) const
{
  std::vector<unsigned int> output {};

  for( const auto &vtkfile : vtkfileList )
  {
    SYS_T::file_check( vtkfile );

    VEC_T::insert_end( output, VEC_T::cast_to_unsigned_int( VTK_T::read_int_PointData( vtkfile, "GlobalNodeID" ) ) );
  }

  return output;
}

std::vector<unsigned int> NodalBC_3D_FSI::CA_benchmark_BC( const std::vector<std::string> &vtkfileList ) const
{
  std::vector<unsigned int> output {};

  for( const auto &vtkfile : vtkfileList )
  {
    SYS_T::file_check( vtkfile );
    std::vector<int> node_y0 {};

    // read all node id
    const std::vector<int> node_id = VTK_T::read_int_PointData( vtkfile, "GlobalNodeID" );

    // read all node xyz
    int num_node, num_cell;
    std::vector<int> ien;
    std::vector<double> node_xyz {};

    VTK_T::read_grid( vtkfile, num_node, num_cell, node_xyz, ien );

    // select the node that y=0
    for(int ii{0}; ii < num_node; ++ii)
    {
      const double y_value = node_xyz[3*ii + 1];

      if(std::abs(y_value) < 1.0e-15)
        node_y0.push_back(node_id[ii]);
    }

    // add them to output
    VEC_T::insert_end( output, VEC_T::cast_to_unsigned_int( node_y0 ) );
  }

  return output;
}

// EOF
