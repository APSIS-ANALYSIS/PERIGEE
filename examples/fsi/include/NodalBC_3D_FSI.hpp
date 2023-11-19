#ifndef NODALBC_3D_FSI_HPP
#define NODALBC_3D_FSI_HPP
// ============================================================================
// NodalBC_3D_FSI.hpp
//
// This is a nodal boundary condition specification class that is specifically
// designed for ALE-FSI problems.
// ============================================================================
#include "INodalBC.hpp"
#include "Tet_Tools.hpp"

class NodalBC_3D_FSI : public INodalBC
{
  public:
    NodalBC_3D_FSI( const std::string &fluid_file,
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
        const int &fsiBC_type );

    NodalBC_3D_FSI( const std::string &fluid_file,
        const int &nFunc,
        const int &fsiBC_type );

    virtual ~NodalBC_3D_FSI() {};

    virtual unsigned int get_dir_nodes(const unsigned int &ii) const
    {return dir_nodes[ii];}

    virtual unsigned int get_per_slave_nodes(const unsigned int &ii) const
    {
      SYS_T::print_fatal("Error: periodic nodes are not defined in NodalBC_3D_CMM.\n");
      return 0;
    }

    virtual unsigned int get_per_master_nodes(const unsigned int &ii) const
    {
      SYS_T::print_fatal("Error: periodic nodes are not defined in NodalBC_3D_CMM.\n");
      return 0;
    }

    virtual unsigned int get_num_dir_nodes() const {return num_dir_nodes;}

    virtual unsigned int get_num_per_nodes() const {return 0;}

  private:
    std::vector<unsigned int> dir_nodes;
    unsigned int num_dir_nodes;

    NodalBC_3D_FSI() {};

    std::vector<unsigned int> get_vtk_nodal_id( const std::vector<std::string> &vtk_file_list_name ) const;

    // For coronary artery benchmark, assume ringBC_type = 0
    std::vector<unsigned int> CA_benchmark_BC( const std::vector<std::string> &vtk_file_list_name ) const;
};

#endif
