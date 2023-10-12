#ifndef NODALBC_3D_WALL_HPP
#define NODALBC_3D_WALL_HPP
// ============================================================================
// NodalBC_3D_wall.hpp
//
// This is a class that stores the wall nodes with the ring nodes excluded.
//
// Author: Ju Liu
// Datei Created: May 13 2021
// ============================================================================
#include "INodalBC.hpp"
#include "VTK_Tools.hpp"

class NodalBC_3D_wall : public INodalBC
{
  public:
    NodalBC_3D_wall( const std::vector<std::string> &inflow_files,
        const std::string &wall_file,
        const std::vector<std::string> &outflow_files,
        const int &nFunc, const int &elemtype );

    virtual ~NodalBC_3D_wall() = default;

    virtual unsigned int get_dir_nodes(const unsigned int &ii) const
    {return dir_nodes[ii];}

    virtual unsigned int get_per_slave_nodes(const unsigned int &ii) const
    {return per_slave_nodes[ii];}

    virtual unsigned int get_per_master_nodes(const unsigned int &ii) const
    {return per_master_nodes[ii];}

    virtual unsigned int get_num_dir_nodes() const {return num_dir_nodes;}

    virtual unsigned int get_num_per_nodes() const {return num_per_nodes;}

  private:
    std::vector<unsigned int> dir_nodes {};
    unsigned int num_dir_nodes {};

    std::vector<unsigned int> per_slave_nodes {}, per_master_nodes {};
    unsigned int num_per_nodes {};

    NodalBC_3D_wall() = delete;
};

#endif
