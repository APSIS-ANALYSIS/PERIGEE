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
#include "Tet_Tools.hpp"

class NodalBC_3D_wall : public INodalBC
{
  public:
    NodalBC_3D_wall( const std::vector<std::string> &inflow_files,
        const std::string &wall_file,
        const std::vector<std::string> &outflow_files,
        const int &nFunc, const int &elemtype = 501 );

    virtual ~NodalBC_3D_wall() {};
};

#endif
