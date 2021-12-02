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
        const int &fsiBC_type )
{}

// EOF
