#ifndef NODALBC_SOLID_HPP
#define NODALBC_SOLID_HPP
// ============================================================================
// NodalBC_Solid.hpp
//
// Solid-specific NodalBC with displacement-driving metadata.
//
// Date: Mar. 19 2026
// ============================================================================
#include "NodalBC.hpp"

class NodalBC_Solid : public NodalBC
{
  public:
    NodalBC_Solid( const std::string &vtkfile,
        const int &nFunc,
        const int &in_is_disp_driven )
    : NodalBC( std::vector<std::string>{vtkfile}, nFunc ),
      is_disp_driven( in_is_disp_driven )
    {}

    virtual ~NodalBC_Solid() = default;

    int get_is_disp_driven() const { return is_disp_driven; }

  private:
    const int is_disp_driven;

    NodalBC_Solid() = delete;
};

#endif
