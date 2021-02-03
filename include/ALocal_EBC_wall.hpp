#ifndef ALOCAL_EBC_WALL_HPP
#define ALOCAL_EBC_WALL_HPP
// ==================================================================
// ALocal_EBC_wall.hpp
//
// Analysis use: local subdomain's elemental boundary condition. This
// is a derived class from ALocal_EBC to add the wall thickness and
// young's modulus, as well as the fluid density used to compute the
// young's modulus.
//
// Author: Ju Liu
// Date: Aug. 10 2020
// ==================================================================
#include "ALocal_EBC.hpp"

class ALocal_EBC_wall : public ALocal_EBC
{
  public:
    ALocal_EBC_wall( const std::string &fileBaseName,
        const int &cpu_rank, const std::string &gname="/ebc" );

    virtual ~ALocal_EBC_wall();

    virtual void print_info() const;

    virtual void get_thickness(const int &eindex, double * const &e_thickness) const;

    virtual void get_youngsmod(const int &eindex, double * const &e_youngsmod) const;

    virtual double get_fluid_density() const {return fluid_density;}

  protected:
    // Fluid densities are used to generate the youngs modulus for arteries at
    // the preprocessing stage.
    double fluid_density;

    // If this partition owns any part of the wall, the thickness and
    // youngsmod vectors are each of length num_local_node[0].
    // Otherwise, these vectors are of length 0.
    std::vector<double> thickness, youngsmod;
};

#endif
