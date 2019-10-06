#ifndef ITMESH_HPP
#define ITMESH_HPP
// ============================================================================
// ITMesh.hpp
// This is an interface file for describing T-spline meshes in multi-dimension.
// This class is pure virtual.
//
// Date: Aug. 14 2015
// ============================================================================
#include <vector>

class ITMesh
{
  public:
    ITMesh(){};
    virtual ~ITMesh(){};

    virtual void print_info() const = 0;

    virtual int get_s_degree() const = 0;
    virtual int get_t_degree() const = 0;
    virtual int get_u_degree() const = 0;

    virtual int get_nFunc() const = 0;
    virtual int get_nElem() const = 0;

    virtual double get_ctrlPts(const int &node, const int &dir) const = 0;
    
    virtual int get_nLocBas( const int &element ) const = 0;
    
    virtual int get_nLocBas_max() const = 0;

    virtual int get_nLocBas_min() const = 0;

    virtual int get_nLocBas_sum() const = 0;

    virtual int get_IEN(const int &element, const int &node) const = 0;

    virtual void get_bExt(const int &element, const int &node, std::vector<double> &ext ) const = 0;
};

#endif
