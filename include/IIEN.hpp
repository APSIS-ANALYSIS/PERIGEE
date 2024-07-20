#ifndef IIEN_HPP
#define IIEN_HPP
// ==================================================================
// IIEN.hpp
// This is an interface file for the IEN array classes.
//
// Author: Ju Liu
// Date: Sept. 24th 2013
// ==================================================================

class IIEN
{
  public:
    IIEN() = default;

    virtual ~IIEN() = default;

    // get the IEN arrray for element e at local node ii
    virtual int get_IEN( const int &ee, const int &ii ) const = 0;

    // get the first 4/8 IEN value of a given element ee
    virtual std::array<int,4> get_IEN_array4( const int &ee ) const
    {
      return {{ get_IEN(ee,0), get_IEN(ee,1), get_IEN(ee,2), get_IEN(ee,3) }};
    }

    virtual std::array<int,8> get_IEN_array8( const int &ee ) const
    {
      return {{ get_IEN(ee,0), get_IEN(ee,1), get_IEN(ee,2), get_IEN(ee,3),
                get_IEN(ee,4), get_IEN(ee,5), get_IEN(ee,6), get_IEN(ee,7) }};
    }

    // get the number of local basis functions (per element)
    virtual int get_nLocBas( const int &ee = 0 ) const = 0;

    // print IEN array
    virtual void print_info() const = 0;
};

#endif
