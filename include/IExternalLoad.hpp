#ifndef IEXTERNALLOAD_HPP
#define IEXTERNALLOAD_HPP

#include "Vector_3.hpp"

class IExternalLoad
{
  public:
    IExternalLoad() = default;

    virtual ~IExternalLoad() = default;

    virtual Vector_3 get_f( Vector_3 xyz, double t ) const = 0;

    virtual Vector_3 get_H( int ebc_id, Vector_3 xyz, double t, 
        Vector_3 normal ) const = 0;
};

#endif
