#ifndef IAEXTRACTOR_HPP
#define IAEXTRACTOR_HPP
// ==================================================================
// IAExtractor.hpp
// Interface for extraction operator. It is designed for reading the
// extraction operator from disk, and distributing the extractor to 
// element routine. A represents fem-analysis.
//
// Date:
// Nov. 8th 2013
// ==================================================================
#include "HDF5_PartReader.hpp"

class IAExtractor
{
  public:
    IAExtractor(){};
    virtual ~IAExtractor(){};

    virtual void get_EXT_x(const int &e, std::vector<double> &ext_x) 
      const {std::cout<<"This function is not implemented. \n";}

    virtual void get_EXT_y(const int &e, std::vector<double> &ext_y) 
      const {std::cout<<"This function is not implemented. \n";}
    
    virtual void get_EXT_z(const int &e, std::vector<double> &ext_z) 
      const {std::cout<<"This function is not implemented. \n";}
    
    // pass extractor by double array. User is responsible for deleting
    // the pointer.
    virtual void get_EXT_x(const int &e, double * &ext_x) 
      const {std::cout<<"This function is not implemented. \n";}
    
    virtual void get_EXT_y(const int &e, double * &ext_y) 
      const {std::cout<<"This function is not implemented. \n";}
    
    virtual void get_EXT_z(const int &e, double * &ext_z) 
      const {std::cout<<"This function is not implemented. \n";}

    // pass full extractor by vector
    virtual void get_EXT( const int &e, const int &ii, std::vector<double> &ext) const 
    { SYS_T::commPrint("Error: get_EXT is not implemented. \n");}

    // pass full extractor by double array. The user is responsible for deleting
    // the dynamic double array
    virtual void get_EXT( const int &e, const int &ii, double * &ext ) const
    { SYS_T::commPrint("Error: get_EXT is not implemented. \n");}


    virtual int get_sdegree() const
    {SYS_T::commPrint("Error: get_sdegree is not implemented. \n"); return 0;}

    virtual int get_tdegree() const
    {SYS_T::commPrint("Error: get_tdegree is not implemented. \n"); return 0;}

    virtual void print_info() const
    { SYS_T::commPrint("Error: print_info is not implemented. \n");}

};


#endif
