#ifndef ALOCAL_IEN_HPP
#define ALOCAL_IEN_HPP
// ==================================================================
// ALocal_IEN.hpp
// Local IEN for Analysis code.
// This class stores the IEN double array for each elements. The IEN
// array is stored in a one-dimensional vector, with stride length
// nLocBas.
//
// Author: Ju Liu
// Date: Nov. 10th 2013
// ==================================================================
#include "HDF5_Reader.hpp"

class ALocal_IEN
{
  public:
    // --------------------------------------------------------------
    // Constructor : read data by specifying the part file base name 
    //               and cpu rank
    // --------------------------------------------------------------
    ALocal_IEN( const std::string &fileBaseName, const int &cpu_rank );

    // Destructor
    virtual ~ALocal_IEN();

    // --------------------------------------------------------------
    // Get stride length
    // This returns the stride length of the IEN array. For certain problems
    // the results is different from the nLocBas stored in AGlobal_Mesh_Info.
    // The logic is this: the one in AGlobal_Mesh_Info is related to the
    // original data in the preprocessor, i.e. the geometry. The one stroed
    // in this class is the number of basis functions for the physics 
    // interpolation. In non-isoparametric elements, the two can be different. 
    // --------------------------------------------------------------
    virtual int get_stride() const {return nLocBas;}

    // --------------------------------------------------------------
    // Get the number of local element. This should be compatible
    // with the one stored in ALocal_Elem.
    // --------------------------------------------------------------
    virtual int get_nlocalele() const {return nlocalele;}

    // Data access functions
    virtual int get_LIEN(const int &elem, const int &node) const
    { return LIEN[elem*nLocBas + node]; }
    
    virtual void get_LIEN_e(const int &elem, std::vector<int> &elem_ien) const
    {
      elem_ien.resize(nLocBas);
      const int offset = elem * nLocBas;
      for(int ii=0; ii<nLocBas; ++ii)
        elem_ien[ii] = LIEN[offset + ii];
      VEC_T::shrink2fit(elem_ien);
    }

    // --------------------------------------------------------------
    // ! get the element e's ien array in an int array: elem_ien.
    //   user is responsible for allocating memory space for elem_ien and
    //   delete it after usage.
    // --------------------------------------------------------------
    virtual void get_LIEN_e(const int &elem, int * const &elem_ien) const
    {
      const int offset = elem * nLocBas;
      for(int ii=0; ii<nLocBas; ++ii)
        elem_ien[ii] = LIEN[offset + ii];
    }

    // --------------------------------------------------------------
    // ! isNode_in_Elem returns a bool value that tells if the node ii
    //   is in the element ee.
    //   Input
    //   ee : the local processor's index for element
    //   ii : the local processor's index for node
    // --------------------------------------------------------------
    virtual bool isNode_in_Elem(const int &elem, const int &node) const
    {
      std::vector<int> eIEN;
      get_LIEN_e(elem, eIEN);
      std::vector<int>::const_iterator it = find(eIEN.begin(), eIEN.end(), node);
      return (it != eIEN.end());
    }

    // Print the info of this class.
    virtual void print_info() const;
  
  protected:
    // The number of local element. This int data is assumed to be also
    // stored in ALocal_Elem.
    int nlocalele; 
    
    // The value can be changed due to the FEM basis function enrichment
    // Users should check if one uses a derived class.
    int nLocBas; 
    
    // The LIEN array may be enriched due to use of non-isoparametric elements.
    // Users should check if one uses a derived class.
    std::vector<int> LIEN; 
};

#endif
