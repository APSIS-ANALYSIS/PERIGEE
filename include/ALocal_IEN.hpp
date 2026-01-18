#ifndef ALOCAL_IEN_HPP
#define ALOCAL_IEN_HPP
// ============================================================================
// ALocal_IEN.hpp
//
// Local IEN array for Analysis code.
// This class stores the IEN array for volumetric elements. The IEN array is 
// stored in a one-dimensional vector, with stride length nLocBas. For the ee-th
// element, its ii-th node has volumetric mesh node index 
//                    LIEN[ ee * nLocBas + ii ].
// Local means the IEN array returns the nodal indices in the local_to_global
// array, that is the local + ghost nodes belonging to the CPU subdomain.
//
// Author: Ju Liu
// Date: Nov. 10th 2013
// ============================================================================
#include "HDF5_Reader.hpp"

class ALocal_IEN
{
  public:
    // ------------------------------------------------------------------------
    // ! Constructor : read data by specifying the part file base name and cpu 
    //   rank
    // ------------------------------------------------------------------------
    ALocal_IEN( const std::string &fileBaseName, int cpu_rank );

    ALocal_IEN( const HDF5_Reader * const &h5r );

    // ------------------------------------------------------------------------
    // ! Destructor
    // ------------------------------------------------------------------------
    virtual ~ALocal_IEN() = default;

    // ------------------------------------------------------------------------
    // ! Get stride length
    //   This returns the stride length of the IEN array. For problems with
    //   non-isoparametric elements, this result is different from the nLocBas
    //   stored in AGlobal_Mesh_Info. The nLocBas stored in AGlobal_Mesh_Info
    //   is related to the original data in the preprocessor, i.e. the geometry.
    //   The one stored in this class is the number of basis functions for the 
    //   physics interpolation.
    // ------------------------------------------------------------------------
    virtual int get_stride() const {return nLocBas;}

    // ------------------------------------------------------------------------
    // ! Get the number of local elements. This value should be compatible
    //   with the one stored in the ALocal_Elem class.
    // ------------------------------------------------------------------------
    virtual int get_nlocalele() const {return nlocalele;}

    // ------------------------------------------------------------------------
    // ! Data access functions
    // ------------------------------------------------------------------------
    virtual int get_LIEN(int elem, int node) const
    { return LIEN[elem*nLocBas + node]; }
    
    virtual std::vector<int> get_LIEN( int ee ) const
    {
      std::vector<int> elem_ien(nLocBas, 0);
      for(int ii=0; ii<nLocBas; ++ii) elem_ien[ii] = LIEN[ee * nLocBas + ii];
    
      return elem_ien;
    }

    // ------------------------------------------------------------------------
    // ! get the element e's ien array in an int array: elem_ien.
    //   User is responsible for allocating and deallocating memory for elem_ien.
    // ------------------------------------------------------------------------
    virtual void get_LIEN(int ee, int * const &elem_ien) const
    {
      for(int ii=0; ii<nLocBas; ++ii) elem_ien[ii] = LIEN[ee * nLocBas + ii];
    }

    // ------------------------------------------------------------------------
    // ! isNode_in_Elem returns a bool value that tells if the node with volume
    //   mesh node index ii is in the element ee or not.
    //   Input
    //   ee : the local processor's index for element
    //   ii : the local processor's index for node
    // ------------------------------------------------------------------------
    virtual bool isNode_in_Elem(int elem, int node) const
    {
      const std::vector<int> eIEN = get_LIEN( elem );
      std::vector<int>::const_iterator it = find(eIEN.begin(), eIEN.end(), node);
      return (it != eIEN.end());
    }

    // ------------------------------------------------------------------------
    // ! Print the info of this class.
    // ------------------------------------------------------------------------
    virtual void print_info() const;
  
  protected:
    // ------------------------------------------------------------------------
    // The number of local element. 
    // This integer data is assumed to be also stored in the ALocal_Elem class.
    // ------------------------------------------------------------------------
    int nlocalele; 
    
    // ------------------------------------------------------------------------
    // The number of local/elemental basis functions.
    // The value can be changed due to the FEM basis function enrichment
    // Users should check if one uses a derived class.
    // ------------------------------------------------------------------------
    int nLocBas; 
    
    // ------------------------------------------------------------------------
    // The LIEN array may be enriched due to use of non-isoparametric elements.
    // Users should check if one uses a derived class.
    // ------------------------------------------------------------------------
    std::vector<int> LIEN {};

    // Disallow default constructor
    ALocal_IEN() = delete;
};

#endif
