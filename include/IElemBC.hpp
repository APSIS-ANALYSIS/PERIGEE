#ifndef IELEMBC_HPP
#define IELEMBC_HPP
// ============================================================================
// IElemBC.hpp
//
// Object:
// This is a base class that provides an interface for elemental
// boundary condition specification. It provides the global element index for
// weak enforcement boundary integral. In specific, it will return the index of
// the element whose front/back/left/right/top/bottom face needs weak boundary
// integral, and the number of these type elements.
//
// Date: Aug 17 2015
// ============================================================================

#include "Vec_Tools.hpp"

class IElemBC
{
  public:
    // Default Constructor : empty xxx_elem vectors
    IElemBC();

    virtual ~IElemBC();

    // ------------------------------------------------------------------------
    // get_num_xxx_elem()
    // returns the number of elements whose xxx face needs weak boundary
    // integration. xxx = front back left right top bottom.
    // ------------------------------------------------------------------------
    virtual unsigned int get_num_front_elem() const
    {return front_elem.size();}

    virtual unsigned int get_num_back_elem() const
    {return back_elem.size();}

    virtual unsigned int get_num_left_elem() const
    {return left_elem.size();}

    virtual unsigned int get_num_right_elem() const
    {return right_elem.size();}

    virtual unsigned int get_num_top_elem() const
    {return top_elem.size();}

    virtual unsigned int get_num_bottom_elem() const
    {return bottom_elem.size();}



    // ------------------------------------------------------------------------
    // get_xxx_elem(ii)
    // returns the xxx-face integral element's global index. The parameter ii
    // runs as 0 <= ii < get_num_xxx_elem(). 
    // ------------------------------------------------------------------------
    virtual unsigned int get_front_elem(const unsigned int &ii) const
    {return front_elem[ii];}

    virtual unsigned int get_back_elem(const unsigned int &ii) const
    {return back_elem[ii];}

    virtual unsigned int get_left_elem(const unsigned int &ii) const
    {return left_elem[ii];}

    virtual unsigned int get_right_elem(const unsigned int &ii) const
    {return right_elem[ii];}

    virtual unsigned int get_top_elem(const unsigned int &ii) const
    {return top_elem[ii];}

    virtual unsigned int get_bottom_elem(const unsigned int &ii) const
    {return bottom_elem[ii];}


    virtual void print_info() const;


  protected:
    std::vector<unsigned int> left_elem, right_elem;
    std::vector<unsigned int> top_elem, bottom_elem;
    std::vector<unsigned int> front_elem, back_elem;

    // --------------------------------------------------------------
    // Generate_BCElem_2D 
    // output the boundary element with different facial orientations. 
    // For an element mesh 
    //  t   6 7 8  | Fro: 2 5 8
    //  ^   3 4 5  | Bac: 0 3 6
    //  |   0 1 2  | Lef: 0 1 2
    // 0----> s    | Rig: 6 7 8
    // --------------------------------------------------------------
    void Generate_BCElem_2D( const int &nElem_x, const int &nElem_y,
        std::vector<unsigned int> &front, std::vector<unsigned int> &back,
        std::vector<unsigned int> &left, std::vector<unsigned int> &right );

    
    // -------------------------------------------------------------------------
    // Similar to the function for nodes, but list the elements that has faces
    // on the boundary. 
    // -------------------------------------------------------------------------
    void Generate_BCElem_3D_A( const int &nElem_x, const int &nElem_y,
        const int &nElem_z,
        std::vector<unsigned int> &front, std::vector<unsigned int> &back, 
        std::vector<unsigned int> &left, std::vector<unsigned int> &right, 
        std::vector<unsigned int> &top, std::vector<unsigned int> &bottom ) const;


    // ------------------------------------------------------------------------
    // Generate_BCElems_B is a similar function to the _A version. It takes into
    // account of the nonzero starting element index, which is useful for
    // multipatch cases.
    // ------------------------------------------------------------------------
    void Generate_BCElem_3D_B( const int &nElem_x, const int &nElem_y,
        const int &nElem_z, const int &es,
        std::vector<unsigned int> &front, std::vector<unsigned int> &back, 
        std::vector<unsigned int> &left, std::vector<unsigned int> &right, 
        std::vector<unsigned int> &top, std::vector<unsigned int> &bottom ) const;


};

#endif
