#ifndef IMESHPART_HPP
#define IMESHPART_HPP
// ============================================================================
// IMeshPart.hpp
// This is a light weight implementation of Mesh Partition interface, which is a
// modification of the original IPart class.
//
// Date: Sept. 7 2015
// ============================================================================
#include "IMesh.hpp"
#include "IGlobal_Part.hpp"
#include "Map_Node_Index.hpp"
#include "IIEN.hpp"
#include "NURBS_Bezier.hpp"

class IMeshPart
{
  public:
    IMeshPart(){};
    virtual ~IMeshPart(){};

    // get functions
    virtual int get_nlocalnode() const
    {std::cerr<<"Error: get_nlocalnode is not implemented. \n"; exit(EXIT_FAILURE); return 0;}

    virtual int get_nghostnode() const  
    {std::cerr<<"Error: get_nghostnode is not implemented. \n"; exit(EXIT_FAILURE); return 0;}

    virtual int get_ntotalnode() const  
    {std::cerr<<"Error: get_ntotalnode is not implemented. \n"; exit(EXIT_FAILURE); return 0;}

    virtual int get_nbadnode() const  
    {std::cerr<<"Error: get_nbadnode is not implemented. \n"; exit(EXIT_FAILURE); return 0;}

    virtual int get_nlocghonode() const
    {std::cerr<<"Error: get_nlocghonode is not implemented. \n"; exit(EXIT_FAILURE); return 0;}

    virtual int get_cpu_rank() const
    {std::cerr<<"Error: get_cpu_rank is not implemented. \n"; exit(EXIT_FAILURE); return 0;}

    virtual int get_local_to_global(const int &pos) const
    {std::cerr<<"Error: get_local_to_global is not implemented. \n"; exit(EXIT_FAILURE); return 0;}

    virtual int get_elemLocIndex(const int &index) const
    {std::cerr<<"Error: get_elemLocIndex is not implmeneted. \n"; exit(EXIT_FAILURE); return 0;}
   


    // Detect if node/elem is in partition subdomain
    virtual bool isElemInPart(const int &eindex) const
    {std::cerr<<"Error: isElemInPart is not implemented. \n"; exit(EXIT_FAILURE); return false;}
    
    virtual bool isNodeInPart(const int &nindex) const
    {std::cerr<<"Error: isNodeInPart is not implemented. \n"; exit(EXIT_FAILURE); return false;}

    
    // ------------------------------------------------------------------------
    // Partition routine for each subdomain
    // This routine read in the cpu_rank and decide the data belong to the
    // subdomain.
    // ------------------------------------------------------------------------
    virtual void set_part_info(const int &in_cpu_rank,
        const IMesh * const &mesh,
        const IGlobal_Part * const &gpart,
        const Map_Node_Index * const &mnindex,
        const IIEN * const &IEN,
        const std::vector<double> &ctrlPts,
        const bool &isPrintInfo )
    {std::cerr<<"Error: set_part_info in IMeshPart is not implemented. \n"; exit(EXIT_FAILURE);}


    // ------------------------------------------------------------------------
    // Writer function for the partitioned info
    // ------------------------------------------------------------------------
    virtual void write( const std::string &inputFileName,
        const IMesh * const &mesh,
        const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_x,
        const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_y,
        const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_z ) const
    {std::cerr<<"Error: write in IMeshPart is not implemented. \n"; exit(EXIT_FAILURE);}

};


#endif
