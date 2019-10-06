#ifndef PARENTELEMENT_HPP
#define PARENTELEMENT_HPP
// ==================================================================
// ParentElement.hpp
// This is an abstract interface for parent element. For a specific
// problem, the parent element should be able to give the following:
// (1) elemDim = 1, 2, or 3;
// (2) nLocBas: the number of local basis functions for this 
//              type element;
// (3) xiDegree, etaDegree, zetaDegree: the degree of each in xi, 
//     eta, zeta directions;
// (4) xiQuaPts, etaQuaPts, zetaQuaPts: the number of quadrature 
//     points in each direction;
// (5) numSides: the number of side lower-dimension element. 
// Date:
// Oct. 31 2013
// ==================================================================

#include <vector>
#include <cstdlib>
#include <iostream>

#include "IQuadPts.hpp"
using namespace std;

class ParentElement
{
  public:
    ParentElement(){};
    virtual ~ParentElement(){};

    // basic data
    virtual int get_nLocBas() const = 0;
    virtual int get_elemDim() const = 0;
    virtual int get_nQuaPts() const = 0;
    virtual int get_numSide() const = 0;
    
    // degree
    virtual int get_xiDegree() const = 0;
    virtual int get_etaDegree() const = 0;
    virtual int get_zetaDegree() const = 0;

    // quadrature info
    virtual int get_xiQuaPts() const = 0;
    virtual int get_etaQuaPts() const = 0;
    virtual int get_zetaQuaPts() const = 0;
    
    // get quadrature weight
    virtual double get_quaWeight( const int &qua_index ) const = 0; 

    // side parent element
    virtual ParentElement* get_BCParent(const int &sideID) const {return NULL;}
    virtual ParentElement* get_BCParent(const int &sideID,
        const class IQuadPts * const &in_quaInfo_s,
        const class IQuadPts * const &in_quaInfo_t ) const {return NULL;}
    virtual ParentElement* get_BCParent(const int &sideID,
        const class IQuadPts * const &in_quaInfo_s,
        const class IQuadPts * const &in_quaInfo_t,
        const class IQuadPts * const &in_quaInfo_u ) const {return NULL;}

    // local basis function evaluation
    virtual double evalBasisFunc(const int &node, const int &quapts) const = 0;
    
    virtual double evalBasisFunc_dXi(const int &node, const int &quapts) const = 0;
    virtual double evalBasisFunc_dEta(const int &node, const int &quapts) const = 0;
    virtual double evalBasisFunc_dZeta(const int &node, const int &quapts) const = 0;
    
    virtual double evalBasisFunc_d2Xi(const int &node, const int &quapts) const = 0;
    virtual double evalBasisFunc_d2Eta(const int &node, const int &quapts) const = 0;
    virtual double evalBasisFunc_d2Zeta(const int &node, const int &quapts) const = 0;
    virtual double evalBasisFunc_dXiEta(const int &node, const int &quapts) const = 0;
    virtual double evalBasisFunc_dXiZeta(const int &node, const int &quapts) const = 0;
    virtual double evalBasisFunc_dEtaZeta(const int &node, const int &quapts) const = 0;

    // print function
    virtual void print() const {cout<<"Print() for parent element: Not implemented. \n";}
};
#endif
