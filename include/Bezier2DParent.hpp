#ifndef BEZIER2DPARENT_HPP
#define BEZIER2DPARENT_HPP
// ==================================================================
// Bezier2DParent.hpp
// It gives a 2d bezier parent reference element as a derived class
// of ParentElement class.
//
// Date:
// Nov. 1 2013
// ==================================================================
#include <vector>
#include <iostream>

#include "IQuadPts.hpp"
#include "Vec_Tools.hpp"
#include "BernsteinBasis.hpp"
#include "ParentElement.hpp"
#include "Bezier1DParent.hpp"
using namespace std;

class Bezier2DParent : public ParentElement
{
  public:
    Bezier2DParent( const int &sdegree, const int &tdegree,
        const class IQuadPts * const &in_quaInfo_s,
        const class IQuadPts * const &in_quaInfo_t );
    virtual ~Bezier2DParent();

    virtual int get_nLocBas() const {return nLocBas;}
    virtual int get_elemDim() const {return 2;}
    virtual int get_nQuaPts() const {return nQuaPts;}
    virtual int get_numSide() const {return 4;}

    virtual int get_xiDegree() const {return xiDegree;}
    virtual int get_etaDegree() const {return etaDegree;}
    virtual int get_zetaDegree() const {return -1;}

    virtual int get_xiQuaPts() const {return nQuaPts_xi;}
    virtual int get_etaQuaPts() const {return nQuaPts_eta;}
    virtual int get_zetaQuaPts() const {return -1;}
    
    virtual double get_quaWeight(const int &qua_index)
      const {return quaWeight[qua_index];}

    virtual ParentElement * get_BCParent(const int &sideID,
        const class IQuadPts * const &in_quaInfo_s,
        const class IQuadPts * const &in_quaInfo_r ) const;
    
    virtual ParentElement * get_BCParent(const int &sideID,
        const class IQuadPts * const &in_quaInfo_s,
        const class IQuadPts * const &in_quaInfo_r,
        const class IQuadPts * const &in_quaInfo_u ) const;

    virtual double evalBasisFunc(const int &node, const int &quapts) 
      const {return basisValue[node * nQuaPts + quapts];}
    
    virtual double evalBasisFunc_dXi(const int &node, const int &quapts) 
      const {return basisValue_dXi[node * nQuaPts + quapts];}
    
    virtual double evalBasisFunc_dEta(const int &node, const int &quapts) 
      const {return basisValue_dEta[node * nQuaPts + quapts];}
    
    virtual double evalBasisFunc_dZeta(const int &node,const int &quapts) 
      const {return 0.0;}

    virtual double evalBasisFunc_d2Xi(const int &node, const int &quapts) 
      const {return basisValue_d2Xi[node * nQuaPts + quapts];}

    virtual double evalBasisFunc_d2Eta(const int &node, const int &quapts)
      const {return basisValue_d2Eta[node * nQuaPts + quapts];}

    virtual double evalBasisFunc_d2Zeta(const int &node, const int &quapts)
      const {return 0.0;}
    
    virtual double evalBasisFunc_dXiEta(const int &node, const int &quapts)
      const {return basisValue_dXiEta[node * nQuaPts + quapts];}
    
    virtual double evalBasisFunc_dXiZeta(const int &node, const int &quapts) 
      const {return 0.0;}
    
    virtual double evalBasisFunc_dEtaZeta(const int &node, const int &quapts) 
      const {return 0.0;}

    virtual void print() const;

  private:
    int nLocBas;
    int nQuaPts, nQuaPts_xi, nQuaPts_eta;

    int xiDegree, etaDegree;

    // The basis function values at quadrature points are stored in a 1d vector.
    // node 0 at qua 0, at qua 1, ..., node 0 at qua nQuaPts, node 1 at qua 0,
    // ... node nLocBas at qua nQuaPts
    vector<double> basisValue, basisValue_dXi, basisValue_d2Xi;
    vector<double> basisValue_dXiEta, basisValue_d2Eta, basisValue_dEta;

    vector<double> quaWeight;

    // pointer to the boundary element, i.e., element of reduced dimension
    mutable ParentElement * BCElem[4];
};
#endif
