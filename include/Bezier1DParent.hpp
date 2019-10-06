#ifndef BEZIER1DPARENT_HPP
#define BEZIER1DPARENT_HPP
// ==================================================================
// Bezier1DParent.hpp
// It gives a 1d Bezier parent/reference element as a derived class
// of ParentElement class.
// 
// Date:
// Oct. 31 2013
// ==================================================================
#include "IQuadPts.hpp"
#include "QuadPts_Gauss.hpp"
#include "Vec_Tools.hpp"
#include "BernsteinBasis.hpp"
#include "ParentElement.hpp"

class Bezier1DParent : public ParentElement
{
  public:
    Bezier1DParent( const int &sdegree,
        const class IQuadPts * const &quaInfo );

    virtual ~Bezier1DParent();

    virtual int get_nLocBas() const {return nLocBas;}
    virtual int get_elemDim() const {return 1;}
    virtual int get_nQuaPts() const {return nQuaPts;}
    virtual int get_numSide() const {return 2;}
    
    virtual int get_xiDegree() const {return xiDegree;}
    virtual int get_etaDegree() const {return -1;}
    virtual int get_zetaDegree() const {return -1;}

    virtual int get_xiQuaPts() const {return nQuaPts;}
    virtual int get_etaQuaPts() const {return -1;}
    virtual int get_zetaQuaPts() const {return -1;}

    virtual double get_quaWeight(const int &qua_index) 
      const {return quaWeight[qua_index]; }
    virtual ParentElement * get_BCParent(const int &sideID) const;

    virtual double evalBasisFunc(const int &node, const int &quapts) const;
    
    virtual double evalBasisFunc_dXi(const int &node, const int &quapts) const;
    virtual double evalBasisFunc_dEta(const int &node, 
        const int &quapts) const {return 0.0;}
    virtual double evalBasisFunc_dZeta(const int &node,
        const int &quapts) const {return 0.0;}
    
    virtual double evalBasisFunc_d2Xi(const int &node, const int &quapts) const;
    virtual double evalBasisFunc_d2Eta(const int &node, 
        const int &quapts) const {return 0.0;}
    virtual double evalBasisFunc_d2Zeta(const int &node,
        const int &quapts) const {return 0.0;}
    virtual double evalBasisFunc_dXiEta(const int &node, 
        const int &quapts) const {return 0.0;}
    virtual double evalBasisFunc_dXiZeta(const int &node,
        const int &quapts) const {return 0.0;}
    virtual double evalBasisFunc_dEtaZeta(const int &node,
        const int &quapts) const {return 0.0;}
 
   virtual void print() const; 
  private:
    int nLocBas;
    int nQuaPts;

    int xiDegree;

    vector<double> basisValue, basisValue_dXi, basisValue_d2Xi;
    vector<double> quaWeight;
};

#endif
