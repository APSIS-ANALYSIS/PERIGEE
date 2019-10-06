#ifndef BEZIER3DPARENT_HPP
#define BEZIER3DPARENT_HPP
// ==================================================================
// Bezier3DParent.hpp
// It gives 3d bezier parent reference element as a derived class of
// ParentElement class.
// 
// Date:
// Nov. 4th 2013
// ==================================================================
#include <vector>
#include <iostream>

#include "IQuadPts.hpp"
#include "Vec_Tools.hpp"
#include "BernsteinBasis.hpp"
#include "ParentElement.hpp"
#include "Bezier2DParent.hpp"

using namespace std;

class Bezier3DParent : public ParentElement
{
  public:
    Bezier3DParent( const int &sdegree,
        const int &tdegree,
        const int &udegree, 
        const class IQuadPts * const &in_quaInfo_s,
        const class IQuadPts * const &in_quaInfo_t,
        const class IQuadPts * const &in_quaInfo_u );
    virtual ~Bezier3DParent();

    virtual int get_nLocBas() const {return nLocBas;}
    virtual int get_elemDim() const {return 3;}
    virtual int get_nQuaPts() const {return nQuaPts;}
    virtual int get_numSide() const {return 6;}

    virtual int get_xiDegree() const {return xiDegree;}
    virtual int get_etaDegree() const {return etaDegree;}
    virtual int get_zetaDegree() const {return zetaDegree;}

    virtual int get_xiQuaPts() const {return nQuaPts_xi;}
    virtual int get_etaQuaPts() const {return nQuaPts_eta;}
    virtual int get_zetaQuaPts() const {return nQuaPts_zeta;}

    virtual double get_quaWeight(const int &qua_index)
      const {return quaWeight[qua_index];}

    virtual ParentElement * get_BCParent(const int &sideID,
        const class IQuadPts * const &in_quaInfo_s,
        const class IQuadPts * const &in_quaInfo_t,
        const class IQuadPts * const &in_quaInfo_u ) const;
    
    virtual ParentElement * get_BCParent(const int &sideID,
        const class IQuadPts * const &in_quaInfo_s,
        const class IQuadPts * const &in_quaInfo_t ) const;

    virtual double evalBasisFunc(const int &node, const int &quapts)
      const {return basisValue[node * nQuaPts + quapts];}

    virtual double evalBasisFunc_dXi(const int &node, const int &quapts)
      const {return basisValue_dXi[node * nQuaPts + quapts];}
    
    virtual double evalBasisFunc_dEta(const int &node, const int &quapts)
      const {return basisValue_dEta[node * nQuaPts + quapts];}
    
    virtual double evalBasisFunc_dZeta(const int &node, const int &quapts)
      const {return basisValue_dZeta[node * nQuaPts + quapts];}

    virtual double evalBasisFunc_d2Xi(const int &node, const int &quapts)
      const {return basisValue_d2Xi[node * nQuaPts + quapts];}
    
    virtual double evalBasisFunc_d2Eta(const int &node, const int &quapts)
      const {return basisValue_d2Eta[node * nQuaPts + quapts];}
    
    virtual double evalBasisFunc_d2Zeta(const int &node, const int &quapts)
      const {return basisValue_d2Zeta[node * nQuaPts + quapts];}

    virtual double evalBasisFunc_dXiEta(const int &node, const int &quapts)
      const {return basisValue_dXiEta[node * nQuaPts + quapts];}
    
    virtual double evalBasisFunc_dXiZeta(const int &node, const int &quapts)
      const {return basisValue_dXiZeta[node * nQuaPts + quapts];}
    
    virtual double evalBasisFunc_dEtaZeta(const int &node, const int &quapts)
      const {return basisValue_dEtaZeta[node * nQuaPts + quapts];}

    virtual void print() const;

  private:
    int nLocBas;
    int nQuaPts, nQuaPts_xi, nQuaPts_eta, nQuaPts_zeta;
    int xiDegree, etaDegree, zetaDegree;

    vector<double> basisValue, basisValue_dXi, basisValue_dEta, basisValue_dZeta;
    vector<double> basisValue_d2Xi, basisValue_d2Eta, basisValue_d2Zeta;
    vector<double> basisValue_dXiEta, basisValue_dXiZeta, basisValue_dEtaZeta;

    vector<double> quaWeight;

    mutable ParentElement * BCElem[6];
};
#endif
