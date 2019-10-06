#ifndef IBERNSTEINBASIS_HPP
#define IBERNSTEINBASIS_HPP
// ============================================================================
// IBernsteinBasis.hpp
//
// Interface for Bernstein basis.
//
// Date: Aug 24 2015
// ============================================================================
#include "Sys_Tools.hpp"

class IBernsteinBasis
{
  public:
    IBernsteinBasis(){};

    virtual ~IBernsteinBasis(){};

    // Get polynomial degree
    virtual int get_deg_s() const 
    {SYS_T::commPrint("Error: get_deg_s() is not implemented. \n"); return 0;}
    
    virtual int get_deg_t() const 
    {SYS_T::commPrint("Error: get_deg_t() is not implemented. \n"); return 0;}
    
    virtual int get_deg_u() const 
    {SYS_T::commPrint("Error: get_deg_u() is not implemented. \n"); return 0;}


    // Get the number of quadrature points
    virtual int get_nQuaPts_s() const
    {SYS_T::commPrint("Error: get_nQuaPts_s() is not implemented. \n"); return 0;}
    
    virtual int get_nQuaPts_t() const
    {SYS_T::commPrint("Error: get_nQuaPts_t() is not implemented. \n"); return 0;}
    
    virtual int get_nQuaPts_u() const
    {SYS_T::commPrint("Error: get_nQuaPts_u() is not implemented. \n"); return 0;}


    // Get basis and its derivatives up to the second order
    virtual double get_B(const int &ii, const int &qua) const
    {SYS_T::commPrint("Error: get_B() is not implemented. \n"); return 0.0;}
    
    virtual void get_B(const int &qua, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_B() is not implemented. \n");}

    virtual void get_B(const int &qua_s, const int &qua_t, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_B() is not implemented. \n");}
    
    virtual void get_B(const int &qua_s, const int &qua_t,
       const int &qua_u, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_B() is not implemented. \n");}


    // first-order derivatives
    virtual double get_dB_ds(const int &ii, const int &qua) const
    {SYS_T::commPrint("Error: get_dB_ds() is not implemented. \n"); return 0.0;}
    
    virtual void get_dB_ds(const int &qua, std::vector<double> &val ) const
    {SYS_T::commPrint("Error: get_dB_ds() is not implemented. \n");}

    virtual void get_dB_ds(const int &qua_s, const int &qua_t, std::vector<double> &val ) const
    {SYS_T::commPrint("Error: get_dB_ds() is not implemented. \n");}

    virtual void get_dB_ds(const int &qua_s, const int &qua_t, const int &qua_u, 
        std::vector<double> &val ) const
    {SYS_T::commPrint("Error: get_dB_ds() is not implemented. \n");}

    virtual double get_dB_dt(const int &ii, const int &qua) const
    {SYS_T::commPrint("Error: get_dB_dt() is not implemented. \n"); return 0.0;}

    virtual void get_dB_dt(const int &qua, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_dB_dt() is not implemented. \n");}

    virtual void get_dB_dt(const int &qua_s, const int &qua_t, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_dB_dt() is not implemented. \n");}

    virtual void get_dB_dt(const int &qua_s, const int &qua_t, 
        const int &qua_u, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_dB_dt() is not implemented. \n");}

    virtual double get_dB_du(const int &ii, const int &qua) const
    {SYS_T::commPrint("Error: get_dB_du() is not implemented. \n"); return 0.0;}

    virtual void get_dB_du(const int &qua, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_dB_du() is not implemented. \n");}

    virtual void get_dB_du(const int &qua_s, const int &qua_t, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_dB_du() is not implemented. \n");}

    virtual void get_dB_du(const int &qua_s, const int &qua_t, 
        const int &qua_u, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_dB_du() is not implemented. \n");}

    // second-order derivatives
    virtual double get_d2B_dss(const int &ii, const int &qua) const
    {SYS_T::commPrint("Error: get_d2B_dss() is not implemented. \n"); return 0.0;}

    virtual double get_d2B_dtt(const int &ii, const int &qua) const
    {SYS_T::commPrint("Error: get_d2B_dtt() is not implemented. \n"); return 0.0;}

    virtual double get_d2B_duu(const int &ii, const int &qua) const
    {SYS_T::commPrint("Error: get_d2B_duu() is not implemented. \n"); return 0.0;}

    virtual double get_d2B_dst(const int &ii, const int &qua) const
    {SYS_T::commPrint("Error: get_d2B_dst() is not implemented. \n"); return 0.0;}

    virtual double get_d2B_dsu(const int &ii, const int &qua) const
    {SYS_T::commPrint("Error: get_d2B_dsu() is not implemented. \n"); return 0.0;}

    virtual double get_d2B_dtu(const int &ii, const int &qua) const
    {SYS_T::commPrint("Error: get_d2B_dtu() is not implemented. \n"); return 0.0;}

    // --- get value for all functions in a vector
    virtual void get_d2B_dss(const int &qua, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dss() is not implemented. \n");}

    virtual void get_d2B_dtt(const int &qua, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dtt() is not implemented. \n");}

    virtual void get_d2B_duu(const int &qua, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_duu() is not implemented. \n");}

    virtual void get_d2B_dst(const int &qua, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dst() is not implemented. \n");}

    virtual void get_d2B_dsu(const int &qua, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dsu() is not implemented. \n");}

    virtual void get_d2B_dtu(const int &qua, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dtu() is not implemented. \n");}


    // --- get value with qua index in s t dirs
    virtual void get_d2B_dss(const int &qua_s, const int &qua_t, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dss() is not implemented. \n");}

    virtual void get_d2B_dtt(const int &qua_s, const int &qua_t, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dtt() is not implemented. \n");}

    virtual void get_d2B_duu(const int &qua_s, const int &qua_t, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_duu() is not implemented. \n");}

    virtual void get_d2B_dst(const int &qua_s, const int &qua_t, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dst() is not implemented. \n");}

    virtual void get_d2B_dsu(const int &qua_s, const int &qua_t, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dsu() is not implemented. \n");}

    virtual void get_d2B_dtu(const int &qua_s, const int &qua_t, std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dtu() is not implemented. \n");}



    // --- get value with qua index in s t u dirs
    virtual void get_d2B_dss(const int &qua_s, const int &qua_t, const int &qua_u,
        std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dss() is not implemented. \n");}

    virtual void get_d2B_dtt(const int &qua_s, const int &qua_t, const int &qua_u,
        std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dtt() is not implemented. \n");}

    virtual void get_d2B_duu(const int &qua_s, const int &qua_t, const int &qua_u,
        std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_duu() is not implemented. \n");}

    virtual void get_d2B_dst(const int &qua_s, const int &qua_t, const int &qua_u,
        std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dst() is not implemented. \n");}

    virtual void get_d2B_dsu(const int &qua_s, const int &qua_t, const int &qua_u,
        std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dsu() is not implemented. \n");}

    virtual void get_d2B_dtu(const int &qua_s, const int &qua_t, const int &qua_u,
        std::vector<double> &val) const
    {SYS_T::commPrint("Error: get_d2B_dtu() is not implemented. \n");}

    // print value
    virtual void print_info() const
    {SYS_T::commPrint("Error: print_info() is not implemented.\n");}

};

#endif
