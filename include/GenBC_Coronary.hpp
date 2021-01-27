#ifndef GENBC_CORONARY_HPP
#define GENBC_CORONARY_HPP
// ==================================================================
// GenBC_Coronary.hpp
// This class is used to provide outflow boundary conditions for an
// open loop coronary lumped parameter model (LPM). The RCR model is also
// handled by this class in case coronary and RCR BCs are used together.
// The coronary artery is modeled by an RCRCR circuit with intramyocardial
// pressures:
//
// Ra___Ramicro___Rv___Pd
//    |         |
//   Ca        Cim
//              |
//             Pim
//
// Ra, Ramicro and Rv are resistors for the coronary arteries,
// coronary microvasculature and coronary veins, respectively.
// Ca and Cim are capacitors for proximal and distal vascularture respectively.
// Pd and Pim are the distal and intramyocardial pressures, respectively.
// Intramyocardial pressure Pim is applied to capacitor Cim to model restricted coronary
// flow during systole.
//
// The GenBC_Coronary input file uses the following format,
// -----------------------
// Coronary <num_outlets>
// <face_id0> <Ra> <Ca> <Ramicro> <Cim> <Rv> <Pd> <num_Pim_data> <Pim_scaling_coef>
// <t1> <P1>
// <t2> <P2>
//   ....
// <tn> <Pn>
// <face_id1> <Ra> <Ca> <Ramicro> <Cim> <Rv> <Pd> <num_Pim_data> <Pim_scaling_coef>
// <t1> <P1>
// <t2> <P2>
//   ....
// <tn> <Pn>
//-----------------------
// where n=num_Pim_data. Pim_scaling_coef can be used to scale original pressure Pim.
// If an RCR is applied to an outlet, set num_Pim_data=0 and ommit Pim data as follows,
//
// <outlet_id> <Rp> <C> <Rd> <C> <Rd> <Pd> 0 1.0
//
// ==================================================================
#include "IGenBC.hpp"
class GenBC_Coronary : public IGenBC
{
  public:
    GenBC_Coronary( const char * const &lpn_filename, const int &in_N,
       const double &dt3d );

    virtual ~GenBC_Coronary();

    virtual void print_info() const;

    virtual int get_num_ebc() const {return num_ebc;}

    virtual double get_m( const int &ii, const double &in_dot_Q,
       const double &in_Q ) const ;

    virtual double get_n( const int &ii, const double &in_dot_Q,
       const double &in_Q ) const
    {
      return 0.0;
    }

    // Obtain P for the ii-th outlet surface (coronary or RCR).
    virtual double get_P( const int &ii, const double &in_dot_Q,
       const double &in_Q ) const;

    // Get initial P for the ii-th outlet face at the begining of LPM ODE integration.
    virtual double get_P0( const int &ii ) const;

    // Set initial values for the LPM ODE integration.
    virtual void reset_initial_sol( const int &ii, const double &in_Q_0,
       const double &in_P_0, const double &curr_time );

  private:
    const int N;

    const double h; // delta t = Nh

    // Parameters used to define difference quotient for get_m.
    const double absTol, relTol;

    int num_ebc;

    // starting time for integrating the coronary LPM
    double tstart;

    // ending time for integrating the coronary LPM
    double tend;

    // Vectors storing the Ra, Ca, Ra_micro, Cim, Rv, Pd, Pim scaling values for all coronary outlet faces
    // the length of the vectors is num_ebc
    std::vector<double> Ra,Ca,Ra_micro,Cim,Rv,Pd,alpha_Pim;

    // number of intramyocardial pressure Pim data points for each outlet face
    // The vector length is num_ebc, num_Pimdata=0 indicates an RCR outlet.
    std::vector<int> num_Pimdata;

    // tdata and Pimdata for user-provided intramyocardial pressure waveform (time-pressure) for each coronary outlet face
    // Pimderdata for the corresponding dPim/dt for each coronary outlet face.
    std::vector<std::vector<double> > tdata,Pimdata,Pimderdata;

    // precomputed dPim/dt needed by RK4
    std::vector<std::vector<double>> dPimdt_k1,dPimdt_k2,dPimdt_k3;

    // prev_0D_sol records solutions when each ODE integration is completed.
    mutable std::vector<std::vector<double> > prev_0D_sol;

    // Vectors for outlet initial flow and capacitor pressures (2 capacitors)
    std::vector<double> Q0;

    std::vector<std::vector<double>> Pi0;

    // Evaluate the coronary LPM (2 first order ODEs) for the ii-th outlet face (which is a coronary outlet) and output the ODE derivatives to K.
    void F( const int &ii, const double * const &pi, const double &q, const double &dPimdt, double * const &K) const;

    // Evaluate the RCR ODE for the ii-th outlet face and return the result.
    double F( const int &ii, const double &pi, const double &q) const;

    // Pre-compute dPim/dt at the begining of the ODE integration for the ii-th outlet face
    void get_dPimdt( const int &ii);

    // Evaluate the derivatives of a piecewise cubic hermite interpolating polynomial (PCHIP) between points x1 and x2 with values f1, f2 and derivatives d1 and d2
    // for points xe with a size ne. Outputs: fe
    void  cubic_hermite_derivative( const double &x1, const double &x2, const double &f1, const double &f2,
       const double &d1, const double &d2, const int &ne, const std::vector<double> &xe, std::vector<double> &fe );

    // Set PCHIP for Pim data for the ii-th outlet face.
    void set_phcip( const int &ii);

    // Set derivatives for a PCHIP with user provided points and values.
    void spline_pchip_set( const int &n, const std::vector<double> &x, const std::vector<double> &f, std::vector<double> &d);

    // PCHIP sign-testing routine.
    double pchst( const double &arg1, const double &arg2 ) const;

};

#endif
