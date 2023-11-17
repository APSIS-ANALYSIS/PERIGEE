#include "QuadPts_Gauss_Triangle.hpp"

QuadPts_Gauss_Triangle::QuadPts_Gauss_Triangle( const int &in_num_pts ) : num_pts( in_num_pts )
{
  qp.resize( 3 * num_pts );
  qw.resize( num_pts );

  int idx = 0;
  double a = 0.0, b = 0.0, c = 0.0, w = 0.0;

  switch( num_pts )
  {
    case 3:
      qp[0] = 0.5;
      qp[1] = 0.5;
      qp[2] = 0.0;
      qw[0] = 1.0 / 3.0;

      qp[3] = 0.5;
      qp[4] = 0.0;
      qp[5] = 0.5;
      qw[1] = 1.0 / 3.0;

      qp[6] = 0.0;
      qp[7] = 0.5;
      qp[8] = 0.5;
      qw[2] = 1.0 / 3.0;
      break;

    case 4:
      qp[0] = 1.0 / 3.0;
      qp[1] = 1.0 / 3.0;
      qp[2] = 1.0 / 3.0;
      qw[0] = -0.5625;

      w = 0.520833333333333;
      qp[3] = 0.6;
      qp[4] = 0.2;
      qp[5] = 0.2;
      qw[1] = w;

      qp[6] = 0.2;
      qp[7] = 0.6;
      qp[8] = 0.2;
      qw[2] = w;

      qp[9] = 0.2;
      qp[10] = 0.2;
      qp[11] = 0.6;
      qw[3] = w;
      break;

    case 6:
      a = 0.816847572980459;
      b = 0.091576213509771;
      w = 0.109951743655322;

      qp[0] = a; qp[1] = b; qp[2] = b; qw[0] = w;
      qp[3] = b; qp[4] = a; qp[5] = b; qw[1] = w;
      qp[6] = b; qp[7] = b; qp[8] = a; qw[2] = w;

      a = 0.108103018168070;
      b = 0.445948490915965;
      w = 0.223381589678011;
      qp[9] = a; qp[10] = b; qp[11] = b; qw[3] = w;
      qp[12] = b; qp[13] = a; qp[14] = b; qw[4] = w;
      qp[15] = b; qp[16] = b; qp[17] = a; qw[5] = w;
      break;

    case 13:
      qp[0] = 1.0 / 3.0;
      qp[1] = 1.0 / 3.0;
      qp[2] = 1.0 / 3.0;
      qw[0] = -0.149570044467670;

      a = 0.479308067841923;
      b = 0.260345966079038;
      w = 0.175615257433204;
      qp[3] = a; qp[4] = b; qp[5] = b; qw[1] = w;
      qp[6] = b; qp[7] = a; qp[8] = b; qw[2] = w;
      qp[9] = b; qp[10] = b; qp[11] = a; qw[3] = w;
      
      a = 0.869739794195568;
      b = 0.065130102902216;
      w = 0.053347235608839;
      qp[12] = a; qp[13] = b; qp[14] = b; qw[4] = w;
      qp[15] = b; qp[16] = a; qp[17] = b; qw[5] = w;
      qp[18] = b; qp[19] = b; qp[20] = a; qw[6] = w;

      a = 0.638444188569809;
      b = 0.312865496004875;
      c = 1.0 - a - b;
      w = 0.077113760890257;
      qp[21] = a; qp[22] = b; qp[23] = c; qw[7] = w;
      qp[24] = a; qp[25] = c; qp[26] = b; qw[8] = w;
      qp[27] = b; qp[28] = a; qp[29] = c; qw[9] = w;
      qp[30] = c; qp[31] = a; qp[32] = b; qw[10] = w;
      qp[33] = b; qp[34] = c; qp[35] = a; qw[11] = w;
      qp[36] = c; qp[37] = b; qp[38] = a; qw[12] = w;

      break;

    case 19:
      qp[0] = 1.0 / 3.0; qp[1] = 1.0 / 3.0; qp[2] = 1.0 / 3.0;
      qw[0] = 0.037861091200315;

      a = 0.797426985353087;
      b = 0.101286507323456;
      w = 0.037620425413183;
      qp[3] = a; qp[4] = b; qp[5] = b; qw[1] = w;
      qp[6] = b; qp[7] = a; qp[8] = b; qw[2] = w;
      qp[9] = b; qp[10] = b; qp[11] = a; qw[3] = w;

      a = 0.059715871789770;
      b = 0.470142064105115;
      w = 0.078357352244117;
      qp[12] = a; qp[13] = b; qp[14] = b; qw[4] = w;
      qp[15] = b; qp[16] = a; qp[17] = b; qw[5] = w;
      qp[18] = b; qp[19] = b; qp[20] = a; qw[6] = w;  

      a = 0.535795346449899;
      b = 0.232102326775050;
      w = 0.116271479656966;
      qp[21] = a; qp[22] = b; qp[23] = b; qw[7] = w;
      qp[24] = b; qp[25] = a; qp[26] = b; qw[8] = w;
      qp[27] = b; qp[28] = b; qp[29] = a; qw[9] = w;

      a = 0.941038278231121;
      b = 0.029480860884440;
      w = 0.013444267375166;
      qp[30] = a; qp[31] = b; qp[32] = b; qw[10] = w;
      qp[33] = b; qp[34] = a; qp[35] = b; qw[11] = w;
      qp[36] = b; qp[37] = b; qp[38] = a; qw[12] = w;

      a = 0.738416812340510;
      b = 0.232102326775050;
      c = 0.029480860884440;
      w = 0.037509722455232;
      qp[39] = a; qp[40] = b; qp[41] = c; qw[13] = w;
      qp[42] = a; qp[43] = c; qp[44] = b; qw[14] = w;
      qp[45] = b; qp[46] = a; qp[47] = c; qw[15] = w;
      qp[48] = c; qp[49] = a; qp[50] = b; qw[16] = w;
      qp[51] = b; qp[52] = c; qp[53] = a; qw[17] = w;
      qp[54] = c; qp[55] = b; qp[56] = a; qw[18] = w;
      break;

    case 37:
      idx = 0;
      qp[3*idx] = 1.0 / 3.0; qp[3*idx+1] = 1.0 / 3.0; qp[3*idx+2] = 1.0 / 3.0;
      qw[idx] = 0.051739766065744133555179145422;

      a = 0.950275662924105565450352089520;
      b = 0.024862168537947217274823955239;
      w = 0.008007799555564801597804123460;
      idx += 1;
      qp[3*idx] = a; qp[3*idx+1] = b; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = a; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = b; qp[3*idx+2] = a; qw[idx] = w;
  
      a = 0.171614914923835347556304795551;
      b = 0.414192542538082326221847602214;
      w = 0.046868898981821644823226732071;
      idx += 1;
      qp[3*idx] = a; qp[3*idx+1] = b; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = a; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = b; qp[3*idx+2] = a; qw[idx] = w;
 
      a = 0.539412243677190440263092985511;
      b = 0.230293878161404779868453507244;
      w = 0.046590940183976487960361770070; 
      idx += 1;
      qp[3*idx] = a; qp[3*idx+1] = b; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = a; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = b; qp[3*idx+2] = a; qw[idx] = w;

      a = 0.772160036676532561750285570113;
      b = 0.113919981661733719124857214943;
      w = 0.031016943313796381407646220131; 
      idx += 1;
      qp[3*idx] = a; qp[3*idx+1] = b; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = a; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = b; qp[3*idx+2] = a; qw[idx] = w;

      a = 0.009085399949835353883572964740;
      b = 0.495457300025082323058213517632; 
      w = 0.010791612736631273623178240136;
      idx += 1;
      qp[3*idx] = a; qp[3*idx+1] = b; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = a; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = b; qp[3*idx+2] = a; qw[idx] = w;

      a = 0.062277290305886993497083640527;
      b = 0.468861354847056503251458179727;
      w = 0.032195534242431618819414482205;
      idx += 1;
      qp[3*idx] = a; qp[3*idx+1] = b; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = a; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = b; qp[3*idx+2] = a; qw[idx] = w;

      a = 0.851306504174348550389457672223;
      b = 0.022076289653624405142446876931;
      c = 0.126617206172027096933163647918;
      w = 0.015445834210701583817692900053;
      idx += 1;
      qp[3*idx] = a; qp[3*idx+1] = b; qp[3*idx+2] = c; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = a; qp[3*idx+1] = c; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = a; qp[3*idx+2] = c; qw[idx] = w;

      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = c; qp[3*idx+2] = a; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = c; qp[3*idx+1] = a; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = c; qp[3*idx+1] = b; qp[3*idx+2] = a; qw[idx] = w;

      a = 0.018620522802520968955913511549;
      b = 0.689441970728591295496647976487;
      c = 0.291937506468887771754472382212;  
      w = 0.017822989923178661888748319485;
      idx += 1;
      qp[3*idx] = a; qp[3*idx+1] = b; qp[3*idx+2] = c; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = a; qp[3*idx+1] = c; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = a; qp[3*idx+2] = c; qw[idx] = w;

      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = c; qp[3*idx+2] = a; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = c; qp[3*idx+1] = a; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = c; qp[3*idx+1] = b; qp[3*idx+2] = a; qw[idx] = w;

      a = 0.267625659273967961282458816185;
      b = 0.635867859433872768286976979827;
      c = 0.096506481292159228736516560903;
      w = 0.037038683681384627918546472190;
      idx += 1;
      qp[3*idx] = a; qp[3*idx+1] = b; qp[3*idx+2] = c; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = a; qp[3*idx+1] = c; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = a; qp[3*idx+2] = c; qw[idx] = w;

      idx += 1;
      qp[3*idx] = b; qp[3*idx+1] = c; qp[3*idx+2] = a; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = c; qp[3*idx+1] = a; qp[3*idx+2] = b; qw[idx] = w;
      
      idx += 1;
      qp[3*idx] = c; qp[3*idx+1] = b; qp[3*idx+2] = a; qw[idx] = w;
      break;

    default:
      SYS_T::print_fatal("Error: QuadPts_Gauss_Triangle: input number of quadrature points is not implemented. \n");
      break;
  }

  VEC_T::shrink2fit(qp);
  VEC_T::shrink2fit(qw);

  // Correct the formula by scaling all weights by 0.5
  for(int ii=0; ii<num_pts; ++ii) qw[ii] *= 0.5;
}

QuadPts_Gauss_Triangle::~QuadPts_Gauss_Triangle()
{
  VEC_T::clean(qp);
  VEC_T::clean(qw);
}

void QuadPts_Gauss_Triangle::print_info() const
{
  SYS_T::commPrint("====== Gauss Points for Triangle =======\n");
  SYS_T::commPrint("Number of points = %d\n", num_pts);
  SYS_T::commPrint("qp.size() = %d\n", qp.size());
  SYS_T::commPrint("qw.size() = %d\n", qw.size());
  for(int ii=0; ii<num_pts; ++ii)
    SYS_T::commPrint("  %.15f %.15f %.15f %.15f\n", 
        qp[3*ii], qp[3*ii+1], qp[3*ii+2], qw[ii]);
  SYS_T::commPrint("========================================\n");
}

// EOF
