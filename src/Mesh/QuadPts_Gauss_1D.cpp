#include "QuadPts_Gauss_1D.hpp"

QuadPts_Gauss_1D::QuadPts_Gauss_1D( const int &in_num_pts, const double &min, 
    const double &max ) : num_pts(in_num_pts)
{
  // Make sure that min < max
  SYS_T::print_fatal_if( min >= max, "Error: QuadPts_Gauss_1D, the given range of quadrature domain is incorrect.\n");

  // Make sure that the containers are empty
  qp.clear(); qw.clear();

  // Generate the rule for [0, 1] domain
  switch( num_pts )
  {
    case 1:
      qp.push_back(0.5);
      qw.push_back(1.0);
      break;
    case 2:
      qp.push_back(0.788675134594813);
      qp.push_back(0.211324865405187);
      qw.push_back(0.5);
      qw.push_back(0.5);
      break;
    case 3:
      qp.push_back(0.887298334620742);
      qp.push_back(0.5);
      qp.push_back(0.112701665379258);
      qw.push_back(0.277777777777777);
      qw.push_back(0.444444444444444);
      qw.push_back(0.277777777777777);
      break;
    case 4:
      qp.push_back(0.930568155797026);
      qp.push_back(0.669990521792428);
      qp.push_back(0.330009478207572);
      qp.push_back(0.069431844202974);
      qw.push_back(0.173927422568727);
      qw.push_back(0.326072577431273);
      qw.push_back(0.326072577431273);
      qw.push_back(0.173927422568727);
      break;
    case 5:
      qp.push_back(0.953089922969332);
      qp.push_back(0.769234655052841);
      qp.push_back(0.500000000000000);
      qp.push_back(0.230765344947158);
      qp.push_back(0.046910077030668);
      qw.push_back(0.118463442528095);
      qw.push_back(0.239314335249683);
      qw.push_back(0.284444444444444);
      qw.push_back(0.239314335249683);
      qw.push_back(0.118463442528095);
      break;
    case 6:
      qp.push_back(0.966234757101576);
      qp.push_back(0.830604693233132);
      qp.push_back(0.619309593041598);
      qp.push_back(0.380690406958402);
      qp.push_back(0.169395306766868);
      qp.push_back(0.033765242898424);
      qw.push_back(0.085662246189585);
      qw.push_back(0.180380786524069);
      qw.push_back(0.233956967286345);
      qw.push_back(0.233956967286345);
      qw.push_back(0.180380786524069);
      qw.push_back(0.085662246189585);
      break;
    case 7:
      qp.push_back(0.974553956171379);
      qp.push_back(0.870765592799697);
      qp.push_back(0.702922575688699);
      qp.push_back(0.500000000000000);
      qp.push_back(0.297077424311301);
      qp.push_back(0.129234407200303);
      qp.push_back(0.025446043828621);
      qw.push_back(0.064742483084435);
      qw.push_back(0.139852695744638);
      qw.push_back(0.190915025252560);
      qw.push_back(0.208979591836735);
      qw.push_back(0.190915025252560);
      qw.push_back(0.139852695744638);
      qw.push_back(0.064742483084435);
      break;
    case 8:
      qp.push_back(0.980144928248768);
      qp.push_back(0.898333238706813);
      qp.push_back(0.762766204958164);
      qp.push_back(0.591717321247825);
      qp.push_back(0.408282678752175);
      qp.push_back(0.237233795041836);
      qp.push_back(0.101666761293187);
      qp.push_back(0.019855071751232);
      qw.push_back(0.050614268145188);
      qw.push_back(0.111190517226687);
      qw.push_back(0.156853322938944);
      qw.push_back(0.181341891689181);
      qw.push_back(0.181341891689181);
      qw.push_back(0.156853322938944);
      qw.push_back(0.111190517226687);
      qw.push_back(0.050614268145188);
      break;

    case 9:
      qp.push_back(0.984080119753813);
      qp.push_back(0.918015553663318);
      qp.push_back(0.806685716350295);
      qp.push_back(0.662126711701905);
      qp.push_back(0.500000000000000);
      qp.push_back(0.337873288298096);
      qp.push_back(0.193314283649705);
      qp.push_back(0.081984446336682);
      qp.push_back(0.015919880246187);

      qw.push_back(0.040637194180787);
      qw.push_back(0.090324080347429);
      qw.push_back(0.130305348201468);
      qw.push_back(0.156173538520001);
      qw.push_back(0.165119677500630);
      qw.push_back(0.156173538520001);
      qw.push_back(0.130305348201468);
      qw.push_back(0.090324080347429);
      qw.push_back(0.040637194180787);
      break;
   
    case 10:
      qp.push_back(0.986953264258586);
      qp.push_back(0.932531683344492);
      qp.push_back(0.839704784149512);
      qp.push_back(0.716697697064624); 
      qp.push_back(0.574437169490816);
      qp.push_back(0.425562830509184);
      qp.push_back(0.283302302935376);
      qp.push_back(0.160295215850488);
      qp.push_back(0.067468316655508);
      qp.push_back(0.013046735741414);

      qw.push_back(0.033335672154344);
      qw.push_back(0.074725674575290);
      qw.push_back(0.109543181257991);
      qw.push_back(0.134633359654998);
      qw.push_back(0.147762112357376);
      qw.push_back(0.147762112357376);
      qw.push_back(0.134633359654998);
      qw.push_back(0.109543181257991);
      qw.push_back(0.074725674575290);
      qw.push_back(0.033335672154344);
      break;

    default:
      compute_npts();
      break;
  }
  VEC_T::shrink2fit(qp);
  VEC_T::shrink2fit(qw);

  // Reverse the points and weights to make them in ascending order
  std::reverse(qp.begin(), qp.end());
  std::reverse(qw.begin(), qw.end());
  
  // Now map the rule to the [min, max] domain via a linear change of interval
  // x = (max-min) xi + min
  // see en.wikipedia.org/wiki/Gaussian_quadrature
  for( int ii =0; ii<num_pts; ++ii )
  {
    qp[ii] = (max - min) * qp[ii] + min;
    qw[ii] = (max - min) * qw[ii];
  }
}

void QuadPts_Gauss_1D::print_info() const
{
  std::cout<<"====== Gauss Points ======="<<std::endl;
  std::cout<<"Num of pt = "<<num_pts<<std::endl;
  std::cout<<"qp.size() = "<<qp.size()<<std::endl;
  std::cout<<"qw.size() = "<<qw.size()<<std::endl;
  std::cout<<"qp \t qw \n";
  for(int ii=0; ii<num_pts; ++ii)
    std::cout<<std::setprecision(16)<<qp[ii]<<'\t'<<qw[ii]<<'\n';
  std::cout<<"=========================="<<std::endl;
}

void QuadPts_Gauss_1D::compute_npts()
{
  const unsigned int n = num_pts;
  const unsigned int m = (num_pts + 1) / 2;
  qp.resize(num_pts);
  qw.resize(num_pts);

  const long double long_double_eps = static_cast<long double>(std::numeric_limits<long double>::epsilon()), double_eps = static_cast<long double>(std::numeric_limits<double>::epsilon());

  // now check whether long double is more accurate than double, and set
  // tolerances accordingly. generate a one that really is generated at run-time
  // and is not optimized away by the compiler. that makes sure that the
  // tolerance is set at run-time with the current behavior, not at compile-time
  // (not doing so leads to trouble with valgrind for example).
  volatile long double runtime_one = 1.0;
  const long double tolerance
    = (runtime_one + long_double_eps != runtime_one
        ?
        std::max (double_eps / 100, long_double_eps * 5)
        :
        double_eps * 5
      );


  for (unsigned int i=1; i<=m; ++i)
  {
    long double z = std::cos(MATH_T::PI * (i-.25)/(n+.5));

    long double pp {1.0}, p1 {1.0};

    // Newton iteration
    do
    {
      // compute L_n (z)
      p1 = 1.;
      long double p2 = 0.;
      for (unsigned int j=0; j<n; ++j)
      {
        const long double p3 = p2;
        p2 = p1;
        p1 = ((2.*j+1.)*z*p2-j*p3)/(j+1);
      }
      pp = n*(z*p1-p2)/(z*z-1);
      z = z-p1/pp;
    }
    while (std::abs(p1/pp) > tolerance);

    const double x = .5*z;
    qp[i-1] = .5-x;
    qp[n-i] = .5+x;

    const double w = 1./((1.-z*z)*pp*pp);
    qw[i-1] = w;
    qw[n-i] = w;
  }
}

// EOF
