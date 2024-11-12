#include "Runge_Kutta_Butcher.hpp"
#include <iostream>

// Constructor implementation
Runge_Kutta_Butcher::Runge_Kutta_Butcher(const int &steps, const int &order, const bool &flag) : ss(steps), mm(order), isImplicit(flag) 
{
  // Initialize the size of c and b
  cc.resize(ss);
  bb.resize(ss);

  // Initialize the size of aij
  aa.resize(ss, std::vector<double>(ss, 0.0));  // Set all aij to 0 by default

  // Set RK coefficient based on steps and order
  setCoefficients();
}

// Set coefficients based on the given number of steps and order
void Runge_Kutta_Butcher::setCoefficients() 
{
  switch (ss) 
  {
    case 3:
      if (mm == 3)
      {
        cc[0] = 0.0; cc[1] = 8.0/15.0; cc[2] = 2.0/3.0;
        aa[0][0] = 0.0; aa[1][0] = 8.0/15.0;
        aa[2][0] = 1.0/4.0; aa[2][1] = 5.0/12.0;
        bb[0] = 1.0/4.0; bb[1] = 0.0; bb[2] = 3.0/4.0;
      }
      break;
    case 5:
      if (mm == 3) 
      {
        const double gamma = isImplicit ? 0.57281606248208 : 0.435866521508482;
        cc[0] = 0.0; cc[1] = isImplicit ? 2*gamma : 2*gamma; cc[2] = isImplicit ? 0.491513338007346 : 0.871733043016825; cc[3] = isImplicit ? 0.659057937017090 : 2.340212575108680; cc[4] = isImplicit ? 1.0 : 1.0;
        bb[0] = isImplicit ? -0.053475117047939 : 0.412898042812474; bb[1] = isImplicit ? 0.0 : 0.0; bb[2] = isImplicit ? 2.325340415003076 : 0.0; bb[3] = isImplicit ? -1.844681360437218 : 0.197339890378869; bb[4] = isImplicit ? gamma : gamma;

        aa[0][0] = 0.0;
        aa[1][0] = isImplicit ? gamma: 2*gamma; aa[1][1] = isImplicit ? gamma : 0.0;
        aa[2][0] = isImplicit ? 0.059016362876503 : gamma; aa[2][1] = isImplicit ? -0.140319087351238 : gamma; aa[2][2] = isImplicit ? gamma : 0.0;
        aa[3][0] = isImplicit ? 0.111666256495189 : -0.800998453065629; aa[3][1] = isImplicit ?  -0.226019907027712 : 0.0; aa[3][2] = isImplicit ? 0.4 : 3.141211028174309; aa[3][3] = isImplicit ? gamma : 0.0;
        aa[4][0] = isImplicit ? -0.053475117047939 : 0.356753207779640; aa[4][1] = isImplicit ? 0.0 : -0.197339890378869; aa[4][2] = isImplicit ? 0.0 : 0.881948841393791; aa[4][3] = isImplicit ? 2.325340415003076 : -0.041362158794624; aa[4][4] = isImplicit ? gamma : 0.0;
            
        // Explicit-RK
        // gama = 424782/974569 = 0.435866521508482
        // c2 = 902905985686/1035759735069 = 0.871733043016825
        // c3 = 2684624/1147171 = 2.340212575108680
        // b0 = 487698502336740678603511/1181159636928185920260208 = 0.412898042812474
        // b3 = 302987763081184622639300143137943089/1535359944203293318639180129368156500 = 0.197339890378869
        // a30 = -475883375220285986033264/594112726933437845704163 = -0.800998453065629
        // a32 = 1866233449822026827708736/594112726933437845704163 = 3.141211028174309
        // a40 = 62828845818073169585635881686091391737610308247/176112910684412105319781630311686343715753056000
        // a41 = -302987763081184622639300143137943089/1535359944203293318639180129368156500 = -0.197339890378869
        // a42 = 262315887293043739337088563996093207/297427554730376353252081786906492000 = 0.881948841393791
        // a43 = -987618231894176581438124717087/23877337660202969319526901856000 = -0.041362158794624

        // Implicit-RK
        // gama = 2051948/3582211 = 0.57281606248208
        // c2 = 12015769930846/24446477850549 = 0.491513338007346
        // c3 = 3532944/5360597 = 0.659057937017090
        // b0 = -2032971420760927701493589/38017147656515384190997416 = -0.053475117047939
        // b2 = 2197602776651676983265261109643897073447/945067123279139583549933947379097184164 = 2.325340415003076
        // b3 = -128147215194260398070666826235339/69468482710687503388562952626424 = -1.844681360437218
        // a20 = 259252258169672523902708425780469319755/4392887760843243968922388674191715336228 = 0.059016362876503
        // a21 = -172074174703261986564706189586177/1226306256343706154920072735579148 = -0.140319087351238
        // a30 = 1103202061574553405285863729195740268785131739395559693754/9879457735937277070641522414590493459028264677925767305837 = 0.111666256495189
        // a31 = -103754520567058969566542556296087324094/459050363888246734833121482275319954529 =  -0.226019907027712
        // a32 = 3863207083069979654596872190377240608602701071947128/19258690251287609765240683320611425745736762681950551 = 0.200595525067531
      }
      break;
      // Other cases
    default:
      SYS_T::print_fatal("Error: Runge_Kutta_Butcher::setCoefficients: Unsupported number of steps or order. \n");
      break;
  }
}

// print coefficients
void Runge_Kutta_Butcher::printCoefficients() const
{
  SYS_T::commPrint("Coefficients:\n");
  SYS_T::commPrint("c: \n");
  for (const auto& ci : cc) SYS_T::commPrint("%e \n", ci);
  SYS_T::commPrint("\nb \n");
  for (const auto& bi : bb) SYS_T::commPrint("%e \n", bi);
  SYS_T::commPrint("\na:\n");
  for (const auto& row : aa) 
  {
    for (const auto& aij : row) SYS_T::commPrint("%e ", aij);
      SYS_T::commPrint("\n");
  }
}
