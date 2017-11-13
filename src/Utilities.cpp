// Include classes
#include "Utilities.h"

#include <cmath>

using namespace std;

// Constructor
//=======================
Utilities::Utilities() {}
//=======================



//========================
Utilities::~Utilities() {}
//========================



//===========================================================================
float Utilities::DeltaR( float eta_1, float phi_1, float eta_2, float phi_2 )
{
   float deta = eta_1 - eta_2;
   float dphi = acos(cos(phi_1 - phi_2)); // To cover larger then PI angle between the paricles
//   float dphi = DeltaPhi(phi_1, phi_2);
   return sqrt(deta*deta + dphi*dphi);
}
//===========================================================================


//===================================================
float Utilities::DeltaPhi( float phi_1, float phi_2 )
{
   const double PI = 2.0*acos(0.);
   const double TWOPI = 2.0*PI;
   
   float PHI = fabs(phi_1 - phi_2);
   return (PHI <= PI) ? PHI : TWOPI-PHI;
}
//===================================================
