//rsmultipath.cpp
//Implementation of multipath propagation
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//9 September 2007

#include "rsmultipath.h"
#include "rsdebug.h"

using namespace rs;

//
// MultipathSurface implementation
//

/// Constructor
MultipathSurface::MultipathSurface(rsFloat a, rsFloat b, rsFloat c, rsFloat d, rsFloat factor):
  factor(factor)
{
  rsFloat *mat = reflection.GetData();
  //Fill the reflection matrix
  mat[0] = -a*a+b*b+c*c;
  mat[4] = a*a-b*b+c*c;
  mat[8] = a*a+b*b-c*c;
  mat[1] = mat[3] = -2*a*b;
  mat[2] = mat[6] = -2*a*c;
  mat[5] = mat[7] = -2*b*c;
  //Get the scale factor
  norm_factor = 1/(a*a+b*b+c*c);
  //Get the translation vector
  translation_vector = Vec3(-2*a*d, -2*b*d, -2*c*d);
}

/// Default destructor
MultipathSurface::~MultipathSurface()
{
}

/// Return a point reflected in the surface
Vec3 MultipathSurface::ReflectPoint(const Vec3 &b) const
{
  Vec3 ans = b;
  // Calculate the reflected position of b
  // Calculation is norm_factor*(reflection*b+translation_vector)
  ans *= reflection;
  ans -= translation_vector;
  ans *= norm_factor;
  return ans;
}

/// Get the reflectance factor
rsFloat MultipathSurface::GetFactor() const
{
  return factor;
}
  
