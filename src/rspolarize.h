// Polarization support classes
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//31 March 2008

#ifndef __RS_POLARIZE_H
#define __RS_POLARIZE_H

#include "config.h"
#include <complex>

namespace rs {

///Polarization Scattering Matrix, implements matrices in the Jones calculus
class PSMatrix
{
 public:
  ///Default constructor creates identity PSM
  PSMatrix();
  /// Constructor
  PSMatrix(rsFloat s11, rsFloat s12, rsFloat s21, rsFloat s22);
  /// Copy constructor
  PSMatrix(const PSMatrix &im);
  /// Assignment operator
  PSMatrix& operator= (const PSMatrix &im);
  /// Matrix values
  std::complex<rsFloat> s[4];
};

//Jones vector class
class JonesVector
{
 public:
  /// Constructor
  JonesVector(std::complex<rsFloat> h, std::complex<rsFloat> v);
  /// Copy constructor
  JonesVector(const JonesVector &iv);
  /// Assignment operator
  JonesVector& operator= (const JonesVector &iv);
  /// Multiplication operator
  JonesVector operator* (const PSMatrix &mat);
  /// The horizontal polarization part
  std::complex<rsFloat> h;
  /// The vertical polarization part
  std::complex<rsFloat> v;
};

//Support functions

/// Dot product of two Jones vectors
std::complex<rsFloat> dot(const JonesVector &a, const JonesVector &b);

}

#endif
