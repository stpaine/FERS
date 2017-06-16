//rsgeometry.h
//Classes for calculating the geometry of the world
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//26 May 2006

#ifndef __RSGEOMETRY_H
#define __RSGEOMETRY_H

#include <config.h>

namespace rsGeometry {
}

namespace rs {

class SVec3;

/// 3x3 Matrix, with simple operations
class Matrix3 {
 public:
  rsFloat elements[9];
  /// Default constructor
  Matrix3();
  /// Default destructor
  ~Matrix3();
  /// Get a pointer to the element array
  const rsFloat *GetData() const;
  rsFloat *GetData();
};

/// The Vec3 class is a rectangular 3 vector
class Vec3 {
 public:
  rsFloat x, y, z;
  /// Default constructor
  Vec3();
  /// Constructor which sets co-ordinates
  Vec3(rsFloat x, rsFloat y, rsFloat z);
  /// Constructor with a spherical vector
  Vec3(const SVec3 &svec);
  /// Default destructor
  ~Vec3();
  //
  //operators
  //

  //Vector operations
  Vec3 &operator+=(const Vec3 &b); //!< Vector Addition
  Vec3 &operator-=(const Vec3 &b); //!< Vector Subtraction
  Vec3 &operator*=(const Vec3 &b); //!< Componentwise multiplication
  Vec3 &operator=(const Vec3 &b); //!< Assignment operator

  //Matrix operations
  Vec3 &operator*=(const Matrix3 &m); //!< Multiplication by a 3x3 matrix

  //Scalar operations
  Vec3 &operator*=(const rsFloat b); //!< Multiplication by a scalar
  Vec3 &operator/=(const rsFloat b); //!< Division by a scalar
  Vec3 &operator+=(const rsFloat b); //!< Addition of a scalar

  /// Return the length of the vector
  rsFloat Length() const;
};

//Vector operations
rsFloat DotProduct(const Vec3 &a, const Vec3 &b); //!< Vector Inner product
Vec3 CrossProduct(const Vec3 &a, const Vec3 &b); //!< Vector Cross product
//Componentwise vector operations
Vec3 operator*(const Vec3 &a, const Vec3 &b); //!< Componentwise product
Vec3 operator+(const Vec3 &a, const Vec3 &b); //!< Componentwise add
Vec3 operator-(const Vec3 &a, const Vec3 &b); //!< Componentwise subtract
Vec3 operator/(const Vec3 &a, const Vec3 &b); //!< Componentwise divide
//Scalar operations
Vec3 operator*(const Vec3 &a, const rsFloat b); //!< Multiply by a scalar
Vec3 operator/(const Vec3 &a, const rsFloat b); //!< Division by a scalar
Vec3 operator/(const rsFloat a, const Vec3 &b); //!< Division by a scalar

/// The SVec3 class is a vector in R^3, stored in spherical co-ordinates
class SVec3 {
 public:
  rsFloat length; //!< The length of the vector
  rsFloat azimuth; //!< Angle in the x-y plane (radians)
  rsFloat elevation; //!< Elevation angle above x-y plane (radians)
  /// Default constructor
  SVec3();
  /// Constructor with initializers
  SVec3(rsFloat length, rsFloat azimuth, rsFloat elevation);
  /// Copy constructor
  SVec3(const SVec3 &svec);
  /// Constructor with a rectangular vector
  SVec3(const Vec3 &vec);
  /// Destructor
  ~SVec3();
  //
  //operators
  //
  SVec3 &operator*=(rsFloat b); //!< multiplication by a scalar
  SVec3 &operator/=(rsFloat b); //!< division by a scalar
};

}
#endif
