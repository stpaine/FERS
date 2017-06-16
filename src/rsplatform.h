//rsplatform.h
//Simulator Platform Object
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//Started: 21 April 2006

#ifndef __RSPLATFORM_H
#define __RSPLATFORM_H

#include <config.h>
#include "rspath.h"
#include <string>
#include <boost/utility.hpp>

namespace rs { 

  //Forward definition of MultipathSurface (rsmultipath.h)
  class MultipathSurface;
  
/// The Platform class controls the motion and rotation of all objects in the scene
class Platform: boost::noncopyable {
public:
  /// Default Constructor
  Platform(const std::string &name);
  /// Destructor
  ~Platform();
  /// Return a pointer to the motion path
  Path *GetMotionPath();
  /// Return a pointer to the rotation path
  RotationPath *GetRotationPath(); 
  /// Get the position of the platform at the specified time
  Vec3 GetPosition(rsFloat time) const;
  /// Get the rotation of the platform at the specified time
  SVec3 GetRotation(rsFloat time) const;
  /// Get the name of the platform
  std::string GetName() const;

private:
  Path *motionPath; //!< Pointer to platform's motion path
  RotationPath *rotationPath; //!< Pointer to platform's rotation path
  std::string name; //!< The name of the platform
  Platform *dual; //!< Multipath dual of this platform
  /// Create a dual of this platform for multipath simulation
  friend Platform* rs::CreateMultipathDual(const Platform *plat, const MultipathSurface *surf);
};

  /// Create a dual of this platform for multipath simulation
  Platform *CreateMultipathDual(const Platform *plat, const MultipathSurface *surf);

}

#endif
