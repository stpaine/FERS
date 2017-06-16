//rsobject.h 
//Base class for all physical objects in the world
//Marc Brooker
//19 July 2006

#ifndef __RSOBJECT_H
#define __RSOBJECT_H

#include <config.h>
#include "rspath.h"
#include <string>
#include <boost/utility.hpp>

namespace rs {
  
  // Forward definition of Platform (see rsplatform.h)
  class Platform;

  /// Object is a physical object in the world, it embodies the properties common to all physical objects
  class Object: boost::noncopyable {
  public:
    Object(const Platform *platform, std::string name);
    ~Object();
    /// Get the position of the object in space
    Vec3 GetPosition(rsFloat time) const;
    /// Get the rotation of the object in space
    SVec3 GetRotation(rsFloat time) const;
    /// Get a pointer to the platform this object is attached to
    const Platform *GetPlatform() const;
    /// Get the name of the object
    std::string GetName() const;
  private:
    const Platform *platform; //<! Platform to which the object is attached
    std::string name; //<! Name of the object
  };
}

#endif
