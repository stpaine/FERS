//rsobject.cpp
//Implementation of Object from rsobject.h
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//19 July 2006

#include "rsobject.h"
#include "rsplatform.h"

using namespace rs;

Object::Object(const Platform *platform, std::string name):
  platform(platform),
  name(name)
{
}

Object::~Object() {
}

// Get the position of the object
Vec3 Object::GetPosition(rsFloat time) const
{
  return platform->GetPosition(time);
}

// Get the rotation of the object from it's platform
SVec3 Object::GetRotation(rsFloat time) const
{
  return platform->GetRotation(time);
}

// Get the object's name
std::string Object::GetName() const
{
  return name;
}

/// Get a pointer to the platform this object is attached to
const Platform* Object::GetPlatform() const
{
  return platform;
}
