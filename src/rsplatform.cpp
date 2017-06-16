//rsplatform.cpp
//Implementation of simulator platform class
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//Started: 26 April 2006

#include "rsplatform.h"
#include "rsmultipath.h"

using namespace rs;

//Default constructor
Platform::Platform(const std::string &name):
  name(name),
  dual(0)
{
  motionPath = new Path();
  rotationPath = new RotationPath();
}

//Default destructor
Platform::~Platform()
{
  delete motionPath;
  delete rotationPath;
}

//Get the current position of the platform
Vec3 Platform::GetPosition(rsFloat time) const
{
  return motionPath->GetPosition(time);
}

//Get the current rotation of the platform
SVec3 Platform::GetRotation(rsFloat time) const
{
  return rotationPath->GetPosition(time);
}

//Return a pointer to the motion path
Path* Platform::GetMotionPath()
{
  return motionPath;
}

//Return a pointer to the rotation path
RotationPath* Platform::GetRotationPath()
{
  return rotationPath;
}

//Return the name of the platform
std::string Platform::GetName() const
{
  return name;
}

/// Create a dual of this platform for multipath simulation
Platform *rs::CreateMultipathDual(const Platform *plat, const MultipathSurface *surf) {
  //If the dual already exists, just return it
  if (plat->dual)
    return plat->dual;
  //Create the new platform
  Platform *dual = new Platform(plat->GetName()+"_dual");
  //Set the platform dual
  (const_cast<Platform*>(plat))->dual = dual;
  //Reflect the paths used to guide the platform
  dual->motionPath = ReflectPath(plat->motionPath, surf);
  dual->rotationPath = ReflectPath(plat->rotationPath, surf);
  //Done, return the created object
  return dual;

}
