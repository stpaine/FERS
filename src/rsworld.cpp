//rsworld.cpp
//Implementation of simulator world object
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//Started: 25 April 2006
#include <algorithm>
#include "rsdebug.h"
#include "rsworld.h"
#include "rsantenna.h"
#include "rsradar.h"
#include "rsantenna.h"
#include "rsradarwaveform.h"
#include "rsplatform.h"
#include "rstarget.h"
#include "rstiming.h"
#include "rsmultipath.h"

using namespace rs; //Import the rs namespace for this implementation

/// Function object used to delete objects from vector
  template <typename T> struct objDel
  {
    //operator() deletes the object that x is a pointer to
    void operator()(T x) {
      delete x;
    }
  };


//Default constructor
World::World():
  multipath_surface(0)
{
  
}

//World object destructor
World::~World()
{
  //Delete all the objects in the world 
  std::map<std::string, RadarSignal*>::iterator iter;
  for(iter = pulses.begin(); iter != pulses.end(); iter++)
    delete (*iter).second;

  std::map<std::string, Antenna*>::iterator aiter;
  for (aiter = antennas.begin(); aiter != antennas.end(); aiter++)
    delete (*aiter).second;

  std::map<std::string, PrototypeTiming*>::iterator titer;
  for (titer = timings.begin(); titer != timings.end(); titer++)
    delete (*titer).second;  

  std::for_each(receivers.begin(), receivers.end(), objDel<Receiver *>());
  std::for_each(transmitters.begin(), transmitters.end(), objDel<Transmitter *>());
  std::for_each(targets.begin(), targets.end(), objDel<Target *>());

  //Platforms are deleted last as they are referred to by the other object types
  std::for_each(platforms.begin(), platforms.end(), objDel<Platform *>());

}

//Add a platform to the world
void World::Add(Platform *platform)
{
  platforms.push_back(platform);
}

//Add a transmitter to the world
void World::Add(Transmitter *trans)
{
  transmitters.push_back(trans);
}

//Add a receiver to the world
void World::Add(Receiver *recv)
{
  receivers.push_back(recv);
}

//Add a target to the world
void World::Add(Target *targ)
{
  targets.push_back(targ);
}

//Add a pulse to the world
void World::Add(RadarSignal *pulse)
{
  if (FindSignal(pulse->GetName()))
    throw std::runtime_error("[ERROR] A pulse with the name "+pulse->GetName()+" already exists. Pulses must have unique names");
  pulses[pulse->GetName()] = pulse;
}

//Add an antenna to the world
void World::Add(Antenna *antenna)
{
  if (FindAntenna(antenna->GetName()))
    throw std::runtime_error("[ERROR] An antenna with the name "+antenna->GetName()+" already exists. Antennas must have unique names");
  antennas[antenna->GetName()] = antenna;
}

//Add a timing source to the world
void World::Add(PrototypeTiming *timing)
{
  if (FindTiming(timing->GetName()))
    throw std::runtime_error("[ERROR] A timing source with the name "+timing->GetName()+" already exists. Timing sources must have unique names");
  timings[timing->GetName()] = timing;
}

//Get a pulse from the map of pulses
RadarSignal* World::FindSignal(const std::string& name)
{
  return pulses[name];
}

//Get an antenna from the map of antennas
Antenna* World::FindAntenna(const std::string& name)
{
    return antennas[name];
}

/// Find a timing source with the specified name
PrototypeTiming* World::FindTiming(const std::string& name)
{
  return timings[name];
}

///Add a multipath surface to the world
void World::AddMultipathSurface(MultipathSurface *surface)
{
  if (multipath_surface)
    throw std::runtime_error("[ERROR] Only one multipath surface per simulation is supported");
   multipath_surface = surface;
}

///Process the scene to add virtual receivers and transmitters
void World::ProcessMultipath()
{
  // In this function "duals" are added for each transmitter and receiver
  // a dual has the same properties of the transmitter and receiver, but is reflected in the multipath plane
  if (multipath_surface) {
    //Add duals for each plaform
    std::vector<Platform*>::iterator plat = platforms.begin();
    std::vector<Platform*>::iterator plat_end = platforms.end();
    for (; plat != plat_end; plat++)
      platforms.push_back(CreateMultipathDual(*plat, multipath_surface));
    //Add duals for each receiver
    std::vector<Receiver*>::iterator recv = receivers.begin();
    std::vector<Receiver*>::iterator recv_end = receivers.end();
    for (; recv != recv_end; recv++)
      receivers.push_back(CreateMultipathDual(*recv, multipath_surface));
    //Add duals for each transmitter
    std::vector<Transmitter*>::iterator trans = transmitters.begin();
    std::vector<Transmitter*>::iterator trans_end = transmitters.end();
    for (; trans != trans_end; trans++)
      transmitters.push_back(CreateMultipathDual(*trans, multipath_surface));
  }
  //Clean up the multipath surface
  delete multipath_surface;
  multipath_surface = 0;
}
