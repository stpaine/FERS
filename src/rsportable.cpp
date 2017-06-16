//rsportable.cpp
//All code which is non-standard C++ must go in here, with reasons why
//There should be very little code in this file
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//26 July 2006

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include "rsportable.h"
#include "rsdebug.h"
#include <boost/thread.hpp>

/// Compare two strings, ignoring case
//There isn't a suitable one in the standard library of C or C++
int rsPortable::stricmp(const char *one, const char *two)
{
  return strcasecmp(one, two); //strcasecmp is a GNU extension
}

/// Compute the first order Bessel function of the first kind
rsFloat rsPortable::BesselJ1(rsFloat x)
{
  return j1(x); //j1 is non standard, but found on many platforms
}

/// Round off a floating point number
//This function isn't in C++03, but is in C++0x and TR1 and C99, so should work on most machines
rsFloat rsPortable::rsRound(rsFloat x)
{
  return round(x);
}

/// Detect the number of CPUs in the machine
int rsPortable::CountProcessors()
{
  //C Tong: This can now be done with boost:

  int iNHardwareThreads = boost::thread::hardware_concurrency();

  if(!iNHardwareThreads)
  {
    rsDebug::printf(rsDebug::RS_IMPORTANT, "[IMPORTANT] Unable to get CPU count, assuming 1.\n");
    return 1;
  }
  else
    return iNHardwareThreads;
}

