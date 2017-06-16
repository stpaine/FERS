//xmlimport.h
//Import a simulator world and simulation parameters from an XML file
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//Started: 26 April 2006

#ifndef __XMLIMPORT_H
#define __XMLIMPORT_H

#include <config.h>
#include "rsworld.h"

namespace xml {
  //Load an XML file containing a simulation description
  //putting the structure of the world into the world object in the process
  void LoadXMLFile(std::string filename, rs::World *world);

  class XMLException {
  };
}

#endif
