//rsthreadedsim.h
//Definitions for threaded simulation code
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//19 July 2006

#ifndef __RSTHREADEDSIM_H
#define __RSTHREADEDSIM_H

#include "rsworld.h"

namespace rs {
  void RunThreadedSim(int thread_limit, World *world);
}

#endif
