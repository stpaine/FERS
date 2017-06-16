//rssim.h
//Declarations and definitions for rssim.cpp and rsthreadedsim.cpp

#ifndef __RSSIM_H
#define __RSSIM_H

#include <config.h>
#include "rsworld.h"
#include "rsradar.h"

namespace rs {

  //Functions in rssim.cpp

  //Run a simulation on the specified receiver/transmitter pair
  void SimulatePair(const Transmitter *trans, Receiver *recv, const World *world);

  //Functions in rsthreadedsim.cpp

  //Run the radar simulation specified by world
  //Limit the number of concurrent threads to thread_limit
  void RunThread(int thread_limit, World *world);

}

#endif
