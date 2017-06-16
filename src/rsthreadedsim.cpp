//rsthreadedsim.cpp
//Thread management for the simulator
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//29 May 2006

//One of the goals for FERS is to support multiple processors. This is achieved
//through multithreading. One simulation is performed in for each transmitter-
//receiver pair. A number of these simulations are run in parallel by multi-
//threading, according to the number of CPUs (or cores) the system has.

#include <vector>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/xtime.hpp>
#include <boost/version.hpp>

#if BOOST_VERSION < 105000
#define TIME_UTC_ TIME_UTC
#endif

#include "rssim.h"
#include "rsworld.h"
#include "rsradar.h"
#include "rsdebug.h"
#include "rsthreadedsim.h"

using namespace rs;

//Counter of currently running threads
volatile int threads;
boost::mutex threadsMutex; //Mutex to protect it

/// Flag to set if a thread encounters an error
volatile int error;
boost::mutex errorMutex;

/// Method to decrease the running thread count
static void decThreads()
{
  boost::mutex::scoped_lock lock(threadsMutex);
  threads--;
}

/// Method to flag if a thread experienced an error
static void setError()
{
  boost::mutex::scoped_lock lock(errorMutex);
  error = 1;
}

//SimThread contains a simulator thread
class SimThread {
public:
  //Constructor
  SimThread(const Transmitter *transmitter, Receiver *receiver, const World *world):
    trans(transmitter), recv(receiver), world(world)
  {
  }

  //Operator () is executed when we create the thread
  void operator()() {
    rsDebug::printf(rsDebug::RS_VERBOSE, "[VERBOSE] Created simulator thread for transmitter '%s' and receiver '%s' ", trans->GetName().c_str(), recv->GetName().c_str());
    try {
      SimulatePair(trans, recv, world);
    }
    catch (std::exception &ex) {
      rsDebug::printf(rsDebug::RS_CRITICAL, "[ERROR] First pass thread terminated with unexpected error:\n\t%s\nSimulator will terminate\n", ex.what());
      setError();
    }
    decThreads();
  }

protected:
  //The transmitter/receiver pair to simulate
  const Transmitter *trans;
  Receiver *recv;
  const World *world;
  
  
};

/// RenderThread contains a thread which performs the second phase of the simulation
class RenderThread {
public:
  /// Constructor
  RenderThread(Receiver *recv):
    recv(recv)
  {
  }

  /// Operator () is executed when we create the thread
  void operator()() {
    rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "[VV] Created render thread for receiver '%s'\n", recv->GetName().c_str());
    try {
      recv->Render();
    }
    catch (std::exception &ex) {
      rsDebug::printf(rsDebug::RS_CRITICAL, "[ERROR] Render thread terminated with unexpected error:\n\t%s\nSimulator will terminate\n", ex.what());
      setError();
    }
    decThreads();
  }
protected:
  Receiver *recv; //!< The Receiver to render
  
};


/// Sleep for the specified number of seconds
static void Sleep(int secs)
{
  //We sleep for one second here
  boost::xtime xt;
  boost::xtime_get(&xt, boost::TIME_UTC_);
  xt.sec += secs;
  boost::thread::sleep(xt);
}

//Increase the count of running threads
static void IncThreads()
{
  boost::mutex::scoped_lock lock(threadsMutex);
  threads++;
}

//Run a sim thread for each of the receiver-transmitter pairs, limiting concurrent threads
void rs::RunThreadedSim(int thread_limit, World *world) {
  std::vector<boost::thread *> running;
  std::vector<Receiver*>::iterator ri;
  std::vector<Transmitter*>::const_iterator ti;
  boost::thread mainthrd();
  rsDebug::printf(rsDebug::RS_INFORMATIVE, "[INFO] Using threaded simulation with %d threads.\n", thread_limit);
  //PHASE 1: Do first pass of simulator
  //Loop through the lists for transmitters and receivers
  for (ri = world->receivers.begin(); ri != world->receivers.end(); ri++) {
    for (ti = world->transmitters.begin(); ti != world->transmitters.end(); ti++)
      {
	IncThreads();
	SimThread sim(*ti, *ri, world);
	boost::thread *thrd = new boost::thread(sim);	
	//Delay until a thread is terminated, if we have reached the limit
	while (threads >= thread_limit)
	  boost::thread::yield();
	//If a thread ended in error, abort the simulation
	if (error)
	  throw std::runtime_error("Thread terminated with error. Aborting simulation");
	//Add the thread pointers to a vector to be freed later
	running.push_back(thrd);
      }

  }
  //Wait for all the first pass threads to finish before continuing
  while (threads)
    {
      boost::thread::yield();
      //      rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "[VV] Main Thread Poll, Waiting on %d first pass threads.\n", threads);
      if (error)
	throw std::runtime_error("Thread terminated with error. Aborting simulation");
    }
  //Clean all the thread pointers
  for (std::vector<boost::thread*>::iterator i = running.begin(); i != running.end(); i++)
    delete *i;
  running.clear();

  // Report on the number of responses added to each receiver
  for (ri = world->receivers.begin(); ri != world->receivers.end(); ri++)
    rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "[VV] %d responses added to receiver '%s'\n", (*ri)->CountResponses(), (*ri)->GetName().c_str());

  //PHASE 2: Do render pass of simulation
  //Loop through the lists of receivers and set each to render
  for (ri = world->receivers.begin(); ri != world->receivers.end(); ri++)
      {
	IncThreads();
	RenderThread sim(*ri);
	boost::thread *thrd = new boost::thread(sim);	
	//Delay until a thread is terminated, if we have reached the limit
	while (threads >= thread_limit)
	  boost::thread::yield();
	//If a thread ended in error, abort the simulation
	if (error)
	  throw std::runtime_error("Thread terminated with error. Aborting simulation");
	//Add the thread pointers to a vector to be freed later
	running.push_back(thrd);
      }

  //Wait for all the render threads to finish
  while (threads)
    {
      boost::thread::yield();
      //rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "[VV] Main Thread Poll, Waiting on %d render threads.\n", threads);
      if (error)
	throw std::runtime_error("Thread terminated with error. Aborting simulation");
    }
  //Clean all the thread pointers
  for (std::vector<boost::thread*>::iterator i = running.begin(); i != running.end(); i++)
    delete *i;
  running.clear();
}
