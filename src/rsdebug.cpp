//Implementation of rsdebug.h - Debug messaging
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//20 March 2006
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdarg>
#include <boost/thread/mutex.hpp> //for boost::mutex
#include "rsdebug.h"

rsDebug::Level debug_level = rsDebug::RS_VERY_VERBOSE; //The current debug level 

//We use a mutex for log information printing, to stop messages from getting mangled
boost::mutex debugMutex;

//Print out a debug message including the line and file name
void rsDebug::print(const rsDebug::Level level, const std::string &str, const std::string file, const int line) {
  if (level >= debug_level) {
    boost::mutex::scoped_lock lock(debugMutex); //Lock the mutex
    std::ostringstream oss;
    oss << "[" << file << " " << line << "] ";
    oss << str << "\n";
    std::cerr << oss.str();
    //Mutex will automatically be unlocked here by scoped_lock
  }
}

//Formatted print of the current debug level, doesn't include filename and line
//Uses the cstdarg variable arguments system and the vfprintf function to handle the arguments
//If your system does not have the standard vfprintf function in it's library, you will have to make a pla
void rsDebug::printf(const rsDebug::Level level, const char *format, ...)
{
  if (level >= debug_level) {
    boost::mutex::scoped_lock lock(debugMutex); //Lock the mutex
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    //mutex will automatically be unlocked here by scoped_lock
  }
}

//See comments for printf(Level, char *)
void rsDebug::printf(const rsDebug::Level level, const std::string &format, ...){
  if (level >= debug_level) {
    boost::mutex::scoped_lock lock(debugMutex); //Lock the mutex
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format.c_str(), ap);
    va_end(ap);
    //mutex will automatically be unlocked here by scoped_lock
  }  
}

//Set the current debug level
void rsDebug::setDebugLevel(rsDebug::Level level) {
  if (level <= rsDebug::RS_EXTREMELY_CRITICAL)
    debug_level = level;
}


