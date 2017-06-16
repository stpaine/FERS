//rsdebug.h - Message support functions and debug levels
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//20 March 2006
#ifndef __RS_DEBUG_H
#define __RS_DEBUG_H

#include <string>

//Macro which calls the debug print function with the current file and line
#define DEBUG_PRINT(level, str) rsDebug::print(level, str, __FILE__, __LINE__)

namespace rsDebug {
  

  //RS_INFORMATIVE - Messages which may be informative to the user
  //RS_IMPORTANT
  //RS_CRITICAL - Stuff that is very important
  //RS_EXTREMELY_CRITICAL - Very important messages that must be printed
  enum Level {RS_VERY_VERBOSE, //!< Messages which are only useful for debugging
	      RS_VERBOSE, //!< Messages which are unlikely to prove important
	      RS_INFORMATIVE, //!< Messages which may be informative to the user
	      RS_IMPORTANT, //!< Important messages
	      RS_CRITICAL, //!< Critical messages, such as errors which may lead to incorrect results
	      RS_EXTREMELY_CRITICAL //!< Extremely important messages
  };

  /// Print the current debug message, file and line, as long as level is greater or equal to the current debug level
  void print(const Level level, const std::string &str, const std::string file, const int line);

  /// Print a formatted debug message at the current level
  //Formatting is as per the C printf function, with all the format specifiers supported
  void printf(const Level level, const char *format, ...);

  /// Overloaded printf which takes format string as std::string
  void printf(const Level level, const std::string &format, ...);

  /// Set the current debug level
  void setDebugLevel(Level level);
}

#endif
