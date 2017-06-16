//rspython.h - Objects for extending FERS with Python
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//7 March 2007

#ifndef __RSPYTHON_H
#define __RSPYTHON_H

#include <config.h>
#include <string>
#include "rsgeometry.h"

namespace rsPython {

  ///Function which initializes python
  void InitPython();

  ///Class which implements Python functionality
  struct PythonExtensionData;

  ///Parent class for all python calculation classes
  class PythonExtension {
  public:
    /// Constructor
    PythonExtension(const std::string& module, const std::string& function);
    /// Destructor
    ~PythonExtension();
  protected:
    PythonExtensionData *data; //!< Pointer to python specific data
    std::string module; //!< Name of the module
    std::string function; //!< Name of the function
  };

  ///Python extension for calculating paths
  class PythonPath: public PythonExtension {
  public:
    /// Constructor
    PythonPath(const std::string& module, const std::string& function);
    /// Destructor
    ~PythonPath();
    /// Get the position at the given time
    rs::Vec3 GetPosition(rsFloat t) const;
  };

  /// Python extension for generating noise
  class PythonNoise: public PythonExtension {
  public:
    ///Constructor
    PythonNoise(const std::string& module, const std::string& function);
    ///Destructor
    ~PythonNoise();
    /// Get a noise sample
    rsFloat GetSample() const;
  };

  /// Python extension for generating noise
  class PythonAntennaMod: public PythonExtension {
  public:
    ///Constructor
    PythonAntennaMod(const std::string& module, const std::string& function);
    ///Destructor
    ~PythonAntennaMod();
    /// Get a noise sample
    rsFloat GetGain(const rs::SVec3& direction) const;
  };


}

#endif
