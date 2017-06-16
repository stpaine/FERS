//rsantenna.h
//Class for Antennas, with different gain patterns, etc
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//20 July 2006

#ifndef __RSANTENNA_H
#define __RSANTENNA_H

#include <config.h>
#include <string>
#include "rsgeometry.h"
#include "boost/utility.hpp"
#include "rspython.h"
#include <limits>

namespace rs {
  //Forward declaration of SVec3 and Vec3 (see rspath.h)
  class SVec3;
  class Vec3;
  //Forward declaration of InterpSet (see rsinterp.h)
  class InterpSet;
  //Forward declaration of Pattern (see rspattern.h)
  class Pattern;

  /// The antenna class defines an antenna, which may be used by one or more transmitters
  class Antenna: public boost::noncopyable { //Antennas are not meant to be copied
  public:
    /// Default constructor
    Antenna(const std::string& name);
    /// Destructor
    virtual ~Antenna();
    /// Returns the current gain at a particular angle
    virtual rsFloat GetGain(const SVec3& angle, const SVec3& refangle, rsFloat wavelength) const = 0;
    /// Returns the noise temperature at a particular angle
    virtual rsFloat GetNoiseTemperature(const SVec3& angle) const;
    /// Set the antenna's loss factor (values > 1 are physically impossible)
    void SetEfficiencyFactor(rsFloat loss);
    /// Gets the Loss Factor
    rsFloat GetEfficiencyFactor() const;
    /// Return the name of the antenna
    std::string GetName() const;
  protected:
    /// Get the angle off boresight
    rsFloat GetAngle(const SVec3 &angle, const SVec3 &refangle) const;
  private:
    /// The loss factor for this antenna
    rsFloat lossFactor; //!< Loss factor
    /// The name of the antenna
    std::string name;
  };

  // Functions to create Antenna objects

  /// Create an Isotropic Antenna
  Antenna* CreateIsotropicAntenna(const std::string &name);

  /// Create an antenna with it's gain pattern stored in an XML file
  Antenna* CreateXMLAntenna(const std::string &name, const std::string &file);

  /// Create an antenna with it's gain pattern stored in an HDF5 file
  Antenna* CreateFileAntenna(const std::string &name, const std::string &file);

  /// Create an antenna with gain pattern described by a Python program
  Antenna* CreatePythonAntenna(const std::string &name, const std::string &module, const std::string &function);

  /// Create a Sinc Pattern Antenna
  // see rsantenna.cpp for meaning of alpha and beta
  Antenna* CreateSincAntenna(const std::string &name, rsFloat alpha, rsFloat beta, rsFloat gamma);

  /// Create a Gaussian Pattern Antenna
  Antenna* CreateGaussianAntenna(const std::string &name, rsFloat azscale, rsFloat elscale);

  /// Create a Square Horn Antenna
  Antenna* CreateHornAntenna(const std::string &name, rsFloat dimension);
  
  /// Create a parabolic reflector dish
  Antenna* CreateParabolicAntenna(const std::string &name, rsFloat diameter);

}

namespace rsAntenna {

  //Antenna with an Isotropic radiation pattern
  class Isotropic: public rs::Antenna {
  public:
    /// Default constructor
    Isotropic(const std::string& name);
    /// Default destructor
    virtual ~Isotropic();
    /// Get the gain at an angle
    virtual rsFloat GetGain(const rs::SVec3 &angle, const rs::SVec3 &refangle, rsFloat wavelength) const;
  };

  //Antenna with a sinc (sinx/x) radiation pattern
  class Sinc: public rs::Antenna {
  public:
    /// Constructor
    Sinc(const std::string& name, rsFloat alpha, rsFloat beta, rsFloat gamma);
    /// Destructor
    virtual ~Sinc();
    /// Get the gain at an angle
    virtual rsFloat GetGain(const rs::SVec3 &angle, const rs::SVec3 &refangle, rsFloat wavelength) const;
  private:
    rsFloat alpha; //!< First parameter (see equations.tex)
    rsFloat beta; //!< Second parameter (see equations.tex)
    rsFloat gamma; //!< Third parameter (see equations.tex)
  };

  //Antenna with a Gaussian radiation pattern
  class Gaussian: public rs::Antenna {
  public:
    /// Constructor
    Gaussian(const std::string& name, rsFloat azscale, rsFloat elscale);
    /// Destructor
    virtual ~Gaussian();
    /// Get the gain at an angle
    virtual rsFloat GetGain(const rs::SVec3 &angle, const rs::SVec3 &refangle, rsFloat wavelength) const;
  private:
    rsFloat azscale; //!< Azimuth scale parameter 
    rsFloat elscale; //!< Elevation scale parameter 
  };

  /// Square horn antenna
  class SquareHorn: public rs::Antenna {
  public:
    /// Constructor
    SquareHorn(const std::string& name, rsFloat dimension);
    /// Default destructor
    ~SquareHorn();
    /// Get the gain at an angle
    rsFloat GetGain(const rs::SVec3 &angle, const rs::SVec3 &refangle, rsFloat wavelength) const;
  private:
    rsFloat dimension; //!< The linear size of the horn
  };

  /// Parabolic dish antenna
  class ParabolicReflector: public rs::Antenna {
  public:
    /// Constructor
    ParabolicReflector(const std::string& name, rsFloat diameter);
    /// Default destructor
    ~ParabolicReflector();
    /// Get the gain at an angle
    rsFloat GetGain(const rs::SVec3 &angle, const rs::SVec3 &refangle, rsFloat wavelength) const;
  private:
    rsFloat diameter;    
  };

  /// Antenna with gain pattern loaded from and XML description file
  class XMLAntenna: public rs::Antenna {
  public:
    /// Constructor
    XMLAntenna(const std::string& name, const std::string &filename);
    /// Default destructor
    ~XMLAntenna();
    /// Get the gain at an angle
    rsFloat GetGain(const rs::SVec3 &angle, const rs::SVec3 &refangle, rsFloat wavelength) const;
  private:
    /// Load data from the antenna description file
    void LoadAntennaDescription(const std::string& filename);
    rsFloat max_gain; //!< Maximum Antenna gain
    rs::InterpSet* azi_samples; //!< Samples in the azimuth direction
    rs::InterpSet* elev_samples; //!< Samples in the elevation direction
  };

  /// Antenna with gain pattern loaded from an HDF5 2D pattern (as made by antennatool)
  class FileAntenna: public rs::Antenna {
  public:
    /// Constructor
    FileAntenna(const std::string& name, const std::string &filename);
    /// Default destructor
    ~FileAntenna();
    /// Get the gain at an angle
    rsFloat GetGain(const rs::SVec3 &angle, const rs::SVec3 &refangle, rsFloat wavelength) const;
  private:
    /// The antenna gain pattern
    rs::Pattern *pattern;
  };

  /// Antenna with gain pattern calculated by a Python module
  class PythonAntenna: public rs::Antenna {
  public:
    /// Constructor
    PythonAntenna(const std::string& name, const std::string &module, const std::string& function);
    /// Default destructor
    ~PythonAntenna();
    /// Get the gain at an angle
    rsFloat GetGain(const rs::SVec3 &angle, const rs::SVec3 &refangle, rsFloat wavelength) const;
  private:
    rsPython::PythonAntennaMod py_antenna;
  };

}

#endif
