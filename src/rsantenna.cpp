//rsantenna.cpp
//Implementation of Antenna Class
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//20 July 2006

#define TIXML_USE_STL //Tell tinyxml to use the STL instead of it's own string class

#include <stdexcept>
#include <algorithm>
#include "rsantenna.h"
#include "rsdebug.h"
#include <cmath>
#include "rsportable.h"
#include "rsinterp.h"
#include <tinyxml.h>
#include "rspattern.h"
#include "rsradarwaveform.h"

using namespace rs;
using namespace rsAntenna;

//One of the xml utility functions from xmlimport.cpp
rsFloat GetNodeFloat(TiXmlHandle &node);

namespace {
  //Return sin(x)/x
  rsFloat sinc(rsFloat theta)
  {
    return std::sin(theta)/(theta+std::numeric_limits<rsFloat>::epsilon());
  }

  //Return the first order, first kind bessel function of x, divided by x
  rsFloat j1c(rsFloat x)
  {
    if (x == 0)
      return 1;
    return rsPortable::BesselJ1(x)/(x);
  }
}

//Default constructor for the antenna
Antenna::Antenna(const std::string& name):
  lossFactor(1), //Antenna efficiency default is unity
  name(name)
{
}

//Antenna destructor
Antenna::~Antenna()
{
}

//Set the efficiency factor of the antenna
void Antenna::SetEfficiencyFactor(rsFloat loss)
{
  if (loss > 1)
    rsDebug::printf(rsDebug::RS_IMPORTANT, "Using greater than unity antenna efficiency, results might be inconsistent with reality.\n");
  lossFactor = loss;
}

//Get the efficiency factor
rsFloat Antenna::GetEfficiencyFactor() const
{
  return lossFactor;
}

//Return the name of the antenna
std::string Antenna::GetName() const
{
  return name;
}

// Get the angle (in radians) off boresight
rsFloat Antenna::GetAngle(const SVec3 &angle, const SVec3 &refangle) const
{
  //Get the angle off boresight
  SVec3 normangle(angle);
  normangle.length = 1;
  Vec3 cangle(normangle);
  Vec3 ref(refangle);
  return std::acos(DotProduct(cangle, ref));
}

/// Get the noise temperature of the antenna in a particular direction
rsFloat Antenna::GetNoiseTemperature(const SVec3 &angle) const
{
  return 0; //TODO: Antenna noise temperature calculation
}

//
//Isotropic Implementation
//

//Default constructor
Isotropic::Isotropic(const std::string& name):
  Antenna(name)
{
}

//Default destructor
Isotropic::~Isotropic()
{
}

//Return the gain of the antenna
rsFloat Isotropic::GetGain(const SVec3 &angle, const SVec3 &refangle, rsFloat wavelength) const
{
  return GetEfficiencyFactor();
}

//
// Gaussian Implementation
//

  
/// Constructor
Gaussian::Gaussian(const std::string& name, rsFloat azscale, rsFloat elscale):
  Antenna(name),
  azscale(azscale),
  elscale(elscale)
{
  
}

/// Destructor
Gaussian::~Gaussian()
{
}

/// Get the gain at an angle
rsFloat Gaussian::GetGain(const rs::SVec3 &angle, const rs::SVec3 &refangle, rsFloat wavelength) const
{
  SVec3 a = angle - refangle;
  rsFloat azfactor = std::exp(-a.azimuth*a.azimuth*azscale);
  rsFloat elfactor = std::exp(-a.elevation*a.elevation*elscale);
  return azfactor*elfactor;
}


//
//Sinc Implemetation
//

//Constructor
Sinc::Sinc(const std::string& name, rsFloat alpha, rsFloat beta, rsFloat gamma):
  Antenna(name),
  alpha(alpha),
  beta(beta),
  gamma(gamma)
{
}

//Default destructor
Sinc::~Sinc()
{
}

// Return the gain of the antenna at an angle
rsFloat Sinc::GetGain(const SVec3 &angle, const SVec3 &refangle, rsFloat wavelength) const
{
  //Get the angle off boresight
  rsFloat theta = GetAngle(angle, refangle);

  //FIX 2015/02/18 CA Tong:
  //std::pow<double> returns NaN for negative bases with certain fractional indices as they create an uneven 
  //root which is, of course, a complex number. Inputs therefore need to be cast to complex numbers before the 
  //calculation as this will return complex number. Then return the magnitude as the beam gain.

  rs::rsComplex complexSinc(::sinc(beta*theta), 0.0);
  rs::rsComplex complexGamma(gamma, 0.0);
 
  //See "Sinc Pattern" in doc/equations.tex for equation used here
  rsComplex complexGain = alpha * std::pow(complexSinc, complexGamma) * GetEfficiencyFactor();

  //rsDebug::printf(rsDebug::RS_IMPORTANT, "Theta = %f Gain = %f, %f\n", theta, complexGain.real(), complexGain.imag());

  return std::abs(complexGain);
}

//
// SquareHorn Implementation
//

//Constructor
SquareHorn::SquareHorn(const std::string& name, rsFloat dimension):
  Antenna(name),
  dimension(dimension)
{
}

//Destructor
SquareHorn::~SquareHorn()
{
}

//Return the gain of the antenna
//See doc/equations.tex for details
rsFloat SquareHorn::GetGain(const SVec3 &angle, const SVec3 &refangle, rsFloat wavelength) const
{
  rsFloat Ge = 4*M_PI*dimension*dimension/(wavelength*wavelength);
  rsFloat x = M_PI*dimension*std::sin(GetAngle(angle, refangle))/wavelength;
  rsFloat gain = Ge*std::pow(::sinc(x), 2);
  return gain*GetEfficiencyFactor();
}


//
// Parabolic Dish Antenna
//

// Constructor
ParabolicReflector::ParabolicReflector(const std::string& name, rsFloat diameter):
  Antenna(name),
  diameter(diameter)
{
}

//Destructor
ParabolicReflector::~ParabolicReflector()
{
}

//Return the gain of the antenna
//See doc/equations.tex for details
rsFloat ParabolicReflector::GetGain(const SVec3 &angle, const SVec3 &refangle, rsFloat wavelength) const
{
  rsFloat Ge = std::pow(M_PI*diameter/wavelength, 2);
  rsFloat x = M_PI*diameter*std::sin(GetAngle(angle, refangle))/wavelength;
  rsFloat gain = Ge*std::pow(2*::j1c(x), 2);
  return gain*GetEfficiencyFactor();
}

//
// FileAntenna implementation
//

/// Constructor
FileAntenna::FileAntenna(const std::string& name, const std::string &filename):
  Antenna(name)
{
  pattern = new Pattern(filename);
}

/// Default destructor
FileAntenna::~FileAntenna()
{
  delete pattern;
}

/// Get the gain at an angle
rsFloat FileAntenna::GetGain(const rs::SVec3 &angle, const rs::SVec3 &refangle, rsFloat wavelength) const
{
  SVec3 a1 = angle;
  SVec3 a2 = refangle;
  SVec3 in_angle = (a1-a2);
  //  rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "az: %g el: %g\t az2: %g el2: %g\t az3: %g el3: %g\n", angle.azimuth, angle.elevation, refangle.azimuth, refangle.elevation, in_angle.azimuth, refangle.elevation);
  return pattern->GetGain(in_angle)*GetEfficiencyFactor();
}

//
// Antenna with pattern loaded from an XML file
//

// Constructor
XMLAntenna::XMLAntenna(const std::string& name, const std::string &filename):
  Antenna(name)
{
  // Classes to interpolate across elevation and azimuth
  azi_samples = new InterpSet();
  elev_samples = new InterpSet();
  //Load the XML antenna description data
  LoadAntennaDescription(filename);
}

//Destructor
XMLAntenna::~XMLAntenna()
{
  // Clean up the interpolation classes
  delete azi_samples;
  delete elev_samples;
}

//Return the gain of the antenna
//See doc/equations.tex for details
rsFloat XMLAntenna::GetGain(const SVec3 &angle, const SVec3 &refangle, rsFloat wavelength) const
{
  SVec3 t_angle = angle-refangle;
  rsFloat azi_gain = azi_samples->Value(std::fabs(t_angle.azimuth));
  rsFloat elev_gain = elev_samples->Value(std::fabs(t_angle.elevation));
  return azi_gain*elev_gain*max_gain*GetEfficiencyFactor();
}

namespace {

//Load samples of gain along an axis (not a member of XMLAntenna)
void LoadAntennaGainAxis(InterpSet *set, TiXmlHandle &axisXML)
{
  rsFloat angle;
  rsFloat gain;
  //Step through the XML file and load all the gain samples
  TiXmlHandle tmp = axisXML.ChildElement("gainsample", 0);
  for (int i = 0; tmp.Element() != 0; i++) {
    //Load the angle of the gain sample
    TiXmlHandle angleXML = tmp.ChildElement("angle", 0);
    if (!angleXML.Element())
      throw std::runtime_error("[ERROR] Misformed XML in antenna description: No angle in gainsample");
    angle = GetNodeFloat(angleXML);
    //Load the gain of the gain sample
    TiXmlHandle gainXML = tmp.ChildElement("gain", 0);
    if (!gainXML.Element())
      throw std::runtime_error("[ERROR] Misformed XML in antenna description: No gain in gainsample");
    gain = GetNodeFloat(gainXML);
    //Load the values into the interpolation table
    set->InsertSample(angle, gain);
    //Get the next gainsample in the file
    tmp = axisXML.ChildElement("gainsample", i);
  }
}

}

//Load the antenna description file
void XMLAntenna::LoadAntennaDescription(const std::string& filename)
{
  TiXmlDocument doc(filename.c_str());
  //Check the document was loaded correctly
  if (!doc.LoadFile())
    throw std::runtime_error("[ERROR] Could not load antenna description "+filename);
  //Get the XML root node
  TiXmlHandle root(doc.RootElement());
  //Load the gain samples along the elevation axis
  TiXmlHandle tmp = root.ChildElement("elevation", 0);
  if (!tmp.Element())
    throw std::runtime_error("[ERROR] Malformed XML in antenna description: No elevation pattern definition");
  LoadAntennaGainAxis(elev_samples, tmp);
  //Load the gain samples along the azimuth axis
  tmp = root.ChildElement("azimuth", 0);
  if (!tmp.Element())
    throw std::runtime_error("[ERROR] Malformed XML in antenna description: No azimuth pattern definition");
  LoadAntennaGainAxis(azi_samples, tmp);
  // Normalize the antenna patterns and calculate the max gain
  max_gain = std::max(azi_samples->Max(), elev_samples->Max());
  elev_samples->Divide(max_gain);
  azi_samples->Divide(max_gain);
  
}

//
// Antenna with gain pattern calculated by a Python program
//

// Constructor
PythonAntenna::PythonAntenna(const std::string& name, const std::string &module, const std::string& function):
  Antenna(name),
  py_antenna(module, function)
{
  
}

//Destructor
PythonAntenna::~PythonAntenna()
{
}

//Return the gain of the antenna
rsFloat PythonAntenna::GetGain(const SVec3 &angle, const SVec3 &refangle, rsFloat wavelength) const
{
  SVec3 angle_bore = angle - refangle; //Calculate the angle off boresight
  rsFloat gain = py_antenna.GetGain(angle_bore);
  return gain*GetEfficiencyFactor();
}


//
// Functions to create Antenna objects with a variety of properties
//

//Create an isotropic antenna with the specified name
Antenna* rs::CreateIsotropicAntenna(const std::string &name)
{
  rsAntenna::Isotropic *iso = new rsAntenna::Isotropic(name);
  return iso;
}

//Create a Sinc pattern antenna with the specified name, alpha and beta
Antenna* rs::CreateSincAntenna(const std::string &name, rsFloat alpha, rsFloat beta, rsFloat gamma)
{
  rsAntenna::Sinc *sinc = new rsAntenna::Sinc(name, alpha, beta, gamma);
  return sinc;
}

//Create a Gaussian pattern antenna
Antenna* rs::CreateGaussianAntenna(const std::string &name, rsFloat azscale, rsFloat elscale)
{
  rsAntenna::Gaussian *gau = new rsAntenna::Gaussian(name, azscale, elscale);
  return gau;
}

//Create a square horn antenna
Antenna* rs::CreateHornAntenna(const std::string &name, rsFloat dimension)
{
  rsAntenna::SquareHorn *sq = new rsAntenna::SquareHorn(name, dimension);
  return sq;
}

//Create a parabolic reflector antenna
Antenna* rs::CreateParabolicAntenna(const std::string &name, rsFloat diameter)
{
  rsAntenna::ParabolicReflector *pd = new rsAntenna::ParabolicReflector(name, diameter);
  return pd;
}

//Create an antenna with it's gain pattern stored in an XML file
Antenna* rs::CreateXMLAntenna(const std::string &name, const std::string &file)
{
  rsAntenna::XMLAntenna *fa = new rsAntenna::XMLAntenna(name, file);
  return fa;
}

//Create an antenna with it's gain pattern stored in an XML file
Antenna* rs::CreateFileAntenna(const std::string &name, const std::string &file)
{
  rsAntenna::FileAntenna *fa = new rsAntenna::FileAntenna(name, file);
  return fa;
}

//Create an antenna with gain pattern described by a Python program
Antenna* rs::CreatePythonAntenna(const std::string &name, const std::string &module, const std::string &function)
{
  rsAntenna::PythonAntenna* pa = new rsAntenna::PythonAntenna(name, module, function);
  return pa;
}
