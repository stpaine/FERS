//rstarget.cpp - Classes for targets and target RCS
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//11 June 2007

#define TIXML_USE_STL //Tell tinyxml to use the STL instead of it's own string class

#include <cmath>
#include "rstarget.h"
#include "rsdebug.h"
#include "rsinterp.h"
#include "rsnoise.h"
#include <tinyxml.h>

using namespace rs;

//One of the xml utility functions from xmlimport.cpp
rsFloat GetNodeFloat(TiXmlHandle &node);

//
// RCSModel Implementation
//

/// Destructor
RCSModel::~RCSModel()
{
}

//
// RCSConst Implementation
//

/// Return a constant RCS
rsFloat RCSConst::SampleModel() 
{
  return 1.0;
}

/// Destructor
RCSConst::~RCSConst()
{
}

//
// RCSChiSquare Implementation

/// Constructor
RCSChiSquare::RCSChiSquare(rsFloat k)
{
  gen = new GammaGenerator(k);
}
/// Destructor
RCSChiSquare::~RCSChiSquare()
{
  delete gen;
}

/// Get an RCS based on the Swerling II model and the mean RCS
rsFloat RCSChiSquare::SampleModel()
{
  return gen->GetSample();
}

//
// Target Implementation
//

//Default constructor for Target object
Target::Target(Platform *platform, const std::string &name):
  Object(platform, name)
{
  model = 0;
}

//Default destructor for Target object
Target::~Target()
{
  delete model;
}

/// Get the target polarization matrix
PSMatrix Target::GetPolarization() const
{
  return psm;
}

/// Set the target polarization matrix
void Target::SetPolarization(const PSMatrix &in)
{
  psm = in;
}

/// Set the target fluctuation model
void Target::SetFluctuationModel(RCSModel *in)
{
  model = in;
}

//
// IsoTarget Implementation
//

/// Constructor
IsoTarget::IsoTarget(Platform *platform, const std::string &name, rsFloat rcs):
  Target(platform, name),
  rcs(rcs)
{
}

/// Destructor
IsoTarget::~IsoTarget()
{
}

/// Return the RCS at the given angle
rsFloat IsoTarget::GetRCS(SVec3 &inAngle, SVec3 &outAngle) const
{
  if (model)
    return rcs*model->SampleModel();
  else
    return rcs;
}

//
// FileTarget Implementation
//

/// Constructor
FileTarget::FileTarget(Platform *platform, const std::string &name, const std::string &filename):
  Target(platform, name)
{
  //Create the objects for azimuth and elevation interpolation
  azi_samples = new InterpSet();
  elev_samples = new InterpSet();
  //Load the data from the description file into the interpolation objects
  LoadRCSDescription(filename);
}

/// Destructor
FileTarget::~FileTarget()
{
  delete azi_samples;
  delete elev_samples;
}


/// Return the RCS at the given angle
rsFloat FileTarget::GetRCS(SVec3 &inAngle, SVec3 &outAngle) const
{
  //Currently uses a half angle approximation, this needs to be improved
  SVec3 t_angle = inAngle+outAngle;
  rsFloat RCS = std::sqrt(azi_samples->Value(t_angle.azimuth/2.0)*elev_samples->Value(t_angle.elevation/2.0));
  if (model)
    return RCS*model->SampleModel();
  else
    return RCS;
}

namespace {

//Load samples of gain along an axis (not a member of FileAntenna)
void LoadTargetGainAxis(InterpSet *set, TiXmlHandle &axisXML)
{
  rsFloat angle;
  rsFloat gain;
  //Step through the XML file and load all the gain samples
  TiXmlHandle tmp = axisXML.ChildElement("rcssample", 0);
  for (int i = 0; tmp.Element() != 0; i++) {
    //Load the angle of the gain sample
    TiXmlHandle angleXML = tmp.ChildElement("angle", 0);
    if (!angleXML.Element())
      throw std::runtime_error("[ERROR] Misformed XML in target description: No angle in rcssample");
    angle = GetNodeFloat(angleXML);
    //Load the gain of the gain sample
    TiXmlHandle gainXML = tmp.ChildElement("rcs", 0);
    if (!gainXML.Element())
      throw std::runtime_error("[ERROR] Misformed XML in target description: No rcs in rcssample");
    gain = GetNodeFloat(gainXML);
    //Load the values into the interpolation table
    set->InsertSample(angle, gain);
    //Get the next gainsample in the file
    tmp = axisXML.ChildElement("rcssample", i);
  }
}

} //End of Anon. namespace

///Load data from the RCS description file
void FileTarget::LoadRCSDescription(const std::string& filename)
{
  TiXmlDocument doc(filename.c_str());
  //Check the document was loaded correctly
  if (!doc.LoadFile())
    throw std::runtime_error("[ERROR] Could not load target description from "+filename);
  //Get the XML root node
  TiXmlHandle root(doc.RootElement());
  //Load the gain samples along the elevation axis
  TiXmlHandle tmp = root.ChildElement("elevation", 0);
  if (!tmp.Element())
    throw std::runtime_error("[ERROR] Malformed XML in target description: No elevation pattern definition");
  LoadTargetGainAxis(elev_samples, tmp);
  //Load the gain samples along the azimuth axis
  tmp = root.ChildElement("azimuth", 0);
  if (!tmp.Element())
    throw std::runtime_error("[ERROR] Malformed XML in target description: No azimuth pattern definition");
  LoadTargetGainAxis(azi_samples, tmp);
}

//
// Functions for creating objects of various target types
//

/// Create an isometric radiator target
Target* rs::CreateIsoTarget(Platform *platform, const std::string &name, rsFloat rcs)
{
  return new IsoTarget(platform, name, rcs);
}

/// Create a target, loading the RCS pattern from a file
Target* rs::CreateFileTarget(Platform *platform, const std::string &name, const std::string &filename)
{
  return new FileTarget(platform, name, filename);
}
