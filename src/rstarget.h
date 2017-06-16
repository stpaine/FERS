//rstarget.h
//Defines the class for targets
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//Started: 26 April 2006
#ifndef __RSTARGET_H
#define __RSTARGET_H

#include <config.h>
#include "rsobject.h"
#include "rspath.h"
#include "rspolarize.h"

namespace rs {

  //Forward declarations 
  class InterpSet; //from rsinterp.h
  class GammaGenerator; //from rsnoise.h


  /// General RCS statistical model class
  class RCSModel {
  public:
    /// Destructor
    virtual ~RCSModel();
    /// Get an RCS based on the statistical model and the mean RCS
    virtual rsFloat SampleModel()  = 0;
  };

  /// RCS Statistical Model class supporting Swerling V
  class RCSConst: public RCSModel {
  public:
    /// Destructor
    virtual ~RCSConst();
    /// Return a constant RCS
    virtual rsFloat SampleModel();
  };

  /// RCS statistical model following Swerling's Chi-square (actually Gamma) model
  // See Swerling, "Radar Probability of Detection for Some Additional Target Cases", IEEE Trans. Aer. Elec. Sys., Vol 33, 1997
  class RCSChiSquare: public RCSModel {
  public:
    /// Constructor
    RCSChiSquare(rsFloat k); //k is the shape parameter for the distribution
    /// Destructor
    virtual ~RCSChiSquare();
    /// Get an RCS based on the Swerling II model and the mean RCS
    virtual rsFloat SampleModel();
  private:
    GammaGenerator *gen;
  };

  /// Target models a simple point target with a specified RCS pattern
  class Target: public Object {
  public:
    /// Constructor
    Target(Platform *platform, const std::string &name);
    /// Destructor
    virtual ~Target();
    /// Returns the Radar Cross Section at a particular angle
    virtual rsFloat GetRCS(SVec3 &inAngle, SVec3 &outAngle) const = 0;
    /// Get the target polarization matrix
    virtual PSMatrix GetPolarization() const;
    /// Set the target polarization matrix
    virtual void SetPolarization(const PSMatrix &in);
    /// Set the target fluctuation model
    virtual void SetFluctuationModel(RCSModel *in);
  protected:
    PSMatrix psm; //!< Polarization scattering matrix for target interaction
    RCSModel *model; //!< Statistical model of target RCS fluctuations
  };

  /// Target with an isotropic (constant with angle) RCS
  class IsoTarget: public Target {
  public:
    /// Constructor
    IsoTarget(Platform *platform, const std::string &name, rsFloat rcs);
    /// Destructor
    virtual ~IsoTarget();
    /// Return the RCS at the given angle
    virtual rsFloat GetRCS(SVec3 &inAngle, SVec3 &outAngle) const;
  private:
    rsFloat rcs; //!< Constant RCS
  };

  /// Target with an RCS interpolated from a table of values
  class FileTarget: public Target {
  public:
    /// Constructor
    FileTarget(Platform *platform, const std::string &name, const std::string &filename);
    /// Destructor
    virtual ~FileTarget();
    /// Return the RCS at the given angle
    virtual rsFloat GetRCS(SVec3 &inAngle, SVec3 &outAngle) const;
  private:
    rs::InterpSet* azi_samples; //!< Samples of RCS in the azimuth plane
    rs::InterpSet* elev_samples; //!< Samples of RCS in the elevation plane
    ///Load data from the RCS description file
    void LoadRCSDescription(const std::string& filename);
  };

  // Functions for creating objects of various target types

  /// Create an isometric radiator target
  Target* CreateIsoTarget(Platform *platform, const std::string &name, rsFloat rcs);
  
  /// Create a target, loading the RCS pattern from a file
  Target* CreateFileTarget(Platform *platform, const std::string &name, const std::string &filename);


}

#endif
