//Timing source for simulation - all objects must be slaved to a timing source
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//16 October 2006

#ifndef __RS_TIMING_H
#define __RS_TIMING_H

#include <config.h>
#include <boost/utility.hpp>
#include <string>
#include <vector>

namespace rs {
  //Forward definitions
  class ClockModelGenerator; //rstiming.h

  ///Timing controls the timing of systems to which it is attached
  class Timing : boost::noncopyable {
  public:
    /// Constructor
    Timing(const std::string &name);
    /// Destructor
    virtual ~Timing();
    /// Get the real time of a particular pulse
    virtual rsFloat GetPulseTimeError() const = 0;
    /// Get the next sample of time error for a particular pulse
    virtual rsFloat NextNoiseSample() = 0;
    /// Skip a sample, computing only enough to preserve long term correlations
    virtual void SkipSamples(long long samples) = 0;
    /// Get the name of the timing source
    std::string GetName() const;
  private:
    std::string name; //!< The name of the prototype this is based on
  };

  /// Prototypes for timing sources used by system devices
  class PrototypeTiming {
  public:
    /// Constructor
    PrototypeTiming(const std::string &name);
    /// Add an alpha and a weight to the timing prototype
    void AddAlpha(rsFloat alpha, rsFloat weight);
    /// Get the alphas and weights from the prototype
    void GetAlphas(std::vector<rsFloat> &get_alphas, std::vector<rsFloat> &get_weights) const;
    /// Get the phase offset
    rsFloat GetPhaseOffset() const;
    /// Get the frequency offset
    rsFloat GetFreqOffset() const;
    /// Get the frequency
    rsFloat GetFrequency() const;
    /// Get the value of the sync on pulse flag
    bool GetSyncOnPulse() const;
    /// Set a constant frequency offset
    void AddFreqOffset(rsFloat offset);
    /// Set a constant phase offset
    void AddPhaseOffset(rsFloat offset);
    /// Set a random frequency offset
    void AddRandomFreqOffset(rsFloat stdev);
    /// Set a random phase offset
    void AddRandomPhaseOffset(rsFloat stdev);
    /// Set the base frequency of the clock model
    void SetFrequency(rsFloat freq);
    /// Get the name of the prototype
    std::string GetName() const;
    /// Set the sync on pulse flag -- timing error resets at the start of the pulse
    void SetSyncOnPulse();

  private:
    std::string name; //!< The name of the prototype timing source
    std::vector<rsFloat> alphas; //!< Alpha parameters for 1/f^alpha clock model
    std::vector<rsFloat> weights; //!< Weights for 1/f^alpha clock model
    rsFloat freq_offset; //!< Constant frequency offset 
    rsFloat phase_offset; //!< Constant phase offset
    rsFloat random_phase; //!< Standard deviation of random phase offset
    rsFloat random_freq; //!< Standard deviation of random frequency offset
    rsFloat frequency; //!< The nominal oscillator frequency
    bool synconpulse; //!< Reset timing error at the start of each pulse
  };

  /// Implementation of clock timing based on the 1/f model with linear interpolation
  class ClockModelTiming: public Timing {
  public:
    /// Constructor
    ClockModelTiming(const std::string &name);
    /// Destructor
    virtual ~ClockModelTiming();
    /// Get the next sample of time error for a particular pulse
    virtual rsFloat NextNoiseSample();
    /// Skip a sample, computing only enough to preserve long term correlations
    virtual void SkipSamples(long long samples);
    /// Reset the clock phase error to zero
    void Reset();
    /// Get the value of the sync on pulse flag
    bool GetSyncOnPulse() const;
    /// Initialize the clock model generator
    virtual void InitializeModel(const PrototypeTiming *timing);
    /// Get the real time of a particular pulse
    virtual rsFloat GetPulseTimeError() const;
    /// Get the carrier frequency of the modelled clock
    rsFloat GetFrequency() const;
    /// Return the enabled state of the clock model
    bool Enabled();
  private:
    bool enabled; //!< Is the clock model going to produce non-zero samples?
    ClockModelGenerator *model; //!< Clock model for intra-pulse samples
    std::vector<rsFloat> alphas; //!< Alpha parameters for 1/f^alpha clock model
    std::vector<rsFloat> weights; //!< Weights for 1/f^alpha clock model
    rsFloat frequency; //!< Carrier frequency of the modelled clock
    bool synconpulse; //!< Reset the timing at the start of each pulse
  };

}

#endif
