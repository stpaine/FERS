//rsnoise.h Functions and classes to generate noise of various types
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//14 August 2006

#ifndef __RSNOISE_H
#define __RSNOISE_H

#include <config.h>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/utility.hpp>
#include <vector>
#include "rsdsp.h"
#include "rspython.h"


namespace rsNoise {
  /// Initialize the random number generator code (must be called once, after the loading of the script)
  void InitializeNoise();
  /// Clean up the noise code
  void CleanUpNoise();
  /// Return a single sample of white gaussian noise
  rsFloat WGNSample(rsFloat stddev);
  /// Return a single uniformly distributed sample in [0, 1]
  rsFloat UniformSample();
  /// Calculate noise power from the temperature
  rsFloat NoiseTemperatureToPower(rsFloat temperature, rsFloat bandwidth);
}

namespace rs {

/// NoiseGenerator - base class for different types of noise generator
class NoiseGenerator: boost::noncopyable {
public:
  /// Constructor
  NoiseGenerator();
  /// Destructor
  virtual ~NoiseGenerator();
  /// Get a single random sample
  virtual rsFloat GetSample() = 0;
};

/// Generator of Gaussian white noise
class WGNGenerator: public NoiseGenerator {
public:
  /// Constructor
  WGNGenerator(rsFloat stddev);
  /// Default Constructor
  WGNGenerator();
  /// Destructor
  virtual ~WGNGenerator();
  /// Get a single random sample
  virtual rsFloat GetSample();
private:
  boost::normal_distribution<rsFloat> dist; //!< PRNG distribution
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > *gen; //!< Variate Generator (see boost::random docs)
  rsFloat temperature; //!< Noise temperature
};

  /// Gamma distributed noise generator
  class GammaGenerator: public NoiseGenerator {
  public:
    /// Constructor
    GammaGenerator(rsFloat k); //x_bar is the 'scale' parameter and k is the 'shape' parameter
    /// Destructor
    virtual ~GammaGenerator();
    /// Get a single random sample
    virtual rsFloat GetSample();
    /// Operator to get a random sample
    rsFloat operator()();
  private:
    boost::gamma_distribution<rsFloat> dist; //!< Gamma distribution function
    boost::variate_generator<boost::mt19937&, boost::gamma_distribution<> > gen; //!< Variate generator
  };


  /// Single branch of the multirate clock generator model
  class FAlphaBranch: boost::noncopyable {
  public:
    /// Constructor
    FAlphaBranch(rsFloat ffrac, unsigned int fint, FAlphaBranch *pre, bool last);
    /// Destructor
    ~FAlphaBranch();
    /// Get a sample from the branch
    rsFloat GetSample();
    /// Flush the buffer, and refull with samples
    void Flush(rsFloat scale);
  private:
    /// Initialize the filters, etc.
    void Init();
    /// Clean up the filters
    void Clean();
    /// Refill the buffer with samples
    void Refill();
    /// Calculate a single sample
    rsFloat CalcSample();

    IIRFilter *shape_filter; //!< The filter for shaping the noise
    rsFloat shape_gain; //!< Gain of shaping filter
    IIRFilter *integ_filter; //!< Integrator filter
    rsFloat integ_gain;  //!< Gain of integration filter

    rsFloat upsample_scale; //!< Scaling factor for upsampling
    IIRFilter *highpass; //!< Highpass filter
    FAlphaBranch *pre; //!< Next lower branch in the structure
    bool last; //!< If this filter is the top branch, don't upsample
    DecadeUpsampler *upsampler; //!< Upsampler for this branch
    rsFloat *buffer; //!< Buffer for storing samples from the upsampler
    unsigned int buffer_samples; //!< Number of samples available in the buffer
    rsFloat ffrac; //!< Fractional part of filter curve
    rsFloat fint; //!< Integer part of filter curve
    rsFloat offset_sample; //!< Sample from the branch below us
    bool got_offset; //!< Are we waiting for the offset
    rsFloat pre_scale; //!< Previous branch scale factor
    friend class MultirateGenerator;
  };


  /// Class to generate 1/f^alpha noise based on multirate approach
  class MultirateGenerator: public NoiseGenerator {
  public:
    /// Constructor
    MultirateGenerator(rsFloat alpha, unsigned int branches);
    /// Destructor
    ~MultirateGenerator();
    /// Get a single noise sample
    rsFloat GetSample();
    /// Skip a number of samples, preserving correlations of period longer than the sample count
    void SkipSamples(long long samples);
    /// Reset the output to zero
    void Reset();
  private:
    rsFloat scale; //!< Scale for normalizing values
    /// Create the branches of the filter structure tree
    void CreateTree(rsFloat falpha, int fint, unsigned int branches);
    /// Get the co-efficients of the shaping filter 
    FAlphaBranch *topbranch; //!< Top branch of the filter structure tree
  };

  /// Class to generate noise based on the weighted sum of 1/f^alpha noise
  class ClockModelGenerator: public NoiseGenerator {
  public:
    /// Constructor
    ClockModelGenerator(const std::vector<rsFloat> &alpha, const std::vector<rsFloat> &in_weights, rsFloat frequency, rsFloat phase_offset, rsFloat freq_offset, int branches);
    /// Destructor
    ~ClockModelGenerator();
    /// Get a single noise sample
    rsFloat GetSample();
    /// Skip noise samples, calculating only enough to preserve long-term correlations
    void SkipSamples(long long samples);
    /// Reset the noise to zero
    void Reset();
    /// Is the generator going to produce non-zero samples
    bool Enabled();
  private:
    std::vector<MultirateGenerator *> generators; // The multirate generators which generate noise in each band
    std::vector<rsFloat> weights; // Weight of the noise from each generator
    rsFloat phase_offset; //!< Offset from nominal phase
    rsFloat freq_offset; //!< Offset from nominal base frequency
    rsFloat frequency; //!< Nominal base frequency
    unsigned long count; //!< Number of samples generated by this generator

  };

  /// Generator of noise using a python module
  class PythonNoiseGenerator: public NoiseGenerator {
  public:
    ///Constructor
    PythonNoiseGenerator(const std::string& module, const std::string& function);
    ///Destructor
    ~PythonNoiseGenerator();
    ///Get a single noise sample
    rsFloat GetSample();
  private:
    rsPython::PythonNoise generator;
  };


}

#endif
