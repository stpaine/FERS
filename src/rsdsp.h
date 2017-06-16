//Digital Signal Processing support functions
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//30 July 2007

#ifndef __RS_DSP_H
#define __RS_DSP_H

#include <config.h>
#include <boost/utility.hpp>
#include <vector>
#include <complex>
#include <map>
#include "rsradarwaveform.h"

namespace rs {

  //
  // Support functions
  //

  /// Upsample size samples stored *in by an integer ratio and store the result in (pre-allocated) out
  void Upsample(const rsComplex *in, int size, rsComplex *out, int ratio);
  /// Downsample size samples stored *in by an integer ratio and store the result in (pre-allocated) out
  void Downsample(const rsComplex *in, int size, rsComplex *out, int ratio);

  /// Filter, parent class for digital filters
  class DSPFilter: boost::noncopyable {
  public:
    /// Constructor
    DSPFilter();
    /// Destructor
    virtual ~DSPFilter();
    /// Pass a single sample through the filter
    virtual rsFloat Filter(rsFloat sample) = 0;
    /// Pass an array of samples through the filter, filtering in place
    virtual void Filter(rsFloat *samples, int size) = 0;
  };

  /// IIR (ARMA) Digital Filter, implemented with Direct Form II
  // Supports filters of the type A(z)/B(z)
  class IIRFilter: public DSPFilter {
  public:
    /// Constructor
    IIRFilter(const std::vector<rsFloat> &den_coeffs, const std::vector<rsFloat> &num_coeffs);
    /// Constuctor
    IIRFilter(const rsFloat *den_coeffs, const rsFloat *num_coeffs, unsigned int order);
    /// Destructor
    virtual ~IIRFilter();
    /// Pass a single sample through the filter
    virtual rsFloat Filter(rsFloat sample);
    /// Pass an array of samples through the filter, filtering in place
    virtual void Filter(rsFloat *samples, int size);
  private:
    rsFloat *w; //!< Past x values
    //    rsFloat *wy; //!< Past y values
    rsFloat *a; //!< Denominator co-efficients
    rsFloat *b; //!< Numerator co-efficients
    unsigned int order; //!< Filter order
  };

  /// FIR (MA) Digital Filter
  //Supports filters of the type B(z)/1
  class FIRFilter: public DSPFilter {
  public:
    /// Constructor
    FIRFilter(const std::vector<rsFloat> &coeffs);
    FIRFilter(const rsFloat* coeffs, int count);
    /// Destructor
    virtual ~FIRFilter();
    /// Pass a single sample through the filter
    virtual rsFloat Filter(rsFloat sample);
    /// Pass an array of samples through the filter, filtering in place
    virtual void Filter(rsFloat *samples, int size);
    /// Pass an array of complex samples through the filter, filtering in place
    void Filter(std::complex<rsFloat> *samples, int size);
  private:
    rsFloat *w; //!< Filter state
    rsFloat *filter; //!< Filter coefficients
    unsigned int order; //!< Filter order
  };

  /// Auto Regressive (AR) Digital Filter
  //Supports filters of the type 1/A(z)
  class ARFilter: public DSPFilter {
  public:
    /// Constructor
    ARFilter(const std::vector<rsFloat> &coeffs);
    /// Destructor
    ~ARFilter();
    /// Pass a single sample through the filter
    virtual rsFloat Filter(rsFloat sample);
    /// Pass an array of samples through the filter, filtering in place
    virtual void Filter(rsFloat *samples, int size);
  private:
    rsFloat *w; //!< Filter state
    rsFloat *filter; //!< Filter coefficients
    unsigned int order; //!< Filter order
  };

  /// Upsamples a signal and applies an anti-imaging filter
  // Implemented using polyphase FIR filter with windowed sinc
  class Upsampler: boost::noncopyable {
  public:
    /// Constructor (ratio of upsampling, co-efficients of anti-imaging filter)
    Upsampler(int ratio);
    /// Destructor
    ~Upsampler();
    /// Upsample the given sample to a pre-allocated target
    void Upsample(rsFloat *insamples, int in_size, rsFloat *outsamples, int out_size);
  private:
    int ratio; //!< Upsampling ratio
    rsFloat *filterbank; //!< FIR polyphase filter bank
    rsFloat *sample_memory; //!< Last samples used, to allow seamless upsampling in blocks
    int filter_size; //!< Length of the interpolation filter
    //Get a sample, from either the provided pointer or sample memory
    inline rsFloat GetSample(rsFloat *samples, int n);
  };

  /// Upsample a signal by a factor of 10
  class DecadeUpsampler: boost::noncopyable {
  public:
    /// Constructor
    DecadeUpsampler();
    /// Destrubtor
    ~DecadeUpsampler();
    /// Upsample one sample at a time, out is array of ten samples
    void Upsample(rsFloat sample, rsFloat *out);
    /// Upsample a large block, out must be ten times bigger than in
    void Upsample(rsFloat *in, int count, rsFloat *out);
  private:
    /// Anti-imaging filter
    IIRFilter *filter;
  };
};

#endif
