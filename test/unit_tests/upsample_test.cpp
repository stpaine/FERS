// Test the polyphase upsampler code from rsdsp.cpp
#include <iostream>
#include <complex>
#include "cycle.h"


using namespace std;

typedef double rsFloat;

  /// Upsamples a signal and applies an anti-imaging filter
  // Implemented using polyphase FIR filter with windowed sinc
  class Upsampler {
  public:
    /// Constructor (ratio of upsampling, co-efficients of anti-imaging filter)
    Upsampler(int ratio);
    /// Destructor
    ~Upsampler();
    /// Upsample the given sample to a pre-allocated target
    void Upsample(std::complex<rsFloat> *insamples, int in_size, std::complex<rsFloat> *outsamples, int out_size);
  private:
    int ratio; //<! Upsampling ratio
    rsFloat *filterbank; //<! FIR polyphase filter bank
    std::complex<rsFloat> *sample_memory; //<! Last samples used, to allow seamless upsampling in blocks
    int filter_size; //!< Length of the interpolation filter
    //Get a sample, from either the provided pointer or sample memory
    inline std::complex<rsFloat> GetSample(std::complex<rsFloat> *samples, int n);
  };

//
// Upsampler implementation
//

/// Constructor
Upsampler::Upsampler(int ratio):
  ratio(ratio)
{
  //Create the FIR interpolation filter
  int filter_size = 8*ratio+1; // 8*ratio should give adequate performance
  //Allocate memory for the filter bank
  filterbank = new rsFloat[filter_size];
  // Simple windowed sinc filter design procedure
  for (int i = 0; i < filter_size; i++) {
    // The Hamming window provides a solid tradeoff between rolloff and stopband attenuation
    rsFloat window_value = 0.54 - 0.46 * std::cos(2*M_PI*i/double(filter_size));
    rsFloat filter_value = Sinc(1.0/double(ratio)*(i-filter_size/2));
    filterbank[i] = filter_value * window_value;
  }
  //Allocate memory for the sample state
  sample_memory = new std::complex<rsFloat>[filter_size];
  //Clear sample memory to zero
  for (int i = 0; i < filter_size; i++)
    sample_memory = 0;
}

/// Destructor
Upsampler::~Upsampler()
{
  // Clean up filter and state
  delete[] filterbank;
  delete[] sample_memory;
}

//Get a sample, from either the provided pointer or sample memory
inline std::complex<rsFloat> Upsampler::GetSample(std::complex<rsFloat> *samples, int n)
{
  if (n >= 0)
    return samples[n];
  else
    return sample_memory[n+filter_size];
}

/// Upsamples a signal and applies an anti-imaging filter
void Upsampler::Upsample(std::complex<rsFloat> *insamples, int in_size, std::complex<rsFloat> *outsamples, int out_size)
{
  //Check the target array size
  if (out_size != (ratio*in_size))
    throw std::runtime_error("Target array size is not correct in Upsample");
  // Polyphase upsampler implementation
  // See fers_upsample_p.m in the documentation for more details
  // Follows the diagram in section 4.7.4 "Polyphase Implementation of Interpolation Filters" of
  // Discrete Time Signal Processing, 2nd ed., Oppenheim and Schafer
  for (int i = 0, branch = 0; i < in_size; i++, branch++)
    {
      if (branch >= ratio)
	branch = 0;
      outsamples[i] = 0;
      // Apply the branch filter to the data
      for (int j = branch; j < filter_size; j += ratio)
	outsamples[i] += filterbank[j] * GetSample(insamples, i-j/ratio);
    }
}

void SimpleTest()
{
  //Create the test vector
  Complex *cmp = 
