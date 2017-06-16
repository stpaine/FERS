//Digital Signal Processing support functions
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//30 July 2007

#include "rsdsp.h"
using namespace rs;

#include <stdexcept>
#include <cmath>
#include <cstring>
#include <string.h>
#include "rsdebug.h"
#include "rsparameters.h"

//
// Support Functions
//
namespace {
  //Calculate sin(pi * x) / (pi*x)
  rsFloat Sinc(rsFloat x) {
    if (x == 0)
      return 1.0;
    return std::sin(x*M_PI)/(x*M_PI);
  }

  /// Create a FIR filter using the Blackman window
  rsFloat *BlackmanFIR(rsFloat cutoff, int &length) {
    // Use double the render filter length, for faster rolloff than the render filter
    length = rsParameters::render_filter_length() * 2;
    rsFloat *coeffs = new rsFloat[length];
    rsFloat N = length / 2.0;
    for (int i = 0; i < length; i++) {
      rsFloat filt = Sinc(cutoff*(i - N));
      // We use the Blackman window, for a suitable tradeoff between rolloff and stopband attenuation
      // Equivalent Kaiser beta = 7.04 (Oppenhiem and Schaffer, Hamming)
      rsFloat window = 0.42 - 0.5*cos(M_PI*i/N) + 0.08*cos(2*M_PI*i/N);
      coeffs[i] = filt*window;
    }
    return coeffs;
  }
}

/// Upsample size samples stored *in by an integer ratio and store the result in (pre-allocated) out
// TODO: this would be better as a multirate upsampler
// In fact, the whole scheme is currently sub-optimal - we could use better filters, better windows and a better approach
// it works ok for now, users seeking higher accuracy can oversample outside FERS until this is fixed
void rs::Upsample(const rsComplex *in, int size, rsComplex *out, int ratio)
{
  /// Design the FIR filter for de-imaging
  int filt_length;
  rsFloat *coeffs = BlackmanFIR(1/rsFloat(ratio), filt_length);

  // Temporary buffer for zero padding and results
  rsComplex *tmp = new rsComplex[size*ratio+filt_length];
  for (int i = 0; i < size*ratio+filt_length; i++)
    tmp[i] = 0;
  /// Stuff the data with a suitable number of zeros
  for (int i = 0; i < size; i++) {
    tmp[i*ratio] = in[i];
    for (int j = 1; j < (ratio-1); j++)
      tmp[i*ratio+j] = 0;
  }
  // Create a FIR filter object
  FIRFilter filt(coeffs, filt_length);
  filt.Filter(tmp, size*ratio+filt_length);
  // Copy results to output buffer
  for (int i = 0; i < size*ratio; i++) {
    out[i] = tmp[i+filt_length/2-1];
  }
  // Clean up
  delete[] tmp;
  delete[] coeffs;
  
}

/// Upsample size samples stored *in by an integer ratio and store the result in (pre-allocated) out
// TODO: This would be better (and much faster) as a multirate downsampler
void rs::Downsample(const rsComplex *in, int size, rsComplex *out, int ratio)
{
  /// Design the FIR filter for anti-aliasing
  int filt_length;  
  rsFloat *coeffs = BlackmanFIR(1/rsFloat(ratio), filt_length);
  // Temporary buffer for zero padding and results
  rsComplex *tmp = new rsComplex[size+filt_length];
  for (int i = size-1; i < size+filt_length; i++)
    tmp[i] = 0;
  // Copy the input into the temporary buffer, leaving zero padding at the end
  for (int i = 0; i < size; i++)
    tmp[i] = in[i];
  // Run the anti aliasing filter on the data
  FIRFilter filt(coeffs, filt_length);
  filt.Filter(tmp, size+filt_length);

    
  //Copy the results to the output buffer
  for (int i = 0; i < size/ratio; i++) {
    out[i] = tmp[i*ratio+filt_length/2]/rsFloat(ratio);
    //    printf("%f+%fi\n", out[i].real(), out[i].imag());
  }
  // Clean up
  delete[] coeffs;
  delete[] tmp;
}

//
// Filter Implementation
//

/// Constructor
DSPFilter::DSPFilter()
{
}

/// Destructor
DSPFilter::~DSPFilter()
{
}

//
// IIRFilter Implementation
//

/// Constructor
IIRFilter::IIRFilter(const std::vector<rsFloat> &den_coeffs, const std::vector<rsFloat> &num_coeffs)
{
  //Get the filter order
  order = den_coeffs.size();
  //Check the filter order
  if (order != num_coeffs.size())
      throw std::logic_error("IIRFilter does not currently support mixed order filters");
  //Allocate memory to store co-efficients and state
  a = new rsFloat[order];
  b = new rsFloat[order];
  w = new rsFloat[order];
  //Load the co-efficients from the vectors into the arrays  
  for (unsigned int i = 0; i < order; i++) {
    a[i] = den_coeffs[i];
    b[i] = num_coeffs[i];
    w[i] = 0;
  }
}

/// Constructor
IIRFilter::IIRFilter(const rsFloat *den_coeffs, const rsFloat *num_coeffs, unsigned int order):
  order(order)
{
  a = new rsFloat[order];
  b = new rsFloat[order];
  w = new rsFloat[order];
  // Load the coefficients into the arrays
  for (unsigned int i = 0; i < order; i++) {
    a[i] = den_coeffs[i];
    b[i] = num_coeffs[i];
    w[i] = 0;
  }
}

/// Destructor
IIRFilter::~IIRFilter()
{
  //Clean up the co-efficients and state
  delete[] a;
  delete[] b;
  delete[] w;
}

/// Pass a single sample through the filter
rsFloat IIRFilter::Filter(rsFloat sample) 
{
  //Shift the filter state
  for (unsigned int j = order-1; j > 0; j--)
    w[j] = w[j-1];
  // Calculate w[0]
  w[0] = sample;
  for (unsigned int j = 1; j < order; j++)
    w[0] -= a[j]*w[j];
  //Calculate y[n]
  rsFloat out = 0;
  for (unsigned int j = 0; j < order; j++)
    out += b[j]*w[j];
  return out;
}

/// Pass an array of samples through the filter, filtering in place
void IIRFilter::Filter(rsFloat *samples, int size)
{
  for (int i = 0; i < size; i++)
    {
      //Shift the filter state
      for (unsigned int j = order-1; j > 0; j--)
	w[j] = w[j-1];
      // Calculate w[0]
      w[0] = samples[i];
      for (unsigned int j = 1; j < order; j++)
	w[0] -= a[j]*w[j];
      //Calculate y[n]
      samples[i] = 0;
      for (unsigned int j = 0; j < order; j++)
	samples[i] += b[j]*w[j];
    }
}

//
// FIRFilter implementation
//

/// Constructor
FIRFilter::FIRFilter(const std::vector<rsFloat> &coeffs)
{
  //Get the filter order
  order = coeffs.size();
  //Allocate memory to store co-efficients and state
  filter = new rsFloat[order];
  w = new rsFloat[order];
  //Load the co-efficients from the vectors into the arrays  
  for (unsigned int i = 0; i < order; i++) {
    filter[i] = coeffs[i];
    w[i] = 0;
  }
}

/// Constructor from coeffs
FIRFilter::FIRFilter(const rsFloat* coeffs, int count) {
  order = count;
  // Allocate memory to store co-efficients and state
  filter = new rsFloat[order];
  w = new rsFloat[order];
  // Load the co-efficients
  for (unsigned int i = 0; i < order; i++) {
    filter[i] = coeffs[i];
    w[i] = 0;
  }  
}

/// Destructor
FIRFilter::~FIRFilter()
{
  // Clean up memory
  delete[] filter;
  delete[] w;
}

/// Pass a single sample through the filter
inline rsFloat FIRFilter::Filter(rsFloat sample)
{


  return 0;
}

/// Pass an array of samples through the filter, filtering in place
  // See Oppenheim and Scaffer, section 6.5 "Basic Network Structures for FIR Systems"
// TODO: Implement one of the more efficient FIR filter forms
inline void FIRFilter::Filter(rsFloat *samples, int size)
{
  // Allocate memory for a delay line, equal to the filter length
  rsFloat* line = new rsFloat[order];
  std::memset(line, 0, sizeof(rsFloat)*order);
  // Perform the inplace convolution with the pulse
  for (int i = 0; i < size; i++) {
    line[0] = samples[i];
    rsFloat res = 0;
    for (unsigned int j = 0; j < order; j++)
      res += line[order-j-1]*filter[j];
    samples[i] = res;
    //Move the line over by one sample
    for (int j = order; j > 0; j--)
      line[j] = line[j-1];
  }
  //Clean up
  delete[] line;
}

/// Pass an array of complex samples through the filter, filtering in place
inline void FIRFilter::Filter(std::complex<rsFloat> *samples, int size) 
{
  // Allocate memory for a delay line, equal to the filter length
  rsComplex* line = new rsComplex[order];
  for (unsigned int i = 0; i < order; i++)
    line[i] = 0;
  // Perform the inplace convolution with the pulse
  for (int i = 0; i < size; i++) {
    line[0] = samples[i];
    rsComplex res = 0;
    for (unsigned int j = 0; j < order; j++)
      res += line[order-j-1]*filter[j];
    samples[i] = res;
    //Move the line over by one sample
    for (int j = order-1; j > 0; j--)
      line[j] = line[j-1];
  }
  //Clean up
  delete[] line;
}

//
// ARFilter implementation
//


/// Constructor
ARFilter::ARFilter(const std::vector<rsFloat> &coeffs)
{
  //Get the filter order
  order = coeffs.size();
  //Allocate memory to store co-efficients and state
  filter = new rsFloat[order];
  w = new rsFloat[order];
  //Load the co-efficients from the vectors into the arrays  
  for (unsigned int i = 0; i < order; i++) {
    filter[i] = coeffs[i];
    w[i] = 0;
  }
}

/// Destructor
ARFilter::~ARFilter()
{
  // Clean up memory
  delete[] filter;
  delete[] w;
}

/// Pass a single sample through the filter
rsFloat ARFilter::Filter(rsFloat sample) 
{
  //Shift the filter state
  for (unsigned int j = order-1; j > 0; j--)
    w[j] = w[j-1];
  // Calculate w[0]
  w[0] = sample;
  for (unsigned int j = 1; j < order; j++)
    w[0] -= filter[j]*w[j];
  //Return the output value of the filter
  return w[0];
}

/// Pass an array of samples through the filter, filtering in place
void ARFilter::Filter(rsFloat *samples, int size)
{
  for (int i = 0; i < size; i++)
    {
      //Shift the filter state
      for (unsigned int j = order-1; j > 0; j--)
	w[j] = w[j-1];
      // Calculate w[0]
      w[0] = samples[i];
      for (unsigned int j = 1; j < order; j++)
	w[0] -= filter[j]*w[j];
      //Calculate y[n]
      samples[i] = w[0];
    }
}

//
// Upsampler implementation
//

/// Constructor
Upsampler::Upsampler(int ratio):
  ratio(ratio)
{
  //Create the FIR interpolation filter
  filter_size = 8*ratio+1; // 8*ratio should give adequate performance  
  //Allocate memory for the filter bank
  filterbank = new rsFloat[filter_size];
  // Simple windowed sinc filter design procedure
  for (int i = 0; i < filter_size; i++) {
    // The Hamming window provides a solid tradeoff between rolloff and stopband attenuation
    rsFloat window_value = 0.54 - 0.46 * std::cos(2*M_PI*i/(rsFloat)(filter_size));
    rsFloat filter_value = Sinc(1.0/(rsFloat)(ratio)*(i-filter_size/2));
    filterbank[i] = filter_value * window_value;
  }
  //Allocate memory for the sample state
  sample_memory = new rsFloat[filter_size/ratio+1];
  //Clear sample memory to zero
  for (int i = 0; i < filter_size/ratio+1; i++)
    sample_memory[i] = 0;
}

/// Destructor
Upsampler::~Upsampler()
{
  // Clean up filter and state
  delete[] filterbank;
  delete[] sample_memory;
}

//Get a sample, from either the provided pointer or sample memory
inline rsFloat Upsampler::GetSample(rsFloat *samples, int n)
{
 if (n >= 0)
    return samples[n];
  else
    return sample_memory[n+filter_size];
}

/// Upsamples a signal and applies an anti-imaging filter
void Upsampler::Upsample(rsFloat *insamples, int in_size, rsFloat *outsamples, int out_size)
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
  //Update the sample history
  int transfer_size = filter_size/ratio+1;
  if (in_size >= transfer_size)
    memcpy(sample_memory, &(insamples[in_size-transfer_size]), transfer_size*sizeof(rsFloat));
  else {
    // Shift existing samples
    for (int i = 0; i < (transfer_size-in_size); i++)
      sample_memory[i] = sample_memory[i+in_size];
    // Add new samples to the end of the buffer
    for (int i = 0; i < in_size; i++)
      sample_memory[i+transfer_size-in_size] = insamples[i];
  }
}

//
// DecadeUpsample Implementation
//

/// Constructor
DecadeUpsampler::DecadeUpsampler() {

  /// 11th order elliptic lowpass at 0.1fs
  rsFloat den_coeffs[12] = {1.0,
			    -10.301102119865,
			    48.5214567642597,
			    -137.934509572412,
			    262.914952985445,
			    -352.788381841481,
			    340.027874008585,
			    -235.39260470286,
			    114.698499845697,
			    -37.4634653062448,
			    7.38208765922137,
			    -0.664807695826097};

  rsFloat num_coeffs[12] = {   2.7301694322809e-06,
			       -1.8508123430239e-05,
			       5.75739466753894e-05,
			       -0.000104348734423658,
			       0.000111949190289715,
			       -4.9384188225528e-05,
			       -4.9384188225522e-05,
			       0.00011194919028971,
			       -0.000104348734423656,
			       5.75739466753884e-05,
			       -1.85081234302388e-05,
			       2.73016943228086e-06  };
  //Initialize the anti-imaging filter
  filter = new IIRFilter(den_coeffs, num_coeffs, 12);
}

/// Destructor
DecadeUpsampler::~DecadeUpsampler() {
  delete filter;
}
  

///Upsample a single sample at a time
void DecadeUpsampler::Upsample(rsFloat sample, rsFloat *out)
{
  // Prepare the output array
  out[0] = sample;
  for (int i = 1; i < 10; i++)
    out[i] = 0;
  // Filter in place
  filter->Filter(out, 10);
}

// Upsample a whole batch of samples
void DecadeUpsampler::Upsample(rsFloat *in, int count, rsFloat *out)
{
  /// Prepare the array for filtering
  for (int i = 0; i < count; i++) {
    out[i*10] = in[i];
    for (int j = 1; j < 10; j++)
      out[i*10+j] = 0;
  }
  /// Filter in place
  filter->Filter(out, count*10);
}

