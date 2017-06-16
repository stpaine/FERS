//signal.cpp
//Interface for the signal class
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//Started 24 May 2006

#ifndef __SIGNAL_H
#define __SIGNAL_H

#include <config.h>
#include <complex>
#include <string>
#include <boost/utility.hpp>
#include <boost/shared_array.hpp>
#include "rsradarwaveform.h"

///Forward declarations
namespace rs {
class Signal; //rssignal.h
}

namespace rsSignal {
  /// Type for storing signal
  typedef std::complex<rsFloat> complex;
  
  /// Add noise to the signal with the given temperature
  void AddNoise(rsFloat *data, rsFloat temperature, unsigned int size, rsFloat fs);
  /// Demodulate a frequency domain signal into time domain I and Q
  complex* IQDemodulate(rsFloat *data, unsigned int size, rsFloat phase);
  /// Simulate the effect of and ADC converter on the signal
  void ADCSimulate(complex *data, unsigned int size, int bits, rsFloat fullscale);
}

namespace rs {



  /// Class to store and process a time domain signal
class Signal: boost::noncopyable {
public:
  Signal(); //!< Default constructor
  ~Signal(); //!< Default destructor

  /// Clear deletes the data currently associated with the signal
  void Clear();
  /// Load data into the signal (time domain, complex)
  void Load(const rsSignal::complex *indata, unsigned int samples, rsFloat samplerate);
  /// Load data into the signal (time domain, real)
  void Load(const rsFloat *indata, unsigned int samples, rsFloat samplerate);
  
  /// Return the sample rate of the signal
  rsFloat Rate() const;
  /// Return the size, in samples of the signal
  unsigned int Size() const;
  /// Get a copy of the signal domain data
  rsFloat* CopyData() const;
  /// Get the number of samples of padding at the beginning of the pulse
  unsigned int Pad() const;
  /// Render the pulse with the specified doppler, delay and amplitude
  boost::shared_array<rs::rsComplex> Render(const std::vector<InterpPoint> &points, rsFloat trans_power, unsigned int &size, rsFloat frac_win_delay) const;
 private:
  /// The signal data
  rsSignal::complex* data;
  /// Size of the signal in samples
  unsigned int size;  
  /// The sample rate of the signal in the time domain
  rsFloat rate;
  /// Number of samples of padding at the beginning of the pulse
  unsigned int pad;
};

}
#endif
