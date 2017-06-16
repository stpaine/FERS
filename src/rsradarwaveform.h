//rsradarwaveform.h
//Classes for different types of radar waveforms
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//8 June 2006

#ifndef __RS_RADARWAVEFORM_H
#define __RS_RADARWAVEFORM_H

#include <config.h>
#include <string>
#include <vector>
#include <complex>
#include <boost/shared_array.hpp>
#include <boost/utility.hpp>
#include "rspolarize.h"

namespace rs {

  ///Forward Declarations
  class Signal; //rssignal.h
  class JonesVector; //rspolarize.h

  ///Complex type for rendering operations
  typedef std::complex<rsFloat> rsComplex;

/// A continuous wave response interpolation point
struct InterpPoint {
  /// Constructor
  InterpPoint(rsFloat power, rsFloat start, rsFloat delay, rsFloat doppler, rsFloat phase, rsFloat noise_temperature);
  rsFloat power; //!< Peak power of the pulse (into 1ohm)
  rsFloat time; //!< Start time (seconds)
  rsFloat delay; //!< Pulse round trip time (seconds)
  rsFloat doppler; //!< Doppler shift (radians)
  rsFloat phase; //!< Phase (radians)
  rsFloat noise_temperature; //!< Noise temperature (kelvin)
};

/// The RadarWaveform class stores information about a pulse shape (or continuous wave wave)
 class RadarSignal: public boost::noncopyable
 {
 public:
   /// Default constructor
   RadarSignal(std::string name, rsFloat power, rsFloat carrierfreq, rsFloat length, Signal* signal);
   /// Destructor
   ~RadarSignal();
   /// Get the signal power
   rsFloat GetPower() const;
   /// Get the signal carrier frequency (Hz)
   rsFloat GetCarrier() const;
   /// Get the name of the signal
   std::string GetName() const;
   /// Return the native sample rate of the waveform
   rsFloat GetRate() const;
   /// Return the length of the pulse
   rsFloat GetLength() const;
   /// Render the pulse with the given parameters
   boost::shared_array<rsComplex> Render(const std::vector<InterpPoint> &points, unsigned int &size, rsFloat frac_win_delay) const;
   /// Get the signal polarization
   JonesVector GetPolarization();
   /// Set the signal polarization
   void SetPolarization(const JonesVector &in);
 private:   
   std::string name; //!< The name of the pulse
   rsFloat power; //!< Power of the signal transmitted (Watts in 1ohm)
   rsFloat carrierfreq; //!< Carrier frequency (Hz)
   rsFloat length; //!< Length of the signal (seconds)
   Signal* signal; //!< Transmitted Signal
   JonesVector polar; //!< Signal Polarization
   /// Default constructor
   RadarSignal();
 };

}

//Functions for creating radar waveforms
namespace rsPulseFactory {
 
/// Load a pulse from a file
rs::RadarSignal* LoadPulseFromFile(const std::string& name, const std::string& filename, rsFloat power, rsFloat carrierfreq);


}


#endif
