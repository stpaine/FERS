//rsradarwaveform.cpp
//Classes for different types of radar waveforms
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//Started: 24 May 2006

#include <cmath>
#include <stdexcept>
#include <fstream>
#include <boost/scoped_array.hpp>
#include <ctype.h>
#include "rsradarwaveform.h"
#include "rsparameters.h"
#include "rsdebug.h"
#include "rssignal.h"
#include "rshdf5.h"

using namespace rs;

//
//RadarWaveform implementation
//

//Default constructor
RadarSignal::RadarSignal(std::string name, rsFloat power, rsFloat carrierfreq, rsFloat length, Signal* signal):
  name(name),
  power(power), 
  carrierfreq(carrierfreq),   
  length(length),
  signal(signal),
  polar(1.0,0) //Default to horiz polarization
{
  if (!signal)
    throw std::logic_error("RadarSignal cannot be constructed with NULL signal");
}

//Destructor
RadarSignal::~RadarSignal()
{
  delete signal;
}

//Return the power of the signal
rsFloat RadarSignal::GetPower() const
{
  return power;
}

//Get the carrier frequency
rsFloat RadarSignal::GetCarrier() const
{
  return carrierfreq;
}

//Get the name of the pulse
std::string RadarSignal::GetName() const
{
  return name;
}

//Get the native sample rate of the pulse
rsFloat RadarSignal::GetRate() const
{
  return signal->Rate();
}

/// Return the length of the pulse
rsFloat RadarSignal::GetLength() const
{
  return length;
}

/// Render the waveform to the target buffer
boost::shared_array<rsComplex> RadarSignal::Render(const std::vector<InterpPoint> &points, unsigned int &size, rsFloat frac_win_delay) const
{
  //Render the return pulse
  boost::shared_array<rsComplex> data = signal->Render(points, power, size, frac_win_delay);
  //Scale the return pulse by the signal power
  rsFloat scale = std::sqrt(power);
  for (unsigned int i = 0; i < size; i++)
    data[i] *= scale;
  return data;
}

/// Get the signal polarization
JonesVector RadarSignal::GetPolarization()
{
  return polar;
}

/// Set the signal polarization
void RadarSignal::SetPolarization(const JonesVector &in)
{
  polar = in;
}

//
// Functions to load signal data from files
//

/// Load the pulse from HDF5 file
RadarSignal* LoadPulseFromHDF5File(const std::string& name, const std::string &filename, rsFloat power, rsFloat carrierfreq)
{
  rsFloat rate;
  unsigned int size;
  rsComplex *data;
  // Load the data from the hdf5 file
  rshdf5::ReadPulseData(filename, &data, size, rate);
  //Create the signal object
  Signal *signal = new Signal();
  // Load the pulse into the signal object
  signal->Load(data, size, rate);
  delete[] data;
  // Create the RadarSignal
  rs::RadarSignal *any = new rs::RadarSignal(name, power, carrierfreq, size/rate, signal);
  return any;
}

/// Load the pulse from a CSV file
RadarSignal* LoadPulseFromCSVFile(const std::string& name, const std::string& filename, rsFloat power, rsFloat carrierfreq)
{
  ///Open the file
  std::ifstream ifile(filename.c_str());
  if (!ifile)
    throw std::runtime_error("Could not open "+filename+" to read pulse waveform");
  /// Read the length and sample rate
  rsFloat rlength, rate;
  ifile >> rlength; //length in samples
  ifile >> rate; //rate
  unsigned int length = static_cast<int>(rlength);
  //Allocate memory for the file contents
  boost::scoped_array<rsComplex> data(new rsComplex[length]);
  //Loop through reading the samples in the file
  unsigned int done = 0;
  while (!ifile.eof() && (done < length))
    ifile >> data[done++];
  if (done != length)
    throw std::runtime_error("Could not read pulse waveform from file "+filename);
  //Create the signal object with the data from the file
  Signal *signal = new Signal();
  signal->Load(data.get(), length, rate);
  //Create the pulse
  rs::RadarSignal *any = new rs::RadarSignal(name, power, carrierfreq, rlength/rate, signal);
  return any;
}

/// Load a pulse from a file and generate an anypulse
rs::RadarSignal* rsPulseFactory::LoadPulseFromFile(const std::string& name, const std::string& filename, rsFloat power, rsFloat carrierfreq)
{
  int ln = filename.length()-1;
  //Identify file types

  //Check for CSV extension
  if ((tolower(filename[ln]) == 'v') && (tolower(filename[ln-1]) == 's') && (tolower(filename[ln-2]) == 'c'))
    return LoadPulseFromCSVFile(name, filename, power, carrierfreq);  
  //Check for H5 extension
  else if ((tolower(filename[ln]) == '5') && (tolower(filename[ln-1]) == 'h'))
    return LoadPulseFromHDF5File(name, filename, power, carrierfreq);    
  // If neither of the above, complain
  else
    throw std::runtime_error("[ERROR] Unrecognised extension while trying to load "+filename);
  return 0;

}

