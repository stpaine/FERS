//rsradar.cpp
//Implementation of classes defined in rsradar.h
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//Started 21 April 2006

#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <limits>
#include "rsradar.h"
#include "rsdebug.h"
#include "rspulserender.h"
#include "rsresponse.h"
#include "rsantenna.h"
#include "rsparameters.h"
#include "rspath.h"
#include "rstiming.h"
#include "rsmultipath.h"
#include "rsplatform.h"

using namespace rs; //Import the rs namespace for clarity

//
// Radar Implementation
//

/// Default Constructor
Radar::Radar(const Platform *platform, const std::string& name):
  Object(platform, name),
  timing(0),
  antenna(0),
  attached(0),
  multipath_dual(false),
  multipath_reflect(0)
{
}

/// Default Destructor
Radar::~Radar()
{
}

/// Attach a receiver to the transmitter for a monostatic configuration
void Radar::MakeMonostatic(Radar* recv)
{
  if (attached)
    throw std::runtime_error("[BUG] Attempted to attach second receiver to transmitter");
  attached = recv;
}

/// Get the attached receiver
//Attached is likely to be 0 (NULL) - which means the transmitter does not share it's antenna
const Radar* Radar::GetAttached() const
{
  return attached;
}

/// Return whether the radar is monostatic
bool Radar::IsMonostatic() const
{
  return attached;
}

/// Set the transmitter's antenna
void Radar::SetAntenna(Antenna* ant)
{
  if (!ant)
    throw std::logic_error("[BUG] Transmitter's antenna set to null");
  antenna = ant;  
}

/// Return the antenna gain in the specified direction
rsFloat Radar::GetGain(const SVec3 &angle, const SVec3 &refangle, rsFloat wavelength) const
{
  return antenna->GetGain(angle, refangle, wavelength);
}

/// Get the noise temperature (including antenna noise temperature)
rsFloat Radar::GetNoiseTemperature(const SVec3 &angle) const
{
  return antenna->GetNoiseTemperature(angle);
}

/// Attach a timing object to the radar
void Radar::SetTiming(Timing* tim) {
  if (!tim)
    throw std::runtime_error("[BUG] Radar timing source must not be set to NULL");
  timing = tim;
}

/// Get a pointer to the radar's timing object
Timing* Radar::GetTiming() const
{
  if (!timing)
    throw std::runtime_error("[BUG] Radar::GetTiming called before timing set");
  return timing;
}

/// Is this object a virtual multipath dual?
bool Radar::IsMultipathDual() const
{
  return multipath_dual;
}

/// Set this object as a virtual multipath dual
void Radar::SetMultipathDual(rsFloat reflect)
{
  multipath_dual = true;
  multipath_reflect = reflect;
  // Sanity check the reflectance factor
  if (multipath_reflect > 1)
    rsDebug::printf(rsDebug::RS_CRITICAL, "[CRITICAL] Multipath reflection factor greater than 1 (=%g) for radar %s, results are likely to be incorrect\n", reflect, GetName().c_str());
}

/// Get the reflecting factor
rsFloat Radar::MultipathDualFactor() const
{
  return multipath_reflect;
}

//
// Transmitter Implementation
//

//Default constructor for Transmitter
Transmitter::Transmitter(const Platform *platform, const std::string& name, bool pulsed):
  Radar(platform, name),
  signal(0),
  pulsed(pulsed),
  dual(0)
{
}

//Default destructor for Transmitter
Transmitter::~Transmitter()
{
  delete GetTiming();
}

/// Set the transmitter's pulse waveform
void Transmitter::SetWave(RadarSignal *wave)
{
  signal = wave;
}

// Return the number of pulses this transmitter produces over the simulation lifetime
int Transmitter::GetPulseCount() const 
{
  rsFloat time = rsParameters::end_time() - rsParameters::start_time();
  if (pulsed) {
    rsFloat pulses = time*prf;
    return static_cast<int>(std::ceil(pulses));
  }
  else
    return 1; //CW systems only have one 'pulse'
}

// Fill the structure with the number'th pulse in the transmitter's pulse list
void Transmitter::GetPulse(TransmitterPulse *pulse, int number) const
{
  
  //Pulse waveform is same as transmitter waveform
  pulse->wave = signal;
  //Calculate start time of pulse
  if (pulsed) 
    pulse->time = static_cast<rsFloat>(number)/prf; //Pulse mode start depends on PRF
  else
    pulse->time = 0; //CW transmitters start at zero for now
  //If there is timing jitter, add it
  if (!timing)
    throw std::logic_error("[BUG] Transmitter "+GetName()+" must be associated with timing source");
  //pulse->time = pulse->time;//+timing->GetPulseTimeError();
  
}

/// Set the Pulse Repetition Frequency of the transmitter
void Transmitter::SetPRF(rsFloat mprf)
{
  rsFloat rate = rsParameters::rate()*rsParameters::oversample_ratio();
  // The PRF must be rounded to an even number of samples
  prf = 1/(std::floor(rate/mprf)/rate);
}

//
// Receiver Implementation
//

//Default constructor for Receiver
Receiver::Receiver(const Platform *platform, std::string name):
  Radar(platform, name),
  noise_temperature(0),
  dual(0),
  flags(0)
{
}

//Default destructor for Receiver
Receiver::~Receiver()
{
  ClearResponses();
  delete timing; //The timing is unique to the receiver
}

//Add a response to the list of responses for this receiver
void Receiver::AddResponse(Response *response)
{
  boost::try_mutex::scoped_lock lock(responses_mutex);
  responses.push_back(response);
}

//Clear the list of system responses
void Receiver::ClearResponses()
{
  std::vector<Response *>::iterator i;
  for (i = responses.begin(); i != responses.end(); i++)
    delete *i;
  responses.clear();
}

/// Comparison function for response*
inline bool CompareTimes(const Response *a, const Response *b)
{
  return (a->StartTime())<(b->StartTime());
}

/// Render the antenna's responses
void Receiver::Render()
{
  try {
    // This mutex should never be locked, enforce that condition
    boost::try_mutex::scoped_try_lock lock(responses_mutex);    
    //Sort the returns into time order
    std::sort(responses.begin(), responses.end(), CompareTimes);
    //Export the pulse descriptions to XML
    if (rsParameters::export_xml())
      ExportReceiverXML(responses, GetName() + "_results");
    //Export a binary containing the pulses
    if (rsParameters::export_binary())
      ExportReceiverBinary(responses, this, GetName(), GetName()+"_results");
    //Export to CSV format
    if (rsParameters::export_csv())
      ExportReceiverCSV(responses, GetName()+"_results");
    //Unlock the mutex
    lock.unlock();
  }
  catch (boost::lock_error &e)
    {
      throw std::runtime_error("[BUG] Responses lock is locked during Render()");
    }
}

/// Get the noise temperature (including antenna noise temperature)
rsFloat Receiver::GetNoiseTemperature(const SVec3 &angle) const
{
  return noise_temperature+Radar::GetNoiseTemperature(angle);
}

/// Get the receiver noise temperature
rsFloat Receiver::GetNoiseTemperature() const
{
  return noise_temperature;
}

/// Set the noise temperature of the receiver
void Receiver::SetNoiseTemperature(rsFloat temp)
{
  if (temp < -std::numeric_limits<rsFloat>::epsilon())
    throw std::runtime_error("Noise temperature set to negative value.");
  noise_temperature = temp;
}

/// Set the length of the receive window
void Receiver::SetWindowProperties(rsFloat length, rsFloat prf, rsFloat skip)
{
  rsFloat rate = rsParameters::rate()*rsParameters::oversample_ratio();
  window_length = length;
  window_prf = prf;
  window_skip = skip;
  // The PRF and skip must be rounded to an even number of samples
  window_prf = 1/(std::floor(rate/window_prf)/rate);
  window_skip = std::floor(rate*window_skip)/rate;
}

/// Return the number of responses
int Receiver::CountResponses() const
{
  return responses.size();
}

/// Get the number of receive windows in the simulation time
int Receiver::GetWindowCount() const
{
  rsFloat time = rsParameters::end_time() - rsParameters::start_time();
  rsFloat pulses = time*window_prf;
  return static_cast<int>(std::ceil(pulses));
}

/// Get the start time of the next window
rsFloat Receiver::GetWindowStart(int window) const
{
  //Calculate start time of pulse
  rsFloat stime = static_cast<rsFloat>(window)/window_prf+window_skip;
  //If there is timing jitter, add it
  if (!timing)
    throw std::logic_error("[BUG] Receiver must be associated with timing source");
  //stime = stime;//+timing->GetPulseTimeError();
  return stime;
}

/// Get the length of the receive window
rsFloat Receiver::GetWindowLength() const
{
  return window_length;
}

/// Get the time skipped before the start of the receive window
rsFloat Receiver::GetWindowSkip() const
{
  return window_skip;
}

/// Get the length of the receive window
rsFloat Receiver::GetPRF() const
{
  return window_prf;
}

/// Set a flag
void Receiver::SetFlag(RecvFlag flag)
{
  flags |= flag;
}

/// Check if a flag is set
bool Receiver::CheckFlag(RecvFlag flag) const
{
  return flags & flag;
}


//
// Multipath dual functions
//

// Create a multipath dual of the given receiver
Receiver* rs::CreateMultipathDual(Receiver *recv, const MultipathSurface *surf)
{
  //If we already have a dual, simply return the pointer to it
  if (recv->dual)
    return recv->dual;
  //Get the dual platform
  Platform *dual_plat = CreateMultipathDual(recv->GetPlatform(), surf);
  //Create a new receiver object
  Receiver *dual = new Receiver(dual_plat, recv->GetName()+"_dual");
  //Assign the new receiver object to the current object
  recv->dual = dual;
  //Copy data from the Radar object
  dual->antenna = recv->antenna;
  if (recv->attached)
    dual->attached = CreateMultipathDual(dynamic_cast<Transmitter*>(const_cast<Radar*>(recv->attached)), surf);
  dual->SetMultipathDual(surf->GetFactor());
  //Copy data from the receiver object
  dual->noise_temperature = recv->noise_temperature;
  dual->window_length = recv->window_length;
  dual->window_prf = recv->window_prf;
  dual->window_skip = recv->window_skip;
  dual->timing = recv->timing;
  //Done, return the created object
  return dual;  
}

// Create a multipath dual of the given transmitter
Transmitter* rs::CreateMultipathDual(Transmitter *trans, const MultipathSurface *surf)
{
  //If we already have a dual, simply return a pointer to it
  if (trans->dual)
    return trans->dual;
  //Get the dual platform
  Platform* dual_plat = CreateMultipathDual(trans->GetPlatform(), surf);
  //Create a new transmitter object
  Transmitter *dual = new Transmitter(dual_plat, trans->GetName()+"_dual", trans->pulsed);
  //Assign the the transmitter object to the current object
  trans->dual = dual;
  //Copy data from the Radar object
  dual->antenna = trans->antenna;
  if (trans->attached)
    dual->attached = CreateMultipathDual(dynamic_cast<Receiver*>(const_cast<Radar*>(trans->attached)), surf);
  dual->SetMultipathDual(surf->GetFactor());
  //Copy data from the transmitter object
  dual->prf = trans->prf;
  dual->pulsed = trans->pulsed;
  dual->signal = trans->signal;
  dual->timing = trans->timing;
  //Done, return the created object
  return dual;
}
