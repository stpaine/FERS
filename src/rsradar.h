//rsradar.h
/// Defines classes for receivers, transmitters and antennas
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//Started: 21 April 2006
#ifndef __RSRADAR_H
#define __RSRADAR_H

#include <config.h>
#include <vector>
#include <string>
#include "rsobject.h"
#include "rsradarwaveform.h"
#include <boost/utility.hpp>
#include <boost/thread/mutex.hpp>

namespace rs {

  //Forward declaration of Response (see rsresponse.h)
  class Response;

  //Forward declaration of Antenna (see rsantenna.h)
  class Antenna;
  
  //Forward declaration of SVec3 and Vec3 (see rspath.h)
  class SVec3;
  class Vec3;

  //Forward declaration of Timing (see rstiming.h)
  class Timing;

  //Forward declaration of Receiver and Transmitter (this file)
  class Receiver;
  class Transmitter;

  //Forward declaration of MultipathSurface (rsmultipath.h)
  class MultipathSurface;

  /// A pulse sent out by a transmitter
  struct TransmitterPulse
  {
    rs::RadarSignal *wave; //!< Base radar waveform
    rsFloat time; //!< The start time of the pulse
  };

  /// Base class for transmitter and receiver
  class Radar: public Object
    {
    public:
      /// Constructor
      Radar(const Platform *platform, const std::string& name);
      /// Destructor
      virtual ~Radar();
      /// Set the transmitter's antenna
      void SetAntenna(Antenna *ant);
      /// Return the antenna gain in the specified direction
      rsFloat GetGain(const SVec3& angle, const SVec3& refangle, rsFloat wavelength) const;
      /// Get the noise temperature (including antenna noise temperature)
      virtual rsFloat GetNoiseTemperature(const SVec3& angle) const;
      /// Make the radar monostatic
      void MakeMonostatic(Radar *recv);
      /// Get the attached radar, if we are monostatic
      const Radar* GetAttached() const;
      /// Return whether the radar has an attached receiver/transmitter
      bool IsMonostatic() const;
      /// Attach a timing object to the radar
      void SetTiming(Timing* tim);
      /// Get a pointer to the radar's timing object
      Timing* GetTiming() const;
      /// Is this object a virtual multipath dual?
      bool IsMultipathDual() const;
      /// Set this object as a virtual multipath dual
      void SetMultipathDual(rsFloat reflect);
      /// Get the reflecting factor
      rsFloat MultipathDualFactor() const;
    protected:
      Timing *timing; //!< The radar's timing source
    private:      
      const Antenna *antenna; //!< The radar's antenna
      const Radar *attached; //!< Other radar which shares antenna (0 if not monostatic)
      bool multipath_dual; //!< This is a virtual radar which exists for multipath simulation
      rsFloat multipath_reflect; //!< The fraction of signal power which is reflected by the multipath surface
      /// Functions which create multipath duals
      friend Receiver* CreateMultipathDual(Receiver *recv, const MultipathSurface *surf);
      friend Transmitter* CreateMultipathDual(Transmitter *trans, const MultipathSurface *surf);
    };

  //Forward declaration of Receiver
  class Receiver;
  
  /// A single radar transmitter, base class for pulsed and CW transmitters
  class Transmitter: public Radar
    {
    public:
      /// Constructor
      Transmitter(const Platform *platform, const std::string& name, bool pulsed);
      /// Destructor
      virtual ~Transmitter();
      /// Set the transmitter's pulse waveform
      void SetWave(RadarSignal *pulse);
      /// Return the number of pulses which will be sent out over the course of the simulation
      int GetPulseCount() const;
      /// Fill the pulse structure with the specified pulse
      void GetPulse(TransmitterPulse *pulse, int number) const;
      /// Set the Pulse Repetition Frequency of the transmitter
      void SetPRF(rsFloat mprf);
    protected:
      rs::RadarSignal *signal; //!< Waveform of transmitted pulses
      rsFloat prf; //!< Transmitter pulse repetition frequency (PRF)
      bool pulsed; //!< Is this a pulsed transmitter?
      Transmitter *dual; //!< Multipath dual of this transmitter
      /// Function to create multipath duals
      friend Transmitter* CreateMultipathDual(Transmitter *trans, const MultipathSurface *surf);
    };

  //Receiver objects contain receivers - the system that receives pulses
  //Receivers model the effects of their antenna
  //They also play a role in the simulation by keeping a list of responses (pulses they
  //recieve) which is built during the first part of the simulation and used to render
  //the result during the second phase.
  class Receiver: public Radar
    {
    public:
      /// Enum type for the receiver behaviour flags
      enum RecvFlag { FLAG_NODIRECT = 1, FLAG_NOPROPLOSS = 2 };
      /// Constructor
      Receiver(const Platform *platform, std::string name = "defRecv");
      virtual ~Receiver();
      /// Add a response to the list of simulation responses
      void AddResponse(Response *response);
      /// Clear the list of simulation responses (deleting responses)
      void ClearResponses();
      /// Export the results of the simulation to file
      void Render();
      /// Get the noise temperature (including antenna noise temperature)
      rsFloat GetNoiseTemperature(const SVec3 &angle) const;
      /// Get the receiver noise temperature
      rsFloat GetNoiseTemperature() const;
      /// Set the noise temperature of the receiver
      void SetNoiseTemperature(rsFloat temp);
      /// Set the length of the receive window
      void SetWindowProperties(rsFloat length, rsFloat prf, rsFloat skip);
      /// Return the number of responses
      int CountResponses() const;
      /// Get the number of receive windows in the simulation time
      int GetWindowCount() const;
      /// Get the start time of the next window
      rsFloat GetWindowStart(int window) const;
      /// Get the length of the receive window
      rsFloat GetWindowLength() const;
      /// Get the time skipped before the start of the receive window
      rsFloat GetWindowSkip() const;
      /// Get the PRF
      rsFloat GetPRF() const;
      /// Set a flag
      void SetFlag(Receiver::RecvFlag flag);
      /// Check if a flag is set
      bool CheckFlag(Receiver::RecvFlag flag) const;
    private:
      /// Vector to hold all the system responses
      std::vector<Response *> responses;
      boost::try_mutex responses_mutex; //!< Mutex to serialize access to responses
      rsFloat noise_temperature; //!< Noise temperature of the receiver
      rsFloat window_length; //!< Length of the receive window (seconds)
      rsFloat window_prf; //!< Window repetition frequency
      rsFloat window_skip; //!< The amount of time at the beginning of an interval to skip before capturing response
      Receiver *dual; //!< Multipath dual of this receiver
      int flags; //!< Flags which control receiver behaviour
      /// Function to create dual
      friend Receiver* CreateMultipathDual(Receiver *recv, const MultipathSurface *surf);
    };

  // Functions for duplicating receivers and transmitters to create duals for multipath simulation
  Receiver* CreateMultipathDual(Receiver *recv, const MultipathSurface *surf);
  Transmitter* CreateMultipathDual(Transmitter *trans, const MultipathSurface *surf);

}

#endif
