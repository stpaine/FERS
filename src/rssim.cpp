//rssim.cpp
//Functions which perform the actual simulations
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//30 May 2006

#include <cmath>
#include <limits>
#include <stdexcept>
#include "rssim.h"
#include "rsworld.h"
#include "rsradar.h"
#include "rsdebug.h"
#include "rsparameters.h"
#include "rsantenna.h"
#include "rstarget.h"
#include "rsresponse.h"
#include "rsnoise.h"

using namespace rs;

namespace {

/// Results of solving the bistatic radar equation and friends
struct REResults {
  rsFloat power;
  rsFloat delay;
  rsFloat doppler;
  rsFloat phase;
  rsFloat noise_temperature;
};

/// Class for range errors in RE calculations
class RangeError {
};

///Solve the radar equation for a given set of parameters
void SolveRE(const Transmitter *trans, const Receiver *recv, const Target *targ, rsFloat time, rsFloat length, RadarSignal *wave, REResults &results)
{
  //Get the positions in space of the three objects
  Vec3 trpos = trans->GetPosition(time);
  Vec3 repos = recv->GetPosition(time);
  Vec3 tapos = targ->GetPosition(time);
  SVec3 transvec = SVec3(tapos-trpos);
  SVec3 recvvec = SVec3(tapos-repos); 
  //Calculate the distances
  rsFloat Rt = transvec.length;
  rsFloat Rr = recvvec.length;
  // From here on, transvec and recvvec need to be normalized
  transvec.length = 1;
  recvvec.length = 1;
  //Sanity check Rt and Rr and throw an exception if they are too small
  if ((Rt <= std::numeric_limits<rsFloat>::epsilon()) || (Rr <= std::numeric_limits<rsFloat>::epsilon()))
    throw RangeError();
  //Step 1, calculate the delay (in seconds) experienced by the pulse
  //See "Delay Equation" in doc/equations/equations.tex
  results.delay = (Rt+Rr)/rsParameters::c();
  //rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Calculated delay as %fs, %fm\n", results.delay, results.delay * rsParameters::c());
  //Get the RCS
  rsFloat RCS = targ->GetRCS(transvec, recvvec);
  //Get the wavelength
  rsFloat Wl = rsParameters::c()/wave->GetCarrier();
  //Get the system antenna gains (which include loss factors)
  rsFloat Gt = trans->GetGain(transvec, trans->GetRotation(time), Wl);
  rsFloat Gr = recv->GetGain(recvvec, recv->GetRotation(results.delay+time), Wl);
    //  rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Gt: %e Gr: %e\n", Gt, Gr);
  //Step 2, calculate the received power using the narrowband bistatic radar equation
  //See "Bistatic Narrowband Radar Equation" in doc/equations/equations.tex  
  results.power = Gt*Gr*RCS/(4*M_PI);
  if (!recv->CheckFlag(Receiver::FLAG_NOPROPLOSS))
    results.power *= (Wl*Wl)/(pow(4*M_PI, 2)*Rt*Rt*Rr*Rr);
  //   rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Pr: %2.9e Rt: %e Rr: %e Gt: %e Gr: %e RCS: %e Wl %e\n", results.power, Rt, Rr, Gt, Gr, RCS, Wl);
  // If the transmitter and/or receiver are multipath duals, multiply by the loss factor
  if (trans->IsMultipathDual())
    results.power *= trans->MultipathDualFactor();
  if (recv->IsMultipathDual())
    results.power *= trans->MultipathDualFactor();
  //Step 3, calculate phase shift
  //See "Phase Delay Equation" in doc/equations/equations.tex
  results.phase = -results.delay*2*M_PI*wave->GetCarrier();
  //Step 4, calculate doppler shift
  //Calculate positions at the end of the pulse
  Vec3 trpos_end = trans->GetPosition(time+length);
  Vec3 repos_end = recv->GetPosition(time+length);
  Vec3 tapos_end = targ->GetPosition(time+length);
  SVec3 transvec_end = SVec3(tapos_end-trpos_end);
  SVec3 recvvec_end = SVec3(tapos_end-repos_end); 
  rsFloat Rt_end = transvec_end.length;
  rsFloat Rr_end = recvvec_end.length;
  //Sanity check Rt_end and Rr_end and throw an exception if they are too small
  if ((Rt_end < std::numeric_limits<rsFloat>::epsilon()) || (Rr_end < std::numeric_limits<rsFloat>::epsilon()))
    throw std::runtime_error("Target is too close to transmitter or receiver for accurate simulation");
  //Doppler shift equation
  //See "Bistatic Doppler Equation" in doc/equations/equations.tex
  rsFloat V_r = (Rr_end-Rr)/length;
  rsFloat V_t = (Rt_end-Rt)/length;
  results.doppler = std::sqrt((1-V_r/rsParameters::c())/(1+V_r/rsParameters::c()))*std::sqrt((1-V_t/rsParameters::c())/(1+V_t/rsParameters::c()));
  //Step 5, calculate system noise temperature
  //We only use the receive antenna noise temperature for now
  results.noise_temperature = recv->GetNoiseTemperature(recv->GetRotation(time+results.delay));
}

//Perform the first stage of CW simulation calculations for the specified pulse and target
void SimulateTarget(const Transmitter *trans, Receiver *recv, const Target *targ, const World *world, const TransmitterPulse *signal) {
  //Get the simulation start and end time
  rsFloat start_time = signal->time;
  rsFloat end_time = signal->time + signal->wave->GetLength();
  //Calculate the number of interpolation points we need to add
  rsFloat sample_time = 1.0 / rsParameters::cw_sample_rate();
  rsFloat point_count = std::ceil(signal->wave->GetLength() / sample_time);
  //Create the response
  Response *response = new Response(signal->wave, trans);
  try {
    //Loop through and add interpolation points
    for (int i = 0; i < point_count; i++) {
      rsFloat stime = i * sample_time + start_time; //Time of the start of the sample
      REResults results;
      SolveRE(trans, recv, targ, stime, sample_time, signal->wave, results);
      InterpPoint point(results.power, stime + results.delay, results.delay, results.doppler, results.phase, results.noise_temperature);
      response->AddInterpPoint(point);
      //rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Added interpolation point with delay: %gs, %fm\n", results.delay, results.delay * rsParameters::c());
      //rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "stime + results.delay = : %gs\n", stime + results.delay);
    }
    //Add one more point at the end
    REResults results;
    SolveRE(trans, recv, targ, end_time, sample_time, signal->wave, results);
    InterpPoint point(results.power, end_time + results.delay, results.delay, results.doppler, results.phase, results.noise_temperature);
    response->AddInterpPoint(point);
  }
  catch (RangeError &re) {
    throw std::runtime_error("Receiver or Transmitter too close to Target for accurate simulation");
  }
  //Add the response to the receiver
  recv->AddResponse(response);
}

/// Solve the Radar Equation and friends (doppler, phase, delay) for direct transmission
void SolveREDirect(const Transmitter *trans, const Receiver *recv, rsFloat time, rsFloat length, const RadarSignal *wave, REResults &results)
{
  //Calculate the vectors to and from the transmitter
  Vec3 tpos = trans->GetPosition(time);
  Vec3 rpos = recv->GetPosition(time);
  SVec3 transvec = SVec3(tpos-rpos);
  SVec3 recvvec = SVec3(rpos-tpos);
  //Calculate the range
  rsFloat R = transvec.length;
  //Normalize transvec and recvvec for angle calculations
  transvec.length = 1;
  recvvec.length = 1;
  //If the two antennas are not in the same position, this can be calculated
  if (R > std::numeric_limits<rsFloat>::epsilon())
    {

      //Step 1: Calculate the delay
      //See "Delay Equation" in doc/equations/equations.tex
      results.delay = R/rsParameters::c();
      
      //Calculate the wavelength
      rsFloat Wl = rsParameters::c()/wave->GetCarrier();
      //Get the antenna gains
      rsFloat Gt = trans->GetGain(transvec, trans->GetRotation(time), Wl);
      rsFloat Gr = recv->GetGain(recvvec, recv->GetRotation(time+results.delay), Wl);
      //Step 2: Calculate the received power 
      //See "One Way Radar Equation" in doc/equations/equations.tex
      results.power = Gt*Gr*Wl*Wl/(4*M_PI);
      if (!recv->CheckFlag(Receiver::FLAG_NOPROPLOSS))
        results.power *= 1/(4*M_PI*pow(R, 2));
      //  rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Pr: %2.9e R: %e Gt: %e Gr: %e Wl %e\n", results.power, R, Gt, Gr, Wl);
      //Step 3: Calculate the doppler shift (if one of the antennas is moving)
      Vec3 tpos_end = trans->GetPosition(time+length);
      Vec3 rpos_end = recv->GetPosition(time+length);
      Vec3 trpos_end = tpos_end - rpos_end;
      rsFloat R_end = trpos_end.Length();
      //Calculate the Doppler shift
      //See "Monostatic Doppler Equation" in doc/equations/equations.tex
      rsFloat vdoppler = (R_end-R)/length;
      results.doppler = (rsParameters::c()+vdoppler)/(rsParameters::c()-vdoppler);
      // For the moment, we don't handle direct paths for multipath duals, until we can do it correctly
      if (trans->IsMultipathDual())
	results.power = 0;
      if (recv->IsMultipathDual())
	results.power = 0;

      //Step 4, calculate phase shift
      //See "Phase Delay Equation" in doc/equations/equations.tex
      results.phase = fmod(results.delay*2*M_PI*wave->GetCarrier(), 2*M_PI);

      //Step 5, calculate noise temperature
      results.noise_temperature = recv->GetNoiseTemperature(recv->GetRotation(time+results.delay));
    }
  else {
    //If the receiver and transmitter are too close, throw a range error
    throw RangeError();
  }
}

/// Model the pulse which is received directly by a receiver from a CW transmitter
void AddDirect(const Transmitter *trans, Receiver *recv, const World *world, const TransmitterPulse *signal) {
  //If receiver and transmitter share the same antenna - there can't be a direct pulse
  if (trans->IsMonostatic() && (trans->GetAttached() == recv))
    return;
  //Get the simulation start and end time
  rsFloat start_time = signal->time;
  rsFloat end_time = signal->time+signal->wave->GetLength();
  //Calculate the number of interpolation points we need to add
  rsFloat sample_time = 1.0/rsParameters::cw_sample_rate();
  rsFloat point_count = std::ceil(signal->wave->GetLength()/sample_time);
  //Create the CW response
  Response *response = new Response(signal->wave, trans);
  try {
    //Loop through and add interpolation points
    for (int i = 0; i < point_count; i++) {
      rsFloat stime = i*sample_time+start_time;
      REResults results;
      SolveREDirect(trans, recv, stime, sample_time, signal->wave, results);
      InterpPoint point(results.power, results.delay+stime, results.delay, results.doppler, results.phase, results.noise_temperature);
      response->AddInterpPoint(point);
    }
    //Add one more point at the end
    REResults results;
    SolveREDirect(trans, recv, end_time, sample_time, signal->wave, results);
    InterpPoint point(results.power, results.delay+end_time, results.delay, results.doppler, results.phase, results.noise_temperature);
    response->AddInterpPoint(point);

  }
  catch (RangeError &re) {
    throw std::runtime_error("Receiver or Transmitter too close to Target for accurate simulation");
  }
  //Add the response to the receiver
  recv->AddResponse(response);
}

}

//Simulate a transmitter-receiver pair with a pulsed transmission
void rs::SimulatePair(const Transmitter *trans, Receiver *recv, const World *world)
{
  //Get the number of pulses
  int pulses = trans->GetPulseCount();
  std::vector<Target*>::const_iterator targ;
  //Build a pulse
  TransmitterPulse* pulse = new TransmitterPulse();
  rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "%d\n", pulses);
  //Loop throught the pulses
  for (int i = 0; i < pulses; i++) {
      trans->GetPulse(pulse, i);
      for (targ = world->targets.begin(); targ != world->targets.end(); targ++) {
	SimulateTarget(trans, recv, *targ, world, pulse);
      }

      // Check if direct pulses are being considered for this receiver
      if (!recv->CheckFlag(Receiver::FLAG_NODIRECT)) {
        //Add the direct pulses
        AddDirect(trans, recv, world, pulse);
      }
    }
  delete pulse;
}
