//rstiming.cpp - Implementation of timing sources

#include "rstiming.h"
#include "rsnoise.h"
#include "rsdebug.h"
#include <algorithm>

using namespace rs;

//
// Timing Implementation
//

/// Default constructor
Timing::Timing(const std::string &name):
  name(name)
{
}

/// Destructor
Timing::~Timing()
{
} 

/// Get the name of the timing source
std::string Timing::GetName() const
{
  return name;
}

//
// PrototypeTiming Implementation
//

/// Constructor
PrototypeTiming::PrototypeTiming(const std::string &name):
  name(name)
{
  freq_offset = 0;
  phase_offset = 0;
  random_phase = 0;
  random_freq = 0;
  frequency = 0;
  synconpulse = false;
}

/// Add an alpha and a weight to the timing prototype
void PrototypeTiming::AddAlpha(rsFloat alpha, rsFloat weight)
{
  alphas.push_back(alpha);
  weights.push_back(weight);
}

/// Get the alphas and weights from the prototype
void PrototypeTiming::GetAlphas(std::vector<rsFloat> &get_alphas, std::vector<rsFloat> &get_weights) const
{
  //Copy the alpha and weight vectors
  get_alphas = alphas;
  get_weights = weights;
}

/// Set a constant frequency offset
void PrototypeTiming::AddFreqOffset(rsFloat offset)
{
  if (random_freq)
    rsDebug::printf(rsDebug::RS_IMPORTANT, "[Important] Random frequency offset and constant frequency offset are set for timing source %s. Only the random offset will be used.", GetName().c_str());
  freq_offset = offset;
}

/// Set a constant phase offset
void PrototypeTiming::AddPhaseOffset(rsFloat offset)
{
  if (random_phase)
    rsDebug::printf(rsDebug::RS_IMPORTANT, "[Important] Random phase offset and constant phase offset are set for timing source %s. Only the random offset will be used.", GetName().c_str());
  phase_offset = offset;
}

/// Set a random frequency offset
void PrototypeTiming::AddRandomFreqOffset(rsFloat stdev)
{
  if (freq_offset)
    rsDebug::printf(rsDebug::RS_IMPORTANT, "[Important] Random frequency offset and constant frequency offset are set for timing source %s. Only the random offset will be used.", GetName().c_str());
  random_freq = stdev;
}

/// Set a random phase offset
void PrototypeTiming::AddRandomPhaseOffset(rsFloat stdev)
{
  if (phase_offset)
    rsDebug::printf(rsDebug::RS_IMPORTANT, "[Important] Random phase offset and constant phase offset are set for timing source %s. Only the random offset will be used.", GetName().c_str());
  random_phase = stdev;
}

/// Get the phase offset
rsFloat PrototypeTiming::GetPhaseOffset() const
{
  if (random_phase != 0)
    return rsNoise::WGNSample(random_phase);
  else
    return phase_offset;
}

/// Get the phase offset
rsFloat PrototypeTiming::GetFreqOffset() const
{
  if (random_freq != 0)
    return rsNoise::WGNSample(random_freq);
  else
    return freq_offset;
}

/// Get the frequency
rsFloat PrototypeTiming::GetFrequency() const
{
  return frequency;
}


/// Get the name of the prototype
std::string PrototypeTiming::GetName() const
{
  return name;
}

/// Set the base frequency of the clock model
void PrototypeTiming::SetFrequency(rsFloat freq) {
  frequency = freq;
}

/// Set the sync on pulse flag -- timing error resets at the start of the pulse
void PrototypeTiming::SetSyncOnPulse()
{
  synconpulse = true;
}

/// Get the value of the sync on pulse flag
bool PrototypeTiming::GetSyncOnPulse() const
{
  return synconpulse;
}

//
// ClockModelTiming Implementation
//

/// Constructor
ClockModelTiming::ClockModelTiming(const std::string &name):
  Timing(name),
  enabled(false),
  model(0)
{
}

/// Destructor
ClockModelTiming::~ClockModelTiming() {
  delete model;
}

/// Initialize the clock model generator
void ClockModelTiming::InitializeModel(const PrototypeTiming *timing)
{
  if (!alphas.empty())
    throw std::logic_error("[BUG] ClockModelTiming::InitializeModel called more than once");
  //Copy the alpha and weight vectors
  timing->GetAlphas(alphas, weights);
  rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "%d\n", alphas.size());
  //Create the generator
  model = new ClockModelGenerator(alphas, weights, timing->GetFrequency(), timing->GetPhaseOffset(), timing->GetFreqOffset(), 15);
  //Warn if frequency is not set
  if (timing->GetFrequency() == 0.0)
    rsDebug::printf(rsDebug::RS_IMPORTANT, "[Important] Timing source frequency not set, results could be incorrect.");
  //Get the carrier frequency
  frequency = timing->GetFrequency();
  // Get the sync on pulse flag
  synconpulse = timing->GetSyncOnPulse();
  //Enable the model
  enabled = true;
   
}

/// Return the enabled state of the clock model
bool ClockModelTiming::Enabled()
{
  return enabled && model->Enabled();
}

/// Get the real time of a particular pulse
rsFloat ClockModelTiming::GetPulseTimeError() const
{
  if (enabled)
    return model->GetSample();
  else
    return 0;
}

/// Skip a sample, computing only enough to preserve long term correlations
void ClockModelTiming::SkipSamples(long long samples)
{
  if (enabled)
       model->SkipSamples(samples);
}

/// Get the value of the sync on pulse flag
bool ClockModelTiming::GetSyncOnPulse() const
{
  return synconpulse;
}

/// Reset the clock phase error to zero
void ClockModelTiming::Reset()
{
  model->Reset();
}

/// Get the next sample of time error for a particular pulse
rsFloat ClockModelTiming::NextNoiseSample()
{
  if (enabled)
    return model->GetSample();
  else
    return 0;
}

/// Get the carrier frequency of the modelled clock
rsFloat ClockModelTiming::GetFrequency() const
{
  return frequency;
}
