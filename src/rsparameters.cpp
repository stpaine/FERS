//rsparameters.cpp
//Implementation of Singleton class to hold common simulation parameters
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//11 June 2006

#include <stdexcept>
#include "time.h"
#include "rsparameters.h"
#include "rsdebug.h"

using namespace rs;

rsParameters* rsParameters::instance = new rsParameters();

namespace {

struct SimParameters {
  rsFloat c; //!< Propagation speed of the wave in the medium
  rsFloat start; //!< The start time of the simulation
  rsFloat end; //!< The end time of the simulation
  rsFloat cw_sample_rate; //<! The number of samples per second to take of changes in the CW state
  rsFloat rate; //!< The sample rate to use for rendering
  unsigned int random_seed; //!< The seed used for random number calculations
  unsigned int adc_bits; //!< The number of bits to use for quantization
  unsigned int filter_length; //!< The length of the filter for rendering purposes
  rsParms::BinaryFileType filetype; //!< The type of binary files produced by binary rendering
  bool export_xml; //!< Export results in XML format
  bool export_csv; //!< Export results in CSV format
  bool export_binary; //!< Export results in binary format
  unsigned int render_threads; //!< Number of threads to use to render each receiver
  unsigned int oversample_ratio; //!< Ratio of oversampling applied to pulses before rendering
};

/// Object which contains all the simulation parameters
SimParameters sim_parms;

}

//Private constructor for rsParameter, should only be called once
rsParameters::rsParameters() {
  //Default value of c, speed of light in a vacuum
  sim_parms.c = 299792458.0; 
  //Simulation defaults to zero length
  sim_parms.start = 0;
  sim_parms.end = 0;
  //CW Interpolation rate defaults to 1000 per second
  sim_parms.cw_sample_rate = 1000;
  // Oversample by default
  sim_parms.rate = 0;
  // Default filter length is 33
  sim_parms.filter_length = 33;
  // Binary file type defaults to CSV
  sim_parms.filetype = rsParms::RS_FILE_FERSBIN;
  // Export xml by default
  sim_parms.export_xml = true;
  // Export csv by default
  sim_parms.export_csv = true;
  // Don't export binary by default
  sim_parms.export_binary = false;
  // The random seed is set the to the current time by default
  sim_parms.random_seed = static_cast<unsigned int>(time(NULL));
  // The default is not to quantize
  sim_parms.adc_bits = 0;
  // Default maximum number of render threads
  sim_parms.render_threads = 1;
  // Default is to disable oversampling
  sim_parms.oversample_ratio = 1;
}

rsParameters *rsParameters::modify_parms()
{
  if (!instance)
    instance = new rsParameters;
  return instance;
}

//Getters for settings
rsFloat rsParameters::c()
{
  if (!instance)
    instance = new rsParameters;
  return sim_parms.c;
}

rsFloat rsParameters::boltzmann_k()
{
  return 1.3806503e-23l;
}

rsFloat rsParameters::start_time()
{
  if (!instance)
    instance = new rsParameters;
  return sim_parms.start;
}

rsFloat rsParameters::end_time()
{
  if (!instance)
    instance = new rsParameters;
  return sim_parms.end;
}

rsFloat rsParameters::cw_sample_rate()
{
  if (!instance)
    instance = new rsParameters;
  return sim_parms.cw_sample_rate;
}

rsParms::BinaryFileType rsParameters::binary_file_type()
{
  if (!instance)
    instance = new rsParameters;
  return sim_parms.filetype;
}

rsFloat rsParameters::rate()
{
  if (!instance)
    instance = new rsParameters;
  return sim_parms.rate;
}

unsigned int rsParameters::random_seed()
{
  if (!instance)
    instance = new rsParameters;
  return sim_parms.random_seed;
}

unsigned int rsParameters::adc_bits()
{
  if (!instance)
    instance = new rsParameters();
  return sim_parms.adc_bits;
}

bool rsParameters::export_xml()
{
  if (!instance)
    instance = new rsParameters;
  return sim_parms.export_xml;
}

bool rsParameters::export_csv()
{
  if (!instance)
    instance = new rsParameters;
  return sim_parms.export_csv;
}

bool rsParameters::export_binary()
{
  if (!instance)
    instance = new rsParameters;
  return sim_parms.export_binary;
}

/// Length to use for the rendering filter
unsigned int rsParameters::render_filter_length()
{
  if (!instance)
    instance = new rsParameters;
  return sim_parms.filter_length;
}

/// Maximum number of threads to use for rendering
unsigned int rsParameters::render_threads() {
  if (!instance)
    instance = new rsParameters;
  return sim_parms.render_threads;
}

unsigned int rsParameters::oversample_ratio()
{
  if (!instance)
    instance = new rsParameters;
  return sim_parms.oversample_ratio;
}

//
//Setters for global parameters
//

void rsParameters::SetC(rsFloat c)
{
  sim_parms.c = c;
  rsDebug::printf(rsDebug::RS_CRITICAL, "[CRITICAL] Propagation speed (c) set to custom value: %8.5f\n", c);
}

void rsParameters::SetTime(rsFloat start, rsFloat end)
{
  sim_parms.start = start;
  sim_parms.end = end;
}

void rsParameters::SetCWSampleRate(rsFloat rate)
{
  sim_parms.cw_sample_rate = rate;
}

void rsParameters::SetRate(rsFloat factor)
{
  sim_parms.rate = factor;
  rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "[VV] System sample rate set to custom value: %8.5f\n", factor);
}

void rsParameters::SetRandomSeed(unsigned int random_seed)
{
  sim_parms.random_seed = random_seed;
}

void rsParameters::SetBinaryFileType(rsParms::BinaryFileType type)
{
  sim_parms.filetype = type;
}

void rsParameters::SetExporters(bool xml, bool csv, bool binary)
{
  sim_parms.export_xml = xml;
  sim_parms.export_csv = csv;
  sim_parms.export_binary = binary;
}

void rsParameters::SetADCBits(unsigned int bits)
{
  sim_parms.adc_bits = bits;
}

void rsParameters::SetRenderFilterLength(unsigned int length)
{
  //Sanity check the render filter length
  if (length < 16)
    throw std::runtime_error("[ERROR] Render filter length must be > 16"); 
  sim_parms.filter_length = length;
  rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "[VV] Render filter length set to custom value: %d\n", length);
}

void rsParameters::SetOversampleRatio(unsigned int ratio)
{
  //Sanity check the ratio
  if (ratio == 0)
    throw std::runtime_error("[ERROR] Oversample ratio must be >= 1");
  sim_parms.oversample_ratio = ratio;
  rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "[VV] Oversampling enabled with ratio %d\n", ratio);
}

void rsParameters::SetThreads(unsigned int threads)
{
  sim_parms.render_threads = threads;
}
