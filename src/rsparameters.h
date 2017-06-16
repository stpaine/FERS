//rsparameters.h
//Singleton class to hold common simulation parameters
// The parameters system holds all global simulation parameters, magic numbers and other global values
//No 'magic numbers' (such as values of c) are to be used in the code - store them here instead
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//11 June 2006

#ifndef __RS_PARAMETERS_H
#define __RS_PARAMETERS_H
#include <config.h>

namespace rsParms {
  enum BinaryFileType {RS_FILE_CSV, RS_FILE_FERSBIN, RS_FILE_RAW};
}

namespace rs {
///'Singleton' class to hold parameters for the simulation
class rsParameters {
public:
  /// Method to return a pointer to the single instance of the class that can be used to modify parameters
  static rsParameters* modify_parms();
  /// Get the value of c (propagation speed in the medium)
  static rsFloat c();
  /// Get the value of boltzmann's constant
  static rsFloat boltzmann_k();
  /// Get the value of start (start time of the simulation)
  static rsFloat start_time();
  /// Get the value of end (end time of the simulation)
  static rsFloat end_time();
  /// Get the CW interpolation sample rate
  static rsFloat cw_sample_rate();
  /// Get the oversample factor
  static rsFloat rate();
  /// Get the current random seed
  static unsigned int random_seed();
  /// Get the number of ADC bits used for quantization
  static unsigned int adc_bits();
  /// Get the binary file type
  static rsParms::BinaryFileType binary_file_type();
  /// Do we export in XML format?
  static bool export_xml();
  /// Do we export in CSV format?
  static bool export_csv();
  /// Do we export in binary format?
  static bool export_binary();
  /// Length to use for the rendering filter
  static unsigned int render_filter_length();
  /// Maximum number of threads to use for rendering
  static unsigned int render_threads();
  /// Number of times to oversample loaded pulses before simulation
  static unsigned int oversample_ratio();
  /// Set the value of c
  void SetC(rsFloat c);
  /// Set the start and end times
  void SetTime(rsFloat start, rsFloat end);
  /// Set the CW sample rate
  void SetCWSampleRate(rsFloat rate);
  /// Set the export sample rate
  void SetRate(rsFloat factor);
  /// Set the random seed
  void SetRandomSeed(unsigned int random_seed);
  /// Set the binary file type
  void SetBinaryFileType(rsParms::BinaryFileType type);
  /// Set the enabled exporters
  void SetExporters(bool xml, bool csv, bool binary);
  /// Set the number of bits used for quantization
  void SetADCBits(unsigned int bits);
  /// Set the render filter length
  void SetRenderFilterLength(unsigned int length);
  /// Set the number of times to oversample loaded pulses before simulation
  void SetOversampleRatio(unsigned int ratio);
  /// Set the number of threads to use
  void SetThreads(unsigned int threads);
protected:
  /// The default constructor is private
  rsParameters();
  /// The copy constructor is private
  rsParameters(const rsParameters &rs);
  /// The assignment operator is private
  rsParameters &operator=(const rsParameters &rs);
  /// Pointer to a single instance of the class
  static rsParameters *instance;
};

}
#endif
