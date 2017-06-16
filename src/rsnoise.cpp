//rsnoise.cpp - Functions for generating different types of noise
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//14 August 2006

#include <cmath>
#include <limits>
#include <boost/random.hpp>
#include "rsnoise.h"
#include "rsdebug.h"
#include "rsparameters.h"
#include "rsdsp.h"

using namespace rs;

namespace {
//Use the Mersenne Twister PRNG with parameter 19937, for why this is a good choice see 
//Mersenne Twister: A 623-dimensionally equidistributed uniform pseudo-random number generator
//Matsumoto et al.
//ACM Transactions on Modeling and Computer Simulation
//January 1998
boost::mt19937* rng;

  //Used to generate single noise samples
  boost::normal_distribution<rsFloat> nd(0,1);
  boost::uniform_real<rsFloat> ud(0, 1.0);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<rsFloat> >* normal_vg;
  boost::variate_generator<boost::mt19937&, boost::uniform_real<rsFloat> >* uniform_vg;
}

//
// Implementation of non-class functions
//

/// Initialize the random number generator code (must be called once, after the loading of the script)
void rsNoise::InitializeNoise()
{
  delete rng;
  delete normal_vg;
  rng = new boost::mt19937(rsParameters::random_seed());
  normal_vg = new boost::variate_generator<boost::mt19937&, boost::normal_distribution<rsFloat> >(*rng, nd);
  uniform_vg = new boost::variate_generator<boost::mt19937&, boost::uniform_real<rsFloat> >(*rng, ud);
}

/// Clean up the noise code
void rsNoise::CleanUpNoise()
{
  delete rng;
  delete normal_vg;
  delete uniform_vg;
}

/// Return a single sample of white gaussian noise
rsFloat rsNoise::WGNSample(rsFloat stddev)
{
  if (stddev > std::numeric_limits<rsFloat>::epsilon())
    return (*normal_vg)()*stddev;
  else
    return 0;
}

/// Return a single uniformly distributed sample in [0, 1]
rsFloat rsNoise::UniformSample()
{
  return (*uniform_vg)();
}

/// Calculate noise amplitude from the temperature
rsFloat rsNoise::NoiseTemperatureToPower(rsFloat temperature, rsFloat bandwidth)
{
  return rsParameters::boltzmann_k()*temperature*bandwidth; //See equations.tex
}

//
// NoiseGenerator Implementation
//

//Constructor
NoiseGenerator::NoiseGenerator()
{
}

//Destructor
NoiseGenerator::~NoiseGenerator()
{
}

//
// Gamma Generator
//

/// Constructor
GammaGenerator::GammaGenerator(rsFloat k):
  dist(k),
  gen(*rng, dist)
{

}

/// Destructor
GammaGenerator::~GammaGenerator()
{
}

/// Get a single random sample
rsFloat GammaGenerator::GetSample()
{
  return gen();
}

/// Operator to get a random sample
rsFloat GammaGenerator::operator()()
{
  return gen();
}
//
// WGNGenerator Implementation
//

//Constructor
WGNGenerator::WGNGenerator(rsFloat stddev)
{
  dist = boost::normal_distribution<rsFloat>(0, stddev);
  gen = new boost::variate_generator<boost::mt19937&, boost::normal_distribution<rsFloat> >(*rng, dist);
}

//Default constructor
WGNGenerator::WGNGenerator()
{
  dist = boost::normal_distribution<rsFloat>(0, 1);
  gen = new boost::variate_generator<boost::mt19937&, boost::normal_distribution<rsFloat> >(*rng, dist);
}

// Destructor
WGNGenerator::~WGNGenerator()
{
  delete gen;
}

// Get a sample from the rng
rsFloat WGNGenerator::GetSample()
{
  return (*gen)();
}

//
// FAlphaBranch Implementation
//

/// Constructor
FAlphaBranch::FAlphaBranch(rsFloat ffrac, unsigned int fint, FAlphaBranch *pre, bool last):
  pre(pre),
  last(last),
  ffrac(ffrac),
  fint(fint)
{

  rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "[VV] Making branch ffrac=%f fint=%d\n", ffrac, fint);
  //Calculate scale for upsampling
  upsample_scale = std::pow(10, ffrac+fint+0.5);
  //Initialize the filters for shaping, highpass and upsampling
  Init();
  // Create a buffer for ten samples
  buffer = new rsFloat[10];
  if (!last)
    Refill();
}

/// Destructor
FAlphaBranch::~FAlphaBranch() {
  delete pre;
  Clean();
}

///Initialize the branch filters
void FAlphaBranch::Init() {
  shape_filter = 0;
  integ_filter = 0;
  highpass = 0;
  //Create the upsampler
  upsampler = new DecadeUpsampler();
  if (pre) {
    /// Numerator coefficients for elliptical highpass
    const rsFloat hp_num[12] = { 3.817871081981451e-01,
				 -4.093384095523618e+00,
				 2.005300512623078e+01,
				 -5.924672881811163e+01,
				 1.172948159891025e+02,
				 -1.633810410083022e+02,
				 1.633810410083034e+02,
				 -1.172948159891052e+02,
				 5.924672881811390e+01,
				 -2.005300512623186e+01,
				 4.093384095523903e+00,
				 -3.817871081981776e-01};
    /// Denominator coefficients for elliptical highpass
    const rsFloat hp_den[12] = {  1.000000000000000e+00,
				  -8.829695665523831e+00,
				  3.583068809011030e+01,
				  -8.811479652970442e+01,
				  1.457874067329429e+02,
				  -1.702715637111961e+02,
				  1.431504350055831e+02,
				  -8.656925883534657e+01,
				  3.687395592491803e+01,
				  -1.052413841411803e+01,
				  1.808292123637038e+00,
				  -1.412932578340511e-01};
    //Initialize the highpass filter
    highpass = new IIRFilter(hp_den, hp_num, 12);
  }
  // Initialize the shaping filter
  if (ffrac == 0.5) {
    /// Numerator co-efficients for 1/f^0.5 rolloff
    const rsFloat sf_num[16] = { 5.210373977738306e-03,
				 -7.694671394585578e-03,
				 1.635979377907092e-03,
				 9.852449140857658e-05,
				 -2.080553126780113e-03,
				 4.088764157029523e-03,
				 -1.549082440084623e-03,
				 9.054734252370680e-04,
				 -3.467369912368729e-04,
				 4.516383087838856e-04,
				 -1.063356106118517e-03,
				 1.330008998057684e-04,
				 6.556909567323943e-04,
				 -4.839476350293955e-04,
				 6.664936170526832e-05,
				 1.528520559763056e-05};
    const rsFloat sf_den[16] = {     1.000000000000000e+00,
				     -2.065565041154101e+00,
				     1.130909190864681e+00,
				     -1.671244644503288e-01,
				     -3.331474931013877e-01,
				     9.952625337612708e-01,
				     -7.123036343635182e-01,
				     3.297062696290504e-01,
				     -1.925691520710595e-01,
				     1.301247006176314e-01,
				     -2.702016290409912e-01,
				     1.455380885858886e-01,
				     1.091921868353888e-01,
				     -1.524953111510459e-01,
				     5.667716332023935e-02,
				     -2.890314873767405e-03};
    //Gain of shaping filter
    shape_gain = 5.210373977738306e-03;
    //Create the shaping filter
    shape_filter = new IIRFilter(sf_den, sf_num, 16);
  }
  else if (ffrac == 0) {
    shape_filter = 0;
  }
  else {
    rsDebug::printf(rsDebug::RS_CRITICAL, "[CRITICAL] Value of ffrac is %f\n", ffrac);
    throw std::runtime_error("Fractional integrator values other than 0.5 not currently supported");
  }
  //Initialize the integration filters
  if (fint > 0) {
    //Gain of integration filter
    integ_gain = 1;//2.755923548698047e-05;
    if (fint == 1) {
      const rsFloat i_den[2] = {1, -1};
      const rsFloat i_num[2] = {1, 0};
      integ_filter = new IIRFilter(i_den, i_num, 2);
    }
    if (fint == 2) {
      const rsFloat i_den[3] = {1, -2, 1};
      const rsFloat i_num[3] = {1, 0, 0};
      integ_filter = new IIRFilter(i_den, i_num, 3);
    }
    if (fint > 2)
      throw std::runtime_error("Only alpha values between 2 and -2 are supported for noise generation");
  }
  //Initialize the offset
  offset_sample = 0;
  got_offset = false;
  // Create a buffer for ten samples
  buffer = new rsFloat[10];
  if (!last)
    Refill();
  pre_scale = 1;
}

/// Get a sample from the branch
rsFloat FAlphaBranch::GetSample() {
  if (!last) {
    rsFloat ret = buffer[buffer_samples];
    buffer_samples++;
    if (buffer_samples == 10)
      Refill();
    return ret;
  }
  else {
    return CalcSample()+offset_sample*upsample_scale;
  }
}

///Clean up the filters, etc
void FAlphaBranch::Clean() {
  delete highpass;
  delete[] buffer;
  delete integ_filter;
  delete shape_filter;
  delete upsampler;
}

/// Calculate a single sample
rsFloat FAlphaBranch::CalcSample() {

  rsFloat sample = rsNoise::WGNSample(1);
  if (shape_filter)
    sample = shape_filter->Filter(sample)/shape_gain;
  if (integ_filter)
    sample = integ_filter->Filter(sample)/integ_gain;
  if (pre) {
    // Apply highpass only if we have branches below us
    sample = highpass->Filter(sample);    
    // If there is a branch below us, add a sample from that
    if (got_offset) {
      sample += pre->GetSample()*pre_scale-offset_sample;
    }   
    else {
      got_offset = true;
      offset_sample = pre->GetSample()*pre_scale;
    }
  } 
  return sample;
}

/// Refill the buffer
void FAlphaBranch::Refill() {
  rsFloat sample = CalcSample();
  // Pass the sample to the upsampler
  upsampler->Upsample(sample, buffer);  
  // Scale the buffer
  for (int i = 0; i < 10; i++) {
     buffer[i] *= upsample_scale;
     buffer[i] += offset_sample;
  }
  // Reset the sample counter
  buffer_samples = 0;
}

/// Refill the buffer
void FAlphaBranch::Flush(rsFloat scale = 1.0) {
  Clean();
  Init();
  pre_scale = scale;
  //  upsample_scale = 0;
}
 
//
// MultirateGenerator Implementation
// 

/// Constructor
MultirateGenerator::MultirateGenerator(rsFloat alpha, unsigned int branches)
{
  rsFloat beta = -(alpha - 2)/2.0;
  //Calculate the integer and fractional parts of beta
  int fint = (int)std::floor(beta);
  rsFloat ffrac = fmod(beta, 1); 
  //Build the multirate filter tree
  CreateTree(ffrac, fint, branches);
  scale = 1.0/std::pow(10.0, (-alpha+2)*2.0);
}

/// Destructor
MultirateGenerator::~MultirateGenerator() 
{
  delete topbranch;
}

/// Get a single noise sample
rsFloat MultirateGenerator::GetSample()
{
  return topbranch->GetSample()*scale;
}

/// Skip a number of samples, preseverving correlations
void MultirateGenerator::SkipSamples(long long samples)
{
  std::vector<FAlphaBranch*> flushbranches;
  int skip_branches = (int)std::floor(std::log10(samples))-1;
  //  rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Skip branches %d\n", skip_branches);
  if (skip_branches > 0) {
    FAlphaBranch *branch = topbranch;
    for (int i = 0; (i < skip_branches) && (branch != 0); i++) {
      flushbranches.push_back(branch);
      branch = branch->pre;
    }
    if (branch) {
      //Now generate the samples of the lower branches
      samples = (int)(samples/std::pow(10.0, (double)skip_branches));
      //rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Skipping %ld samples in %ld branches \n", samples, skip_branches);
      for (int i = 0; i < samples; i++)
        branch->GetSample();      
    }
    // Flush the buffers of the upper branches
    int size = flushbranches.size();
    flushbranches[size-1]->Flush(std::pow(10.0, skip_branches-2.0));
    for (int i = size-2; i >= 0; i--)
      flushbranches[i]->Flush();
  }
  else {
    for (int i = 0; i < samples; i++)
      topbranch->GetSample();
  }

}  

/// Create the branches of the filter structure tree
// The tree is stored as a linked list, which each link a branch
void MultirateGenerator::CreateTree(rsFloat ffrac, int fint, unsigned int branches) 
{
  if (branches == 0)
    throw std::runtime_error("Cannot create multirate noise generator with zero branches");
  // If ffrac and fint are both zero, we only need a single branch
  if ((ffrac == 0) && (fint == 0)) {
    topbranch = new FAlphaBranch(0, 0, 0, true);
  }
  else {
    topbranch = 0;
    for (unsigned int i = 0; i < branches-1; i++)
      topbranch = new FAlphaBranch(ffrac, fint, topbranch, false);
    topbranch = new FAlphaBranch(ffrac, fint, topbranch, true);
  }
}

/// Reset the output to zero
void MultirateGenerator::Reset()
{
  std::vector<FAlphaBranch*> flushbranches;
  //Build a reverse order list of branches
  FAlphaBranch *branch = topbranch;
  while (branch) {
    flushbranches.push_back(branch);
    branch = branch->pre;
  }
  // Flush the branch buffers in reverse order
  int size = flushbranches.size();
  for (int i = size-1; i >= 0; i--)
    flushbranches[i]->Flush();
}

//
// ClockModelGenerator Implementation
//

/// Constructor
ClockModelGenerator::ClockModelGenerator(const std::vector<rsFloat> &alpha, const std::vector<rsFloat> &in_weights, rsFloat frequency, rsFloat phase_offset, rsFloat freq_offset, int branches):
  phase_offset(phase_offset),
  freq_offset(freq_offset),
  frequency(frequency)
{
  weights = in_weights;
  std::vector<rsFloat>::const_iterator iter = alpha.begin();
  std::vector<rsFloat>::iterator witer = weights.begin();
  // Create the generators for each band
  for (; iter != alpha.end(); iter++, witer++) {
      MultirateGenerator *mgen = new MultirateGenerator(*iter, branches);
      generators.push_back(mgen);
      //Calibrate the weights using the measured calibration numbers
      if (*iter == 2) {
	*witer *= std::pow(10.0, 1.2250); 
      }
      else if (*iter == 1) {
	*witer *= std::pow(10.0, 0.25);
      }
      else if (*iter == 0) {
	*witer *= std::pow(10.0, -0.25);
      }
      else if (*iter == -1) {
	*witer *= std::pow(10.0, -0.5);
      }
      else if (*iter == -2) {
	*witer *= std::pow(10.0, -1);
      }
    }
  count = 0;
}

/// Destructor
ClockModelGenerator::~ClockModelGenerator()
{
  std::vector<MultirateGenerator *>::iterator iter;
  for (iter = generators.begin(); iter != generators.end(); iter++)
    delete *iter;
}

/// Get a single noise sample
rsFloat ClockModelGenerator::GetSample()
{
  rsFloat sample = 0;
  // Get noise from the multirate generators for each band
  int size = generators.size();
  for (int i = 0; i < size; i++) {
    sample += generators[i]->GetSample()*weights[i];
  }
  // Add the phase and frequency offsets
  sample += phase_offset;
  //Calculate the count in clock frequencies
  sample += 2*M_PI*freq_offset*count/rsParameters::rate();
  count++;
  return sample;  
}

/// Skip some noise samples, calculating only the branches required to preserve correlations
void ClockModelGenerator::SkipSamples(long long samples)
{
  int gens = generators.size();
  for (int i = 0; i < gens; i++)
    generators[i]->SkipSamples(samples);
  count += samples;
}

/// Reset the noise to zero
void ClockModelGenerator::Reset()
{
  int gens = generators.size();
  for (int i = 0; i < gens; i++)
    generators[i]->Reset();
  count = 0;
}

/// Is the generator going to produce non-zero samples
bool ClockModelGenerator::Enabled() {
  if ((!generators.empty()) || (freq_offset != 0) || (phase_offset != 0))
    return true;
  else
    return false;
}

//
//PythonNoiseGenerator Implementation
//

///Constructor
PythonNoiseGenerator::PythonNoiseGenerator(const std::string& module, const std::string& function):
  generator(module, function)
{

}

///Destructor
PythonNoiseGenerator::~PythonNoiseGenerator()
{
}

///Get a single noise sample
rsFloat PythonNoiseGenerator::GetSample()
{
  return generator.GetSample();
}

