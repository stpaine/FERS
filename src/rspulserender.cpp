//rspulserender.cpp
//Performs the second phase of the simulation - rendering the result
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//7 June 2006

#include <stdexcept>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <cstdio> //The C stdio functions are used for binary export
#include <boost/scoped_array.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/xtime.hpp>

#include "rsradar.h"
#include "rspulserender.h"
#include "rsresponse.h"
#include "rsdebug.h"
#include "rsparameters.h"
#include "rsnoise.h"
#include "rshdf5.h"
#include "rssignal.h"
#include "rstiming.h"
#include "rsdsp.h"
#include "rsportable.h"

#define TIXML_USE_STL
#include <tinyxml.h>

using namespace rs;

namespace {

  ///Open the binary file for response export
  int OpenHDF5File(const std::string &recv_name)
  {
    int hdf5_file = 0;
    if (rs::rsParameters::export_binary()) {
      //Build the filename for the binary file
      std::ostringstream b_oss;
      b_oss.setf(std::ios::scientific);    
      b_oss << recv_name << ".h5";
      //Open the binary file for writing
      hdf5_file = rshdf5::CreateFile(b_oss.str().c_str());
    }
    return hdf5_file;    
  }

  /// Add noise to the window, to simulate receiver noise
  void AddNoiseToWindow(rs::rsComplex *data, unsigned int size, rsFloat temperature)
  {
    if (temperature == 0)
      return;  
    //Calculate the noise power
    rsFloat power = rsNoise::NoiseTemperatureToPower(temperature, rsParameters::rate()*rsParameters::oversample_ratio()/2);
    WGNGenerator generator(sqrt(power)/2.0);
    //Add the noise to the signal
    for (unsigned int i = 0; i < size; i++)
      data[i] += rsComplex(generator.GetSample(), generator.GetSample());
  }

  /// Normalize a window, and quantize the signal as required
  rsFloat QuantizeWindow(rs::rsComplex *data, unsigned int size)
  {
    rsFloat max = 0;
    //Get the full scale amplitude of the pulse
    for (unsigned int i = 0; i < size; i++) {
      if (std::fabs(data[i].real()) > max)
        max = std::fabs(data[i].real());
      if (std::fabs(data[i].imag()) > max)
        max = std::fabs(data[i].imag());
        if (isnan(data[i].real()) || (isnan(data[i].imag())))
          throw std::runtime_error("NAN in QuantizeWindow -- early");
    }
    if (rsParameters::adc_bits() > 0) {
      //Simulate quantization and normalize
      rsSignal::ADCSimulate(data, size, rsParameters::adc_bits(), max);
    }
    else {
      //Just normalize
      if (max != 0) {
        for (unsigned int i = 0; i < size; i++) {
          data[i] /= max;
          if (isnan(data[i].real()) || (isnan(data[i].imag())))
            throw std::runtime_error("NAN in QuantizeWindow -- late");
        }
      }
    }
    //Return the full scale value
    return max;
  }

  /// Add a response array to the receive window
  void AddArrayToWindow(rsFloat wstart, rs::rsComplex *window, unsigned int wsize, rsFloat rate, rsFloat rstart, rs::rsComplex *resp, unsigned int rsize)
  {
    int start_sample = static_cast<int>(rsPortable::rsRound(rate*(rstart-wstart)));
    //Get the offset into the response
    unsigned int roffset = 0;
    if (start_sample < 0) {
      roffset = -start_sample;
      start_sample = 0;
    }
    //Loop through the response, adding it to the window
    for (unsigned int i = roffset; (i < rsize) && ((i+start_sample) < wsize); i++)
      window[i+start_sample] += resp[i];      
  }

  /// Generate the phase noise samples for the receive window
  rsFloat* GeneratePhaseNoise(const rs::Receiver *recv, unsigned int wsize, rsFloat rate, rsFloat &carrier, bool &enabled) 
  {
    //Get a pointer to the receiver's timing object
    ClockModelTiming *timing = dynamic_cast<ClockModelTiming *>(recv->GetTiming()); 
    if (!timing)
      throw std::runtime_error("[BUG] Could not cast receiver->GetTiming() to ClockModelTiming");
    //Allocate memory for the phase noise samples
    rsFloat *noise = new rsFloat[wsize];
    enabled = timing->Enabled();
    if (enabled) {
      //Generate phase noise if timing is enabled
      for (unsigned int i = 0; i < wsize; i++) {
	noise[i] = timing->NextNoiseSample();        
      }
      std::printf("%g\n", noise[0]);
      //Calculate the number of samples to skip
      if (timing->GetSyncOnPulse()) {
	timing->Reset();
	int skip = (int)std::floor(rate*recv->GetWindowSkip());
	timing->SkipSamples(skip);
      }
      else {
        long skip = std::floor(rate/recv->GetPRF()-rate*recv->GetWindowLength());
        //for (long i = 0; i < skip; i++)
        //std::printf("%g\n", timing->NextNoiseSample());
        timing->SkipSamples(skip);
      }      
      carrier = timing->GetFrequency();
    }
    else {
      // Clear samples if not
      for (unsigned int i = 0; i < wsize; i++)
	noise[i] = 0;
      carrier = 1;
    }
    return noise;
  }

  /// Add clock phase noise to the receive window
  void AddPhaseNoiseToWindow(const rsFloat* noise, rs::rsComplex *window, unsigned int wsize)
  {    
      for (unsigned int i = 0; i < wsize; i++) {
	//Modify the phase
        if (isnan(noise[i]))
          throw std::runtime_error("[BUG] Noise is NAN in AddPhaseNoiseToWindow");	
	//	std::printf("%g\n", noise[i]);
	rsComplex mn = exp(rsComplex(0.0, 1.0)*noise[i]);
	window[i] *= mn;
        if (isnan(window[i].real()) || (isnan(window[i].imag())))
          throw std::runtime_error("[BUG] NAN encountered in AddPhaseNoiseToWindow");
      }
  }

  /// Export to the FersBin file format
  void ExportResponseFersBin(const std::vector<rs::Response*>& responses, const rs::Receiver* recv, const std::string& recv_name)
  {
    //Bail if there are no responses to export
    if (responses.empty())
      return;

    int out_bin = OpenHDF5File(recv_name);

    // Create a threaded render object, to manage the rendering process
    ThreadedRenderer thr_renderer(&responses, recv, rsParameters::render_threads());

    // Now loop through the responses and write them to the file
    int window_count = recv->GetWindowCount();
    for (int i = 0; i < window_count; i++) {
      rsFloat length = recv->GetWindowLength();
      rsFloat rate = rsParameters::rate()*rsParameters::oversample_ratio();
      unsigned int size = static_cast<unsigned int>(std::ceil(length*rate));
      //      rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Length: %g Size: %d\n", length, size);
      //Generate the phase noise samples for the window
      rsFloat carrier;
      bool pn_enabled;
      rsFloat* pnoise = GeneratePhaseNoise(recv, size, rate, carrier, pn_enabled);
      //      printf("%g\n", pnoise[0]);
      // Get the window start time, including clock drift effects
      rsFloat start = recv->GetWindowStart(i) + (pnoise[0] / (2 * M_PI * carrier));
      // Get the fraction of a sample of the window delay
      rsFloat frac_delay = start * rate - rsPortable::rsRound(start * rate);
      //      std::printf("frac: %g %g %g\n", frac_delay, start*rate, rsPortable::rsRound(start*rate));
      start = rsPortable::rsRound(start*rate)/rate;
      //rsFloat start = recv->GetWindowStart(i);
      // Allocate memory for the entire window
      rsComplex* window = new rsComplex[size];
      //Clear the window in memory
      memset(window, 0, sizeof(rsComplex)*size);
      //Add Noise to the window
      AddNoiseToWindow(window, size, recv->GetNoiseTemperature());
      // Render to the window, using the threaded renderer
      //std::printf("%g\n", start);
      thr_renderer.RenderWindow(window, length, start, frac_delay);
      //Downsample the contents of the window, if appropriate
      if (rsParameters::oversample_ratio() != 1) {
	// Calculate the size of the window after downsampling
	unsigned int new_size = size / rsParameters::oversample_ratio();
	// Allocate memory for downsampled window
	rsComplex* tmp = new rsComplex[new_size];
	//Downsample the data into tmp
	Downsample(window, size, tmp, rsParameters::oversample_ratio());
	// Set tmp as the new window
	size = new_size;
	delete[] window;
	window = tmp;
      }
      //Add Phase noise to the window
      if (pn_enabled)
        AddPhaseNoiseToWindow(pnoise, window, size);
      //Clean up the phase noise array
      delete[] pnoise;
      //Normalize and quantize the window
      rsFloat fullscale = QuantizeWindow(window, size);
      //Export the binary format
      if (rs::rsParameters::export_binary()) {
	rshdf5::AddChunkToFile(out_bin, window, size, start, rsParameters::rate(), fullscale, i);
      }
      // Clean up memory
      delete[] window;
    } // for (i = 1:windows)
    // Close the binary and csv files
    if (out_bin) rshdf5::CloseFile(out_bin);
  }

}

/// Export the responses received by a receiver to an XML file
void rs::ExportReceiverXML(const std::vector<rs::Response*> &responses, const std::string filename)
{
  //Create the document
  TiXmlDocument doc;
  TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "", "");
  doc.LinkEndChild(decl);
  //Create a root node for the document
  TiXmlElement *root = new TiXmlElement("receiver");
  doc.LinkEndChild(root);

  //dump each response in turn
  std::vector<rs::Response*>::const_iterator ri;
  for (ri = responses.begin(); ri != responses.end(); ri++)
    (*ri)->RenderXML(root);

  //write the output to the specified file
  doc.SaveFile(filename+".fersxml");  
}

/// Export the responses in CSV format
void rs::ExportReceiverCSV(const std::vector<rs::Response*> &responses, const std::string filename)
{
  std::map<std::string, std::ofstream *> streams; //map of per-transmitter open files
  std::vector<rs::Response*>::const_iterator iter;
  for (iter = responses.begin(); iter != responses.end(); iter++)
    {
      std::ofstream *of;
      // See if a file is already open for that transmitter
      std::map<std::string, std::ofstream *>::iterator ofi = streams.find((*iter)->GetTransmitterName());
      // If the file for that transmitter does not exist, add it
      if (ofi == streams.end())	{
        std::ostringstream oss;  
        oss << filename << "_" << (*iter)->GetTransmitterName() << ".csv";
	  //Open a new ofstream with that name
	  of = new std::ofstream(oss.str().c_str());
	  of->setf(std::ios::scientific); //Set the stream in scientific notation mode
	  //Check if the open succeeded
	  if (!(*of))
	    throw std::runtime_error("[ERROR] Could not open file "+oss.str()+" for writing");
	  //Add the file to the map
	  streams[(*iter)->GetTransmitterName()] = of;
	}
      else { //The file is already open
	of = (*ofi).second;
      }
      // Render the response to the file
      (*iter)->RenderCSV(*of);
    }
  //Close all the files that we opened
  std::map<std::string, std::ofstream *>::iterator ofi;
  for (ofi = streams.begin(); ofi != streams.end(); ofi++)
    delete (*ofi).second;
}

/// Export the receiver pulses to the specified binary file, using the specified quantization
void rs::ExportReceiverBinary(const std::vector<rs::Response *> &responses, Receiver* recv, const std::string recv_name, const std::string filename)
{
  ExportResponseFersBin(responses, recv, recv_name);
}

//
// ThreadedRenderer Implementation
//

/// Constructor
ThreadedRenderer::ThreadedRenderer(const std::vector<rs::Response*> *responses, const rs::Receiver* recv, int max_threads):
  responses(responses),
  recv(recv),
  max_threads(max_threads)
{
}
 
/// Destructor
ThreadedRenderer::~ThreadedRenderer()
{
}

/// Render all the responses in a single window
void ThreadedRenderer::RenderWindow(rsComplex *window, rsFloat length, rsFloat start, rsFloat frac_delay)
{
  rsFloat end = start+length; // End time of the window
  //Put together a list of responses seen by this window
  std::queue<Response *> work_list;
  std::vector<Response *>::const_iterator iter = responses->begin();
  for (; iter != responses->end(); iter++) {
    rsFloat resp_start = (*iter)->StartTime();
    rsFloat resp_end = (*iter)->EndTime();
    if ((resp_start <= end) && (resp_end >= start))
      work_list.push(*iter);
  }
  //Manage threads with boost::thread_group
  boost::thread_group group;
  // Create a mutex to protect the worklist
  boost::mutex work_list_mutex;
  //Create a mutex to protect the window
  boost::mutex window_mutex;
  //Create the number of threads we are allowed
  std::vector<RenderThread *> threads;
  {
    boost::mutex::scoped_lock lock(work_list_mutex);
    //   rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Spawning %d render threads\n", max_threads);
    for (int i = 0; i < max_threads; i++) {
      //      rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Spawning %d\n", i);
      RenderThread *thr = new RenderThread(i, &window_mutex, window, length, start, frac_delay, &work_list_mutex, &work_list);
      threads.push_back(thr);      
      group.create_thread(*thr);
    }
  }
  // Wait until all our threads are complete
  group.join_all();  
  //Clean up the thread objects
  std::vector<RenderThread *>::iterator thread_iter = threads.begin();
  for (; thread_iter < threads.end(); thread_iter++)
    delete *thread_iter;
}

//
// RenderThread Implementation
//

/// Constructor
RenderThread::RenderThread(int serial, boost::mutex *window_mutex, rsComplex *window, rsFloat length, rsFloat start, rsFloat frac_delay, boost::mutex *work_list_mutex, std::queue<Response *> *work_list):
  serial(serial),
  window_mutex(window_mutex),
  window(window),
  length(length),
  start(start),
  frac_delay(frac_delay),
  work_list_mutex(work_list_mutex),
  work_list(work_list)
{
}

/// Destructor
RenderThread::~RenderThread()
{
}

/// Step through the worklist, rendering the required responses
void RenderThread::operator()()
{
  rsFloat rate = rsParameters::rate()*rsParameters::oversample_ratio();
  unsigned int size = (unsigned int)(std::ceil(length*rate));
  // Allocate memory for the local window
  local_window = new rsComplex[size];
  for (unsigned int i = 0; i < size; i++)
    local_window[i] = 0;
  // Loop until the work list is empty
  Response *resp = GetWork();
  while (resp) {      
    unsigned int psize;
    rsFloat prate;
    // Render the pulse into memory
    boost::shared_array<rs::rsComplex> array = resp->RenderBinary(prate, psize, frac_delay);
    // Add the array to the window
    AddWindow(array.get(), resp->StartTime(), psize);
    // Get more work, if it's available
    resp = GetWork();
  }
  {
    // Lock the window mutex, to ensure that we are the only thread accessing the window
    boost::mutex::scoped_lock lock(*window_mutex);
    // Add local window to the global window
    for (unsigned int i = 0; i < size; i++)
      window[i] += local_window[i];
  }
  delete[] local_window;
}

/// Add the array to the window, locking the window lock in advance
void RenderThread::AddWindow(rsComplex *array, rsFloat start_time, unsigned int array_size) {
  //Calculate required window parameters
  rsFloat rate = rsParameters::rate()*rsParameters::oversample_ratio();
  unsigned int size = (unsigned int)(std::ceil(length*rate));
  // Add the array to the correct place in the local window
  AddArrayToWindow(start, local_window, size, rate, start_time, array, array_size);  
}

/// Get a response from the worklist, returning NULL on failure
Response *RenderThread::GetWork() {
  // rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Thread %d getting work\n", serial);
  Response *ret;
  { //Scope for scoped lock
    //Lock the work list mutex
    boost::mutex::scoped_lock lock(*work_list_mutex);    
    if (work_list->empty()) //If the work list is empty, return NULL
      return 0;
    // Get the response at the front of the queue
    ret = work_list->front();
    work_list->pop();
  }
  return ret; 
  //scoped_lock unlocks mutex here
}
