//rspulserender.h
//Definitions for pulse rendering functions
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//7 June 2006

#ifndef __RS_PULSE_RENDER
#define __RS_PULSE_RENDER

#include <config.h>
#include <vector>
#include <queue>
#include <string>

//Forward definition of boost threads classes (see boost threads)
namespace boost {
  class mutex;
}

namespace rs {
  //Forward declaration of Response (see rsresponse.h)
  class Response;
  //Forward declaration of Receiver (see rsradar.h)
  class Receiver;




  /// Export the responses received by a receiver to an XML file
  void ExportReceiverXML(const std::vector<rs::Response*> &responses, const std::string filename);

  /// Export the receiver pulses to the specified binary file, using the specified quantization
  void ExportReceiverBinary(const std::vector<rs::Response*> &responses, rs::Receiver* recv, const std::string recv_name, const std::string filename);

  /// Export the receiver responses to the specified CSV value files
  void ExportReceiverCSV(const std::vector<rs::Response*> &responses, const std::string filename);

  /// Management class for threaded rendering
  class ThreadedRenderer {
  public:
    /// Constructor
    ThreadedRenderer(const std::vector<rs::Response*> *responses, const rs::Receiver* recv, int max_threads);
    /// Destructor
    ~ThreadedRenderer();
    /// Render all the responses in a single window
    void RenderWindow(rsComplex *window, rsFloat length, rsFloat start, rsFloat frac_delay);
  private:
    const std::vector<rs::Response*> *responses; //!< Vector of target responses seen by this receiver
    const rs::Receiver* recv; //!< Receiver we are rendering for
    int max_threads; //!< The maximum allowed thread count for rendering
  };

  /// Single thread for rendering
  class RenderThread {
  public:
    /// Constructor
    RenderThread(int serial, boost::mutex *window_mutex, rsComplex *window, rsFloat length, rsFloat start, rsFloat frac_delay, boost::mutex *work_list_mutex, std::queue<Response *> *work_list);
    /// Destructor
    ~RenderThread();
    /// Step through the worklist, rendering the required responses
    void operator()();
  private:
    /// Get a response from the worklist, returning NULL on failure
    Response *GetWork();
    /// Add the array to the window, locking the window lock in advance
    void AddWindow(rsComplex *array, rsFloat start_time, unsigned int array_size);
    int serial; //!< Serial number of this thread
    boost::mutex *window_mutex; //!< Mutex to protect window
    rsComplex *window; //!< Pointer to render window
    rsFloat length; //!< Length of render window (seconds)
    rsFloat start; //!< Start time of render window (seconds)
    rsFloat frac_delay; //!< Fractional window start time (< 1 sample, samples)
    boost::mutex *work_list_mutex; //!< Mutex to protect work list
    std::queue<Response*> *work_list; //!< List of responses to render
    rsComplex *local_window;
    
  };
    
}


#endif
