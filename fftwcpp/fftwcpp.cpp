//Implementation of fft_manager.h and fft_plan.h
//See that file for information on what it does and how
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za

//07/04/2006: Initial Revision

#include <stdexcept>
#include <map>
#include <fftw3.h>
#include "fftwcpp.h"
#include <boost/thread/mutex.hpp> //Included for boost::mutex
#include "fftwcpp_templates.h"


//Mutexes for accessing the planner, local to this module
boost::mutex plannerMutex;

FFTManager* FFTManager::instance = 0;

///Aligned memory deallocation
void FFTAlignedFree(void *ptr)
{
  FFTManager::AlignedFree(ptr);
}

///Initialize FFTW
void FFTInit(int threads)
{
#ifdef FFTW_WITH_THREADS
  fftw_init_threads();
  fftw_plan_with_nthreads(threads);
#endif
}

///Clean up the FFTW code, freeing memory
void FFTCleanUp()
{
  FFTManager::Instance()->Clean();
  //Clean up fftw
#ifdef FFTW_WITH_THREADS
  fftw_cleanup_threads();
#else
  fftw_cleanup();
#endif
  delete FFTManager::Instance();  
}



//Constructor for FFTManager
//This constructor is private and is only called once
FFTManager::FFTManager() {
}

//Function which returns a pointer to the only instance of FFTManager
FFTManager *FFTManager::Instance() {
  if (instance == 0)
    instance = new FFTManager;
  return instance;
}

void * FFTManager::AlignedMalloc(int size) {
  return fftw_malloc(size);
}

void FFTManager::AlignedFree(void *ptr) {
  fftw_free(ptr);
}

//Get a reference to a complex plan, optionally creating it if it does not exist 
FFTComplex * FFTManager::GetComplexPlan(int size, bool create, Complex *in, Complex *out) {
  FFTComplex *plan;
  plan = complex_plans[size];
  //If the plan does not exist, we can optionally create it
  if ((plan == 0) && (create)) {
    plan = new FFTComplex(size, in, out, FFTComplex::FFT_FORWARD);
    complex_plans[size] = plan;
  }
  //Return the plan, or NULL if there is no plan
  return plan;
}

FFTComplex *FFTManager::GetComplexPlanInv(int size, bool create, Complex *in, Complex *out) {
  FFTComplex *plan;
  plan = complex_inv_plans[size];
  //If the plan does not exist, we can optionally create it
  if ((plan == 0) && (create)) {
    plan = new FFTComplex(size, in, out, FFTComplex::FFT_BACKWARD);
    complex_inv_plans[size] = plan;
  }
  //Return the plan
  if (plan)
    return plan;
  else
    throw FFTException("Could not create complex plan");
}

/// Get a reference to a real->complex plan, and create on if it doesn't exist
FFTRealForward* FFTManager::GetRealForwardPlan(int size, rsFloat* in, Complex* out) {
  FFTRealForward *plan;
  plan = real_forward_plans[size];
  // If the plan does not exist, create it
  if (plan == 0) {
    plan = new FFTRealForward(size, in, out);
    real_forward_plans[size] = plan;
  }
  //Return the plan
  return plan;
}

/// Get a reference to a complex->real plan, and create on if it doesn't exist
FFTRealInverse* FFTManager::GetRealInversePlan(int size, Complex* in, rsFloat* out) {
  FFTRealInverse *plan;
  plan = real_inverse_plans[size];
  // If the plan does not exist, create it
  if (plan == 0) {
    plan = new FFTRealInverse(size, in, out);
    real_inverse_plans[size] = plan;
  }
  //Return the plan
  return plan;
}

//Clean up the manager and destroy all plans
void FFTManager::Clean() {
  std::map <int, FFTComplex *>::iterator iter;
  for (iter = complex_plans.begin(); iter != complex_plans.end(); iter++)
    delete (*iter).second;
  complex_plans.clear();
  for (iter = complex_inv_plans.begin(); iter != complex_inv_plans.end(); iter++)
    delete (*iter).second;
  std::map <int, FFTRealInverse *>::iterator ri_iter;
  for (ri_iter = real_inverse_plans.begin(); ri_iter != real_inverse_plans.end(); ri_iter++)
    delete (*ri_iter).second;
  std::map <int, FFTRealForward *>::iterator rf_iter;
  for (rf_iter = real_forward_plans.begin(); rf_iter != real_forward_plans.end(); rf_iter++)
    delete (*rf_iter).second;
  complex_inv_plans.clear();
}

//
// Implementation of FFTComplex
//

//Perform an FFT on two arrays of data
void FFTComplex::transform(int size, Complex *in, Complex *out)
{
  if (plan == 0)
    throw FFTException("Can not perform transform on NULL plan.");
  fftw_execute_dft((fftw_plan_s *)plan, reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out));
}

//Constructor which creates either a forward or reverse transform
FFTComplex::FFTComplex(int size, Complex *in, Complex *out, FFTComplex::fft_direction direction)
{
  if (in == out)
    throw FFTException("[BUG] In place transforms not supported");
  boost::mutex::scoped_lock lock(plannerMutex); //Lock the mutex
  if (direction == FFT_FORWARD)
    plan = fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out), FFTW_FORWARD, FFTW_ESTIMATE);
  else
    plan = fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out), FFTW_BACKWARD, FFTW_ESTIMATE);
  //Scoped lock automatically unlocks the mutex at the end of the scope
}

FFTComplex::~FFTComplex() {
  boost::mutex::scoped_lock lock(plannerMutex); //Lock the planner mutex
  fftw_destroy_plan(reinterpret_cast<fftw_plan>(plan));
  //Scoped lock unlocks mutex here
}

//
// Implementation of FFTRealForward
// 

void FFTRealForward::transform(int size, rsFloat *in, Complex *out)
{
  if (plan == 0)
    throw FFTException("[BUG] Can not transform on NULL plan.");
  fftw_execute_dft_r2c((fftw_plan_s *)plan, in, reinterpret_cast<fftw_complex *>(out));
}

// Constructor which creates a forward transform
FFTRealForward::FFTRealForward(int size, rsFloat *in, Complex *out)
{
  //Lock the planner Mutex
  boost::mutex::scoped_lock lock(plannerMutex);
  plan = fftw_plan_dft_r2c_1d(size, in, reinterpret_cast<fftw_complex *>(out), FFTW_ESTIMATE);  
  //scoped_lock will unlock planner mutex here
}

// Destructor
FFTRealForward::~FFTRealForward()
{
  boost::mutex::scoped_lock lock(plannerMutex); //Lock the planner mutex
  fftw_destroy_plan(reinterpret_cast<fftw_plan>(plan));
  //Scoped lock unlocks mutex here
}

//
// Implementation of FFTRealInverse
// 

void FFTRealInverse::transform(int size, Complex *in, rsFloat *out)
{
  if (plan == 0)
    throw FFTException("[BUG] Can not transform on NULL plan.");
  fftw_execute_dft_c2r((fftw_plan_s *)plan, reinterpret_cast<fftw_complex *>(in), out);
}

// Constructor which creates an inverse transform
FFTRealInverse::FFTRealInverse(int size, Complex *in, rsFloat *out)
{
  //Lock the planner Mutex
  boost::mutex::scoped_lock lock(plannerMutex);
  plan = fftw_plan_dft_c2r_1d(size, reinterpret_cast<fftw_complex *>(in), out, FFTW_ESTIMATE);  
  //scoped_lock will unlock planner mutex here
}

// Destructor
FFTRealInverse::~FFTRealInverse()
{
  boost::mutex::scoped_lock lock(plannerMutex); //Lock the planner mutex
  fftw_destroy_plan(reinterpret_cast<fftw_plan>(plan));
  //Scoped lock unlocks mutex here
}

//
// Implementation of FFTExceptiom
//

/// Constructor
FFTException::FFTException(const std::string &str):
  std::runtime_error(str)
{
}
