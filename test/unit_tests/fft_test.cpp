#include <iostream>
#include "fftwcpp.h"
#include "cycle.h"
#include <ctime>
#include <cstdlib>

using namespace std;

typedef complex<double> cComplex;


void simpleTest() {
  complex<double> *in, *out;
  FFTManager *fftw = FFTManager::Instance();
  FFTComplex *fw, *rev;
  //Allocate memory
  in = (cComplex *)(fftw->AlignedMalloc(128*sizeof(complex<double>)));
  out = (cComplex *)(fftw->AlignedMalloc(128*sizeof(complex<double>)));
  //Create a plan
  fw = fftw->GetComplexPlan(128, true, in, out);
  rev = fftw->GetComplexPlanInv(128, true, in, out);
  //Fill the array
  memset(in, 0, sizeof(complex<double>)*128);
  in[64] = complex<double>(128, 0);
  fw->transform(128, in, out);
  cout << out << " " << in << endl;
  for (int i = 0; i < 20; i++)
    cout << out[i] << ",";
  cout << endl;
  
  fftw->AlignedFree(in);
  fftw->AlignedFree(out);
}

void lookupTest() {
  complex<double> *in, *out;
  const int N = 2;
  int i, sizes[N];
  FFTManager *fftw = FFTManager::Instance();
  ticks start, end;
  srand(time(NULL));
  for (i = 0; i < N; i++)
    sizes[i] = rand() % 1024;
  //Allocate memory
  in = (cComplex *)(fftw->AlignedMalloc(1024*sizeof(complex<double>)));
  out = (cComplex *)(fftw->AlignedMalloc(1024*sizeof(complex<double>)));
  //Create all the plans
  cout << "Creating " << N << " plans" << endl;
  start = getticks();
  for (i = 0; i < N; i++)
    fftw->GetComplexPlan(sizes[i], true, in, out);
  end = getticks();
  cout << "Took " << elapsed(end, start) << " ticks" << endl;
  //Look up the plans
  cout << "Looking up " << N << " plans" << endl;
  start = getticks();
  for (i = 0; i < N; i++)
    fftw->GetComplexPlan(sizes[i], true, in, out);
  end = getticks();
  cout << "Took " << elapsed(end, start) << " ticks" << endl;

  fftw->AlignedFree(in);
  fftw->AlignedFree(out);
}

void computeTest() {
  complex <double> *in, *out;
  FFTComplex *fw, *rev;
  const int size = 102400;
  FFTManager *fftw = FFTManager::Instance();
  ticks start, end;
  //Allocate memory
  in = (cComplex *)(fftw->AlignedMalloc(size*sizeof(complex<double>)));
  out = (cComplex *)(fftw->AlignedMalloc(size*sizeof(complex<double>)));

  //Make the plan
  cout << "Making plan for large transform" << endl;
  start = getticks();
  fw = fftw->GetComplexPlan(size, true, in, out);
  rev = fftw->GetComplexPlanInv(size, true, in, out);
  end = getticks();
  cout << "Took " << elapsed(end, start) << " ticks" << endl;

  //Clear memory
  memset(in, 0, sizeof(complex<double>)*size);  
  in[size/2] = cComplex(size, 0);

  //Run the transform
  cout << "Running Forward Transform of size " << size << endl;
  start = getticks();
  fw->transform(size, in, out);
  end = getticks();
  cout << "Took " << elapsed(end, start) << " ticks" << endl;

  cout << "Normalize of size " << size << endl;
  start = getticks();
  for (int i = 0; i < size; i++)
    out[i] /= size;
  end = getticks();
  cout << "Took " << elapsed(end, start) << " ticks" << endl;

  

  for (int i = 0; i < 10; i++)
    cout << out[i] << ",";
  cout << endl;

  cout << "Running Reverse Transform of size " << size << endl;
  start = getticks();
  rev->transform(size, out, in);
  end = getticks();
  cout << "Took " << elapsed(end, start) << " ticks" << endl;

  for (int i = size/2-5; i < size/2+5; i++)
    cout << in[i] << ",";
  cout << endl;

  fftw->AlignedFree(in);
  fftw->AlignedFree(out);

  fftw->Clean();
}

int main()
{
  simpleTest();
  lookupTest();
  computeTest();
}

  
  
  
