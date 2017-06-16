//Test the noise generation code from fers
#include <config.h>
#include <iostream>
#include <boost/thread.hpp>
#include "rsnoise.h"
#include "rsdsp.h"

using namespace rs;
using namespace std;

unsigned int processors = 4;

//test noise generation in threaded code

/*void NoiseTest() {
  std::vector<rsFloat> alphas;
  std::vector<rsFloat> weights;
  //Initialize for standard five parameter clock model
  for (int i = 0; i < 5; i++) {
    alphas.push_back(i);
    weights.push_back(1);
  }
  ClockModelGenerator *gen = new ClockModelGenerator(alphas, weights, 0, 0, (int)1e6, true);
  for (int i = 0; i < 1e6; i++) {
    rsFloat j = gen->GetSample();
  }
  delete gen;
  }*/

int main()
{
  rsNoise::InitializeNoise();
  /*
   *rsFloat a[12] = {1.000000000000000e+00,
			   -1.009905033992710e+01,
			   4.667711400358178e+01,
			   -1.303037481274771e+02,
			   2.440751993455667e+02,
			   -3.220568002997372e+02,
			   3.054264698325060e+02,
			   -2.081636321246607e+02,
			   9.991255272037478e+01,
			   -3.216159801323370e+01,
			   6.248589053855025e+00,
			   -5.550959230019682e-01};
  rsFloat b[12] = {   1.354744024510072e-04,
			       -1.034084955129721e-03,
			       3.530891199610175e-03,
			       -6.866941761673958e-03,
			       7.725944948126016e-03,
			       -3.491219910193682e-03,
			       -3.491219910193594e-03,
			       7.725944948125953e-03,
			       -6.866941761673926e-03,
			       3.530891199610165e-03,
			       -1.034084955129719e-03,
			       1.354744024510070e-04};
                   */
  //  IIRFilter* iir = new IIRFilter(a, b, 12);
  //for (int i = 0; i < 10240; i++)
  // cout << iir->Filter(rsNoise::WGNSample(1)) << endl;

  DecadeUpsampler upsample;
  for (int i = 0; i < 1e4; i++) {
    rsFloat buf[10];
    upsample.Upsample(rsNoise::WGNSample(1), buf);
    for (int j = 0; j < 10; j++)
      cout << buf[j] << endl;
  }
}
