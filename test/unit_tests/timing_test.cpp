//Test the clock model code from fers
#include <config.h>
#include <iostream>
#include <boost/thread.hpp>
#include <sstream>
#include "rsnoise.h"
#include "rstiming.h"
#include "rsparameters.h"

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

int main(int argc, char *argv[])
{
  //Get command line args
  if (argc != 4) {
    cout << "Usage: timing_test alpha weight count" << endl;
  }
  else {
    istringstream iss1(argv[1]);
    double alpha;
    iss1 >> alpha;
    double weight;
    istringstream iss2(argv[2]);
    iss2 >> weight;
    double count;
    istringstream iss3(argv[3]);
    iss3 >> count;
    cerr << "Alpha: " << alpha << " weight: " << weight << " count: " << count << endl;
    //Initialize noise generation
    rsNoise::InitializeNoise();
    rsParameters::modify_parms()->SetRate(1e3);
    PrototypeTiming proto("test");
    proto.AddAlpha(alpha, weight);
    proto.SetFrequency(1e6);
    ClockModelTiming timing("test");
    timing.InitializeModel(&proto);
      for (int i = 0; i < count; i++)
	cout << timing.NextNoiseSample() << endl;
  }
}
