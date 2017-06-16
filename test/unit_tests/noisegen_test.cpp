//Test the noise generation code from fers
#include <config.h>
#include <iostream>
#include <boost/thread.hpp>
#include "rsnoise.h"

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
  //Initialize noise generation
  rsNoise::InitializeNoise();
  //Create the generator (for alpha = 2)
  MultirateGenerator *gen = new MultirateGenerator(1, 4);
  for (int i = 0; i < 1e7; i++)
    /*cout <<*/ gen->GetSample(); /*<< endl;*/
  delete gen;
  // Spawn threads to test thread safety of ClockModelGenerator
  //  boost::thread_group man;
  //for (unsigned int i = 0; i < processors; i++) {
  // man.create_thread(NoiseTest);
  //}
  // man.join_all();

}
