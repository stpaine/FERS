// Test clock modeling classes
#include <config.h>
#include <iostream>
#include "fftwcpp.h"
#include "rstiming.h"
#include "rsnoise.h"

using namespace rs;

unsigned int processors = 4;

//Create the prototype timing class
PrototypeTiming* MakeProtoTiming()
{
  //Two parameter clock model
  PrototypeTiming *proto = new PrototypeTiming("test");
  proto->AddAlpha(0, 0.05); //White PM
  proto->AddAlpha(2, 0.95); //White FM
  proto->SetFrequency(1e9);
  //proto->AddAlpha(3, 0); //Flicker FM
  return proto;
}

// Test the pulsetiming class
int TestClockModelTiming(const PrototypeTiming *proto)
{
  int pulses = 3;
  int pulse_length = 1000;
  ClockModelTiming *timing = new ClockModelTiming("test_pulse");
  timing->InitializeModel(proto);
  
  for (int i = 0; i < pulses; i++) {
    rsFloat j = timing->GetPulseTimeError();
    //    std::cout << j << std::endl;
  }
  for (int k = 0; k < pulses; k++) {
    for (int i = 0; i < pulse_length; i++) {
      rsFloat j = timing->NextNoiseSample();
      //      std::cout << j << std::endl;
    }
  }
  delete timing;
  return 0;
}


int main()
{
  rsNoise::InitializeNoise();
  PrototypeTiming *timing = MakeProtoTiming();

  //TestPulseTiming(timing);
    TestClockModelTiming(timing);

  rsNoise::CleanUpNoise();
  delete timing;
}
