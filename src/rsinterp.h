//rsinterp.cpp - Classes for interpolation of sets of data
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//11 June 2007

#ifndef __RS_INTERP_H
#define __RS_INTERP_H

#include <vector>
#include <config.h>

namespace rs {

  class InterpSetData;

  /// This class interpolates y, given x based on a number of samples of y
  /// The samples do not need to be equally spaced in x
  class InterpSet {
    //InterpSet is just a thin wrapper around InterpSetData, the class which implements interpolation
    //this was done to allow interpolation to be performed with different algorithms, as needed
  public:
    /// Constructor
    InterpSet();
    /// Destructor
    ~InterpSet();
    /// Load a number of samples into the set
    void LoadSamples(const std::vector<rsFloat> &x, const std::vector<rsFloat> &y);
    /// Load a single sample into the set
    void InsertSample(rsFloat x, rsFloat y);
    /// Get the interpolated value at the given point
    rsFloat Value(rsFloat x);
    /// Get the maximum value in the set
    rsFloat Max() const;
    /// Divide every sample in the set by a given number
    void Divide(rsFloat a);
  private:
    InterpSetData *data; //<! The data is hidden behind a pointer to the implementation
  };


}

#endif
