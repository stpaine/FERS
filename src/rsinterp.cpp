//rsinterp.cpp - Implements interpolation class
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//11 June 2007

#include <map>
#include <utility>
#include <stdexcept>
#include <cmath>
#include "rsinterp.h"
#include "rsdebug.h"

using std::map;
using std::vector;
using std::pair;

using namespace rs;

//
// InterpSetData Implementation
//

///Data storage class for the InterpSet class
class rs::InterpSetData {
public:
  ///Load samples into the set
  void LoadSamples(const vector<rsFloat> &x, const vector<rsFloat> &y);
  ///Load a single sample into the set
  void InsertSample(rsFloat x, rsFloat y);
  ///Get the interpolated value at a given point
  rsFloat Value(rsFloat x);
  /// Get the maximum value in the set
  rsFloat Max() const;
  /// Divide the set by a given number
  void Divide(rsFloat a);
private:
  std::map<rsFloat, rsFloat> data;
};


///Load samples into the set
void InterpSetData::LoadSamples(const vector<rsFloat> &x, const vector<rsFloat> &y)
{
  vector<rsFloat>::const_iterator ix = x.begin();
  vector<rsFloat>::const_iterator iy = y.begin();
  for (; (ix != x.end()) && (iy != y.end()); ix++, iy++)
    {
      data.insert(pair<rsFloat, rsFloat>(*ix, *iy));
    }
}

///Load a single sample into the set
void InterpSetData::InsertSample(rsFloat x, rsFloat y)
{
  data.insert(pair<rsFloat, rsFloat>(x, y));
}

///Get the interpolated value for the given point
rsFloat InterpSetData::Value(rsFloat x)
{
  //Use linear interpolation, for now
  
  //If the map is empty, throw an exception and whine
  if (data.empty())
    throw std::logic_error("[BUG] Interpolation on an empty list in InterpSet");
  //Get the first element with a key greater than k
  map<rsFloat, rsFloat>::const_iterator iter = data.lower_bound(x);
  //If we are at the beginning of the set, return the value
  if (iter == data.begin())
    return (*iter).second;
  map<rsFloat, rsFloat>::const_iterator prev = iter;
  prev--;

  //If we are over the end, return the last value
  if (iter == data.end())
    return (*(prev)).second;
  //If we hit a sample exactly, return the value
  else if ((*iter).first == x)
    return (*iter).second;
  //Perform linear interpolation
  else {
    rsFloat x1 = (*prev).first;
    rsFloat x2 = (*iter).first;
    rsFloat y1 = (*prev).second;
    rsFloat y2 = (*iter).second;    
    return y2*(x-x1)/(x2-x1)+y1*(x2-x)/(x2-x1);
  }
}

/// Get the maximum value in the set
rsFloat InterpSetData::Max() const
{
  map<rsFloat, rsFloat>::const_iterator iter;
  rsFloat max = 0;
  // Scan through the map, updating the maximum
  for (iter = data.begin(); iter != data.end(); iter++) {
    if (std::fabs((*iter).second) > max)
      max = std::fabs((*iter).second);
  }
  return max;  
}

/// Divide the set by a given number
void InterpSetData::Divide(rsFloat a)
{
  map<rsFloat, rsFloat>::iterator iter;
  for (iter = data.begin(); iter != data.end(); iter++)
    (*iter).second /= a;
}

//
// InterpSet Implementation
//

/// Constructor
InterpSet::InterpSet()
{
  data = new InterpSetData();
}


/// Destructor
InterpSet::~InterpSet()
{
  delete data;
}

/// Load a number of samples into the set
void InterpSet::LoadSamples(const std::vector<rsFloat> &x, const std::vector<rsFloat> &y)
{
  data->LoadSamples(x,y);
}

/// Load a single sample into the set
void InterpSet::InsertSample(rsFloat x, rsFloat y)
{
  data->InsertSample(x,y);
}

/// Get the interpolated value at the given point
rsFloat InterpSet::Value(rsFloat x)
{
  return data->Value(x);
}

/// Get the maximum value in the set
rsFloat InterpSet::Max() const
{
  return data->Max();
}

/// Divide every sample in the set by a given number
void InterpSet::Divide(rsFloat a) {
  data->Divide(a);
}
