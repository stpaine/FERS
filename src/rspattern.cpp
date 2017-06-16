//rspattern.cpp
// Implementation of interpolated 2D arrays for gain patterns and RCS patterns
// Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//19 September 2007

#include <config.h>
#include <cmath>
#include "rspattern.h"
#include "rshdf5.h"

using namespace rs;

//
// Pattern Implementation
//

/// Constructor
Pattern::Pattern(const std::string &filename) 
{
  //Load the pattern from the HDF5 file
  pattern = rshdf5::ReadPattern(filename, "antenna", size_azi, size_elev);
}

/// Destructor
Pattern::~Pattern() 
{
  //Clean up
  for (unsigned int i = 0; i < size_azi; i++)
    delete[] pattern[i];
  delete[] pattern;
}

/// Get the gain at the given angle
rsFloat Pattern::GetGain(const rs::SVec3 &angle) const
{
  //Get the nearest points in elevation and azimuth
  double x1 = std::floor((angle.azimuth+M_PI)/(2*M_PI)*(size_azi-1))/(double)(size_azi-1);
  double ex1 = (angle.azimuth+M_PI)/(2*M_PI);
  double x2 = x1+1.0/(double)(size_azi);
  double y1 = std::floor((angle.elevation+M_PI)/(2*M_PI)*(size_elev-1))/(double)(size_elev-1);
  double ey1 = (angle.elevation+M_PI)/(2*M_PI);
  double y2 = y1+1.0/(double)(size_elev);
  //Get the interpolation constants
  double t = (ex1-x1)/(x2-x1);
  double u = (ey1-y1)/(y2-y1);
  
  //Get the offsets into the array
  int arr_x = std::floor(x1*size_azi);
  int arr_y = std::floor(y1*size_elev);
  
  //Get the interpolated value, using bilinear interpolation
  double interp = (1.0-t)*(1.0-u)*pattern[arr_x][arr_y];
  interp += t*(1.0-u)*pattern[(arr_x+1)%size_azi][arr_y];
  interp += t*u*pattern[(arr_x+1)%size_azi][(arr_y+1)%size_elev];
  interp += (1.0-t)*u*pattern[arr_x][(arr_y+1)%size_elev];
  return interp;
}
