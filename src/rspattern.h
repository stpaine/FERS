//rspattern.h
// Interpolated 2D arrays for gain patterns and RCS patterns
// Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//11 September 2007

#ifndef __RS_PATTERN_H
#define __RS_PATTERN_H

#include <string>
#include "rsgeometry.h"

namespace rs {

  ///Class to manage and interpolate a 2D pattern
  class Pattern {
  public:
    /// Constructor
    Pattern(const std::string &filename);
    /// Destructor
    ~Pattern();
    /// Get the gain at the given angle
    rsFloat GetGain(const rs::SVec3 &angle) const;
  private:
    unsigned int size_elev, size_azi; //!< Number of samples in elevation and azimuth
    rsFloat **pattern; //!< 2D Array to store the gain pattern
  };

}

#endif
