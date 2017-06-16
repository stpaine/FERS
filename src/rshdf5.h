/// rshdf5.h - Header file for HDF5 export functions
/// Marc Brooker mbrooker@rrsg.ee.uct.ac.za
/// 03 November 2006

#ifndef __RS_HDF5_H
#define __RS_HDF5_H

#include <complex>

namespace rshdf5 {

///Open the HDF5 file for writing
int CreateFile(const std::string& name);

///Add a dataset to the HDF5 file
void AddChunkToFile(int file, std::complex<rsFloat> *data, unsigned int size, rsFloat time, rsFloat rate, rsFloat fullscale, unsigned int count);

///Close the HDF5 file
 void CloseFile(int file);

/// Read the pulse data from the specified file
 void ReadPulseData(const std::string &name, std::complex<rsFloat> **data, unsigned int &size, rsFloat &rate);

 /// Read an antenna gain pattern or RCS pattern from a file
 rsFloat **ReadPattern(const std::string &name, const std::string &dataset_name, unsigned int &azi_size, unsigned int &elev_size);

}

#endif
