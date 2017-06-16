// Turn FERS HDF5 output into a raw binary file with 32 bit interleaved float samples
// This was used for integeration with the G2 SAR processor, and probably isn't a good idea in most cases - the HDF5 version is much easier to work with
// Marc Brooker mbrooker@rrsg.ee.uct.ac.za
// 18 April 2008

#include <sstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <iomanip>

#include <stdio.h>

extern "C" {
#include <hdf5.h>
#include <H5LTpublic.h>
}

///Open the HDF5 file for reading
hid_t OpenFile(const std::string &name)
{
  hid_t file = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0)
    throw std::runtime_error("[ERROR] Could not open HDF5 file "+name+" to read pulse");
  return file;
}


/// Read the HDF5 datasets and dump them out the the correct format
void ReadAndDump(hid_t file, FILE *outfile)
{
  size_t size = 0;
  unsigned int i;
  double maxI = 0, maxQ = 0;
  //Open the / group
  hid_t slash = H5Gopen1(file, "/");
  if (slash < 0)
    throw std::runtime_error("[ERROR] HDF5 file does not have top level group \"/\"");
  // Count the chunks
  size_t count = 0;
  for (i = 0; i < 100000; i++) {
    std::ostringstream oss;
    oss << "chunk_" << std::setw(6) << std::setfill('0') << i;
    std::string I_chunk_name = oss.str()+"_I";
    std::string Q_chunk_name = oss.str()+"_Q";
    if (H5LTfind_dataset(slash, I_chunk_name.c_str()) < 1)
      break;
        double Iscale, Qscale;
    H5LTget_attribute_double(slash, I_chunk_name.c_str(), "fullscale", &Iscale);
    H5LTget_attribute_double(slash, Q_chunk_name.c_str(), "fullscale", &Qscale);
    if (Iscale > maxI)
      maxI = Iscale;
    if (Qscale > maxQ)
    maxQ = Qscale;
    if (i == 0) {
      hsize_t* dims;
      // Get the dataset ranks
      int rank;
      if (H5LTget_dataset_ndims(slash, I_chunk_name.c_str(), &rank) < 0)
        break;
      dims = new hsize_t[rank];
      // Get the dataset info
      size_t type_size;
      H5T_class_t class_id;
      if (H5LTget_dataset_info(slash, I_chunk_name.c_str(), &(dims[0]), &class_id, &type_size) < 0) {
        std::cerr << "Could not get dataset info for " << I_chunk_name << std::endl;
        break;
      }
      // Allocate memory for the I and Q data
      size = dims[0];
    }    
  }
  std::cout << "MaxI " << maxI << " maxQ " << maxQ << std::endl;
  count = i;
  // Allocate memory for the reshuffle
  float* buffer = new float[count*size*2];
  for (i = 0; i < count; i++) {
    // Get the dataset names
    std::ostringstream oss;
    oss << "chunk_" << std::setw(6) << std::setfill('0') << i;
    std::string I_chunk_name = oss.str()+"_I";
    std::string Q_chunk_name = oss.str()+"_Q";
    if (H5LTfind_dataset(slash, I_chunk_name.c_str()) < 1)
      break;
    double* buffer_I = new double[size];
    double* buffer_Q = new double[size];
    H5LTread_dataset_double(slash, I_chunk_name.c_str(), buffer_I);
    H5LTread_dataset_double(slash, Q_chunk_name.c_str(), buffer_Q);
    // Get the fullscale for I and Q
    double Iscale, Qscale;
    H5LTget_attribute_double(slash, I_chunk_name.c_str(), "fullscale", &Iscale);
    H5LTget_attribute_double(slash, Q_chunk_name.c_str(), "fullscale", &Qscale);
    // Write out the data
    for (unsigned int j = 0; j < size; j++) {
      float I = buffer_I[j]*Iscale/maxI;
      float Q = buffer_Q[j]*Qscale/maxQ;
      //fwrite(&I, 1, 1, outfile);
      //fwrite(&Q, 1, 1, outfile);
      buffer[(j+i*size)*2] = I;
      buffer[(j+i*size)*2+1] = Q;
    }
    //Clean up
    delete[] buffer_I;
    delete[] buffer_Q;    
  }
  fwrite(buffer, count*size*2, sizeof(float), outfile);
  std::cout << "Read " << i-1 << " windows of length " << size << std::endl;
  delete[] buffer;
}

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: h52raw infile outfile\n";
    return 1;
  }
  hid_t infile = OpenFile(argv[1]);
  if (infile < 0) {
    std::cerr << "Could not open HDF5 file " << argv[1] << std::endl;
    return 1;
  }
  FILE *outfile = fopen(argv[2], "w");
  if (!outfile) {
    std::cerr << "Could not open file " << argv[2] << std::endl;
    return 1;
  }
  //Perform the reading and dumping
  ReadAndDump(infile, outfile);
  //Close the files
  H5Fclose(infile);
  fclose(outfile);
}
