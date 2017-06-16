/// rshdf5.cpp - Export binary data to the HDF5 file format
/// Marc Brooker mbrooker@rrsg.ee.uct.ac.za
/// 03 November 2006

#include <config.h>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iomanip>
#include "rshdf5.h"
#include "rsparameters.h"
#include "rsdebug.h"

using namespace rshdf5;

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

/// Read the pulse data from the specified file
void rshdf5::ReadPulseData(const std::string &name, std::complex<rsFloat> **data, unsigned int &size, rsFloat &rate)
{
  rate = rs::rsParameters::rate();
  //Open the HDF5 file
  hid_t file = OpenFile(name);
  //Get the size of the dataset "pulse"
  size_t type_size;
  H5T_class_t class_id;
  hsize_t* dims;
  //Open the / group
  hid_t slash = H5Gopen1(file, "/");
  if (slash < 0)
    throw std::runtime_error("[ERROR] HDF5 file "+name+" does not have top level group \"/\"");
  //Open the I group
  hid_t Igroup = H5Gopen1(slash, "I");
  if (Igroup < 0)
    throw std::runtime_error("[ERROR] HDF5 file "+name+" does not have group \"I\"");
  // Get the rank of the groups
  int rank;
  H5LTget_dataset_ndims(Igroup, "value", &rank); 
  dims = new hsize_t[rank];
  //Get the data set information
  herr_t res = H5LTget_dataset_info(Igroup, "value", &(dims[0]), &class_id, &type_size);
  if (res < 0)
    throw std::runtime_error("[ERROR] HDF5 file "+name+" does not have dataset \"value\" in group \"I\"");
  // Allocate memory for the pulse
  size = dims[0];
  double* buffer_I = new double[size];
  // Read in the I dataset
  res = H5LTread_dataset_double(Igroup, "value", buffer_I);
  if (res < 0)
    throw std::runtime_error("[ERROR] Error reading dataset I of file "+name);
  //Close the I group
  H5Gclose(Igroup);
  //Open the Q group
  hid_t Qgroup = H5Gopen1(slash, "Q");
  if (Qgroup < 0)
    throw std::runtime_error("[ERROR] HDF5 file "+name+" does not have group \"Q\"");
  // Check the Q dataset is the same size
  res = H5LTget_dataset_info(Qgroup, "value", &(dims[0]), &class_id, &type_size);
  if (res < 0)
    throw std::runtime_error("[ERROR] HDF5 file "+name+" does not have dataset \"Q\"");
  if (size != dims[0])
    throw std::runtime_error("[ERROR] Dataset \"Q\" is not the same size as dataset \"I\" in file " + name);
  //Allocate memory for the Q set
  double* buffer_Q = new double[size];
  //Read in the Q datase
  res = H5LTread_dataset_double(Qgroup, "value", buffer_Q);
  if (res < 0)
    throw std::runtime_error("[ERROR] Error reading dataset Q of file "+name);
  // Close group q
  H5Gclose(Qgroup);
  H5Gclose(slash);
  //Close the HDF5 file
  H5Fclose(file);
  //Copy the data to the target array
  *data = new std::complex<rsFloat>[size];
  for (unsigned int i = 0; i < size; i++)
    (*data)[i] = std::complex<rsFloat>(buffer_I[i], buffer_Q[i]);
  // Clean up memory
  delete[] buffer_I;
  delete[] buffer_Q;  
  delete[] dims;
}

///Open the HDF5 file for writing
int rshdf5::CreateFile(const std::string& name)
{
  hid_t file = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0)
    throw std::runtime_error("[ERROR] Could not create HDF5 file "+name+" for export");
  return static_cast<int>(file);
}

///Add a dataset to the HDF5 file
void rshdf5::AddChunkToFile(int file, std::complex<rsFloat> *data, unsigned int size, rsFloat time, rsFloat rate, rsFloat fullscale, unsigned int count)
{
  //Create the name of the dataset
  std::ostringstream oss;
  oss << "chunk_" << std::setw(6) << std::setfill('0') << count;
  std::string I_chunk_name = oss.str()+"_I";
  std::string Q_chunk_name = oss.str()+"_Q";
  //Create the size variable needed by the lite api
  hsize_t datasize = size;
  //Write out the I data
  double *I = new double[size];
  double *Q = new double[size];
  //Seperate I and Q
  for (unsigned int i = 0; i < size; i++) {
    I[i] = data[i].real();
    Q[i] = data[i].imag();
  }
  //Create the dataset, using the HDF5 Lite API
  if (H5LTmake_dataset_double(static_cast<hid_t>(file), I_chunk_name.c_str(), 1, &datasize, I) < 0)
    throw std::runtime_error("[ERROR] Error while writing data to HDF5 file");  
  if (H5LTmake_dataset_double(static_cast<hid_t>(file), Q_chunk_name.c_str(), 1, &datasize, Q) < 0)
    throw std::runtime_error("[ERROR] Error while writing data to HDF5 file");

  //Add attributes to the data set, with the attributes of the response
  if (H5LTset_attribute_double(static_cast<hid_t>(file), I_chunk_name.c_str(), "time", &time, 1) < 0)
    throw std::runtime_error("[ERROR] Error while setting attribute \"time\" on chunk " + I_chunk_name);
  if (H5LTset_attribute_double(static_cast<hid_t>(file), I_chunk_name.c_str(), "rate", &rate, 1) < 0)
    throw std::runtime_error("[ERROR] Error while setting attribute \"rate\" on chunk " + I_chunk_name);
  if (H5LTset_attribute_double(static_cast<hid_t>(file), I_chunk_name.c_str(), "fullscale", &fullscale, 1) < 0)
    throw std::runtime_error("[ERROR] Error while setting attribute \"fullscale\" on chunk " + I_chunk_name);
  if (H5LTset_attribute_double(static_cast<hid_t>(file), Q_chunk_name.c_str(), "time", &time, 1) < 0)
    throw std::runtime_error("[ERROR] Error while setting attribute \"time\" on chunk " + Q_chunk_name);
  if (H5LTset_attribute_double(static_cast<hid_t>(file), Q_chunk_name.c_str(), "rate", &rate, 1) < 0)
    throw std::runtime_error("[ERROR] Error while setting attribute \"rate\" on chunk " + Q_chunk_name);
  if (H5LTset_attribute_double(static_cast<hid_t>(file), Q_chunk_name.c_str(), "fullscale", &fullscale, 1) < 0)
    throw std::runtime_error("[ERROR] Error while setting attribute \"fullscale\" on chunk " + Q_chunk_name);

  //Free the buffers
  delete[] I;
  delete[] Q;  
}

///Close the HDF5 file
void rshdf5::CloseFile(int file) {
  if (H5Fclose(static_cast<hid_t>(file)) < 0)
    throw std::runtime_error("[ERROR] Error while closing HDF5 file");
}

 /// Read an antenna gain pattern or RCS pattern from a file
rsFloat **rshdf5::ReadPattern(const std::string &name, const std::string &dataset_name, unsigned int &azi_size, unsigned int &elev_size)
{
  hid_t file_id; 
  int rank;
  hsize_t dims[2];
  size_t type_size;
  H5T_class_t data_class;
  //Load the HDF5 file
  file_id = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0)
    throw std::runtime_error("[ERROR] Cannot open HDF5 file "+name + " to read antenna data");
  // Get the rank of the dataset
  herr_t err = H5LTget_dataset_ndims(file_id, dataset_name.c_str(), &rank);
  if (err < 0)
    throw std::runtime_error("[ERROR] Could not get rank of dataset \""+dataset_name+"\" in file " + name);
  else if (rank != 2)
    throw std::runtime_error("[ERROR] Dataset \""+dataset_name+"\" in file "+name+" does not have rank 2");
  //Get the dimensions of the file
  err = H5LTget_dataset_info(file_id, dataset_name.c_str(), &(dims[0]), &data_class, &type_size);
  if (err < 0)
    throw std::runtime_error("[ERROR] Could not get dimensions of dataset \""+dataset_name+"\" in file " + name);
  if (type_size != sizeof(float))
    throw std::runtime_error("[ERROR] Type size incorrect in dataset \""+dataset_name+"\" in file " + name);
  // Allocate memory for the pattern
  float *data = new float[dims[0]*dims[1]];
  /// Load the pattern into memory
  err = H5LTread_dataset_float(file_id, dataset_name.c_str(), data);
  if (err < 0) {
    delete[] data;
    throw std::runtime_error("[ERROR] Could not read float data from dataset \""+dataset_name+"\" in file" + name);
  }
  ///Close the HDF5 file
  if (H5Fclose(file_id) < 0)
    throw std::runtime_error("[ERROR] Error while closing HDF5 file "+name);
  /// Break the data down into a 2D array
  azi_size = dims[0];
  elev_size = dims[1];
  rsFloat **ret = new rsFloat*[azi_size];
  for (unsigned int i = 0; i < azi_size; i++) {
    ret[i] = new rsFloat[elev_size];
    for (unsigned int j = 0; j < elev_size; j++) 
      ret[i][j] = data[i*azi_size+j];
  }
  //Clean up
  delete[] data;
  return ret;
}
