//rspython.cpp - Implementation of Python extensions
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//7 March 2007


#include <Python.h>
#include "rspython.h"
#include <stdexcept>
#include "rsdebug.h"

using namespace rsPython;

/// Keep track of whether python has been globally initialized
bool initialized;


///Initialize the python environment
void rsPython::InitPython()
{
  if (!initialized) {
    Py_Initialize();
    rsDebug::printf(rsDebug::RS_VERBOSE, "Using Python version %s\n", Py_GetVersion());
    PyRun_SimpleString("import sys; sys.path.append('.');");
    initialized = true;
  }
}


// Define the PythonExtensionImpl class
struct rsPython::PythonExtensionData {
  PyObject *pModule, *pFunc;
};

//
// PythonExtension Implementation
//


/// Constructor
PythonExtension::PythonExtension(const std::string& module, const std::string& function):
  module(module),
  function(function)
{
  data = new PythonExtensionData;
  //Import the module
  PyObject *modname = PyString_FromString(module.c_str());
  data->pModule = PyImport_Import(modname);
  Py_DECREF(modname);

  if (!data->pModule) {
    PyErr_Print();
    throw std::runtime_error("Could not import Python module "+module);
  }
  //Import the function
  PyObject *funcname = PyString_FromString(function.c_str());
  data->pFunc = PyObject_GetAttr(data->pModule, funcname);
  Py_DECREF(funcname);
  if (!data->pFunc) {
    PyErr_Print();
    throw std::runtime_error("Could not import Python function "+function+" from module "+module);
  }
  if (!PyCallable_Check(data->pFunc))
    throw std::runtime_error("Python object "+function+" from module "+module+" is not callable");
}

/// Destructor
PythonExtension::~PythonExtension()
{
  Py_DECREF(data->pFunc);
  Py_DECREF(data->pModule);
  delete data;
}

//
// PythonPath Implementation
//

/// Constructor
PythonPath::PythonPath(const std::string& module, const std::string& function):
  PythonExtension(module, function)
{
}

/// Destructor
PythonPath::~PythonPath()
{
}

/// Get the position at a given time
rs::Vec3 PythonPath::GetPosition(rsFloat t) const
{
  //Insert t into a tuple for passing to the function
  PyObject *pargs = PyTuple_Pack(1, PyFloat_FromDouble(t));
  //Check that the Tuple was created successfully
  if (!pargs)
    throw std::runtime_error("Could not create new Python Tuple in PythonPath");
  //Call the function
  PyObject *result = PyObject_CallObject(data->pFunc, pargs);
  //Release the arguments
  Py_DECREF(pargs);
  //Check the call was completed
  if (!result) {
    PyErr_Print();
    throw std::runtime_error("Call of function "+function+" from module "+module+" failed");
  }
  //Unpack the results into Vec3
  PyObject *px = PyTuple_GetItem(result, 0);
  PyObject *py = PyTuple_GetItem(result, 1);
  PyObject *pz = PyTuple_GetItem(result, 2);
  if ((!px) || (!py) || (!pz)) {
    PyErr_Print();
    throw std::runtime_error("Python function "+function+" did not return expected tuple");
  }
  rs::Vec3 vec(PyFloat_AsDouble(px), PyFloat_AsDouble(py), PyFloat_AsDouble(pz));
  Py_DECREF(result);
  return vec; 
}


//
//PythonNoise Implementation
//

///Constructor
PythonNoise::PythonNoise(const std::string& module, const std::string& function):
  PythonExtension(module, function)
{

}

///Destructor
PythonNoise::~PythonNoise()
{

}

/// Get a noise sample
rsFloat PythonNoise::GetSample() const
{
  //Call the function
  PyObject *result = PyObject_CallObject(data->pFunc, 0);
  //Check the call was completed
  if (!result) {
    PyErr_Print();
    throw std::runtime_error("Call of function "+function+" from module "+module+" failed");
  }
  //Unpack the results
  rsFloat sample = PyFloat_AsDouble(result);
  Py_DECREF(result);
  return sample; 
}

//
// PythonAntennaMod Implementation
//

///Constructor
PythonAntennaMod::PythonAntennaMod(const std::string& module, const std::string& function):
  PythonExtension(module, function)
{

}

///Destructor
PythonAntennaMod::~PythonAntennaMod()
{

}

/// Get the antenna gain in the specified direction
rsFloat PythonAntennaMod::GetGain(const rs::SVec3& direction) const
{
  //Insert t into a tuple for passing to the function
  PyObject *pargs = PyTuple_Pack(2, PyFloat_FromDouble(direction.azimuth), PyFloat_FromDouble(direction.elevation));
  //Check that the Tuple was created successfully
  if (!pargs)
    throw std::runtime_error("Could not create new Python Tuple in PythonPath");
  //Call the function
  PyObject *result = PyObject_CallObject(data->pFunc, pargs);
  //Release the arguments
  Py_DECREF(pargs);
  //Check the call was completed
  if (!result) {
    PyErr_Print();
    throw std::runtime_error("Call of function "+function+" from module "+module+" failed");
  }
  //Translate the result
  rsFloat sample = PyFloat_AsDouble(result);
  Py_DECREF(result);
  return sample; 
}
