# FERS - The Flexible Extensible Radar Simulator

## NOTE: PLEASE REFER TO THE GIT WIKI FOR FURTHER INFORMATION

---

FERS is a simulator for simulating the performance and output of a variety of radar systems. It is designed to support a wide range of traditional and modern radar system designs. FERS is currently under active development - this version is likely to contain many bugs, incomplete or missing features, incomplete or missing documentation and major inaccuracies.

FERS has been used to generate useful results for real-world projects, so it might be useful for your projects. Documentation is currently very sparse - it will be improved soon.

The features which are currently implemented are:

* Creation of returns from arbitrary pulse shapes
* Simple propagation, doppler and phase model
* Modelling of Multistatic and Monostatic radar systems
* Modelling of CW and Pulsed radars
* Export to CSV, XML and HDF5 file formats
* Proper range-gate and timing for pulsed radars
* Modelling of 1/f noise on local oscillators
* Effects of multipath propagation

#### AUTHORS:

FERS was written by:

Marc Brooker (marcbrooker@gmail.com)

FERS is currently maintained by the following people:

Craig Tong (craig.tong@uct.ac.za)

#### BUILDING FERS

FERS depends on a number of external libraries, which you need to install before attempting to build FERS. The libraries you need to have installed are:

* Boost (at least boost::threads, boost::system and boost::random header) which is available freely from http://www.boost.org
* FFTW3 which is available freely from http://www.fftw3.org
* libhdf5
* TinyXML which is available freely from http://sourceforge.net/projects/tinyxml/

On a Debian or Ubuntu system (or pretty much any other decent GNU/Linux distribution) these libraries should be available pre-packaged for your installing pleasure.

FERS also depends on the cmake system. This system can generate makefiles, VC++ projects, KDevelop projects and more from the FERS sources. Please install cmake before attempting to build FERS.

* Download and extract the tar ball.
* Navigate to the extracted directory (by default: `cd fers/`)
* Enter the following commands:
    ```bash
    mkdir build
    cd build
    cmake ../
    make
    ```

-------------------------------------------------------------------------------------
## NOTE: PLEASE REFER TO THE GIT WIKI FOR FURTHER INFORMATION ON SOLVING COMPILE ERRORS
-------------------------------------------------------------------------------------

A "fers" binary will be then be placed in the "src" directory which can be copied to a location of your choice.

FERS can be build in debug mode by replacing the "cmake ../" with cmake -DCMAKE_BUILD_TYPE=Debug ../

Note it is advisable use separate build folders for release and debug builds.

FERS is written in standard C++ and should compile and run on many architectures and operating systems. Please report successes and failures of running this software on non-Linux and non-x86 platforms to the authors, so we can improve the software and make it more portable.

#### DOCUMENTATION

Documentation is available in the doc/ directory. Highlights include:

doc/equations/equations.tex - All equations used by FERS in convenient LaTeX form
fersxml.dtd - Document Type Definition for the XML script file format

#### THANKS

The authors of tinyXml, FFTW and Boost for making excellent software freely available.

#### COPYRIGHT NOTICE

FERS is covered by the following copyright notice. Should you wish to acquire a copy of FERS not covered by these terms, please contact the Department of Electrical Engineering at the University of Cape Town.

Please note that this copyright only covers the source code, program binaries and build system of FERS. Any input files you create and results created by the simulator are not covered by this notice and remain copyright of their original author.

This copyright notice does not cover the source code of tinyxml (found in the tinyxml directory). Please see the file readme.txt in that directory for a copyright notice for that code. TinyXML can, however, be freely distributed along with the code of FERS.

FERS - The Flexible Extensible Radar Simulator
Copyright (C) 2006 Marc Brooker and Michael Inggs

This program is free software; you can redistribute it and/or modify
it under the terms of version 2 of the GNU General Public License as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA


```
cmake -D FERS_LIB_HDF5="/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so" -D FERS_LIB_HDF5_HL="/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.so" -D CMAKE_CXX_FLAGS="-I/usr/include/hdf5/serial/" ../
```