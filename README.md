# Ichnos
Ichnos is a kind of general purpose particle tracking code. As such Ichnos can trace virtual particles in simulated velocity fields. However the development of the code is heavily influenced by the projects that Ichnos developed for which is particle tracking in groundwater fields. In particular in groundwater fields that derived from finite element groundwater models. Nevertheless groundwater velocity fields from finite difference and finite volume codes can alos be used. A key consideration during the development of the code is to handle simulations multi-million velocity nodes. Therefore the code can be used with velocitiy fields that are split accross many processors.

In particle tracking codes there are essentially two main functionalities. i) Tracing the particles in a velocity field. ii) Interpolating a velocity field. The code has been designed so that tracing functions are unaware of the details of the velocity field. To achieve that there is one velocity base class that tracing function is aware. Then each different type of field has to derive from the base velocity class and overwrite few functions. 

At the moment we provide one type of velocity field i.e. point cloud type and three modes

  1. Steady State flow fields
  2. Transient State flow fields (Coming soon...)
  3. Stochastic flow fields (This is experimental)

# Getting the code
For Windows there are executable available here. However those executables may not contain the latest updates, fixes and new bugs.
For linux or getting the latest version under windows one need to build the project    


# Building Ichnos
To build ichnos we use [cmake](https://cmake.org/). The code contains the required [CMakeList.txt](https://github.com/giorgk/ichnos/blob/master/ichnos_code/CMakeLists.txt). However it is possible that it may need adjustments for each system.

Ichnos depends on a few libraries

* [MPI] All mpi frameworks should work. 
* [Boost](https://www.boost.org/) witn mpi support. Any relatively new version of boost is good. In the CMakeList.txt we require 1.58 but the code works with lower versions.
* [CGAL](https://www.cgal.org/) Similarly any relatively recent version should work

## Bulding the dependencies
Gettting the dependencies right sometimes can be quite tricky. To build Ichnos we have used three options for different systems.

* [vcpkg](https://github.com/microsoft/vcpkg) is probably the easiset method to build the dependencies. It has always worked under windows without issues. Under linux can work but not without issues.
* [Spack](https://spack.io/) Using spack is probably the best option for linux, yet is a bit more involved compared th vcpkg. There are more steps about spack in the following
* **Build everything** Thats requires some experience and is recommnded on system without administration provilages, e.g clusters

## Building with vcpkg
All that is needed is to pass the following argument
```
-DCMAKE_TOOLCHAIN_FILE=path\to\vcpkg\scripts\buildsystems\vcpkg.cmake
```

## Building with Spack
### Prepare the environment
First create an environment that we name it ichnos
```
spack env create ichnos
```
and make it active
```
spack env activate ichnos
```
### CMake
First we are going to install cmake. Any version should work however I'm specifing here 3.9.4 because this has already been built in other spack environments
 ```
 spack install cmake@3.9.4 target=x86_64
 ```
 Now the command `spack find` should find that the cmake is installed anc the `cmake --version` should return the version information.
 
 ### MPI
 Next we will install openmpi. The `spack info openmpi` will return the available openmpi versions, while at the top specifies which one is the prefered. At the time of writing this that was the 3.1.4 so we build it as
 ```
 spack install openmpi@3.1.4 target=x86_64
 ```
 Usually clusters provide mpi frameworks and it is always a good practise to use these. To instruct spack to use an external package ones has to create a `packages.yaml` file with the content similar to the following
 ```
packages:
  mpi:
    buildable: False
  openmpi:
    modules:
      openmpi@4.0.1: openmpi/4.0.1
    buildable: False
 ```
 It is important to set for the external packages the buildable flag as false. 
 
 *Note: It is possible to build your own mpi via spack however there is the chance that you get errors during runtime that mpi is linked with a different version of mpi*
  
### Boost
 The next piece ichnos needs is boost. I'll install 1.69 as this is the one already exists from previous projects. 
 ```
 spack install boost@1.69.0+mpi ^openmpi@3.1.4 target=x86_64
 ```
 ### CGAL
 CGAL is the next library that is required by ichnos. CGAL depends on boost so we will specify that we want to build cgal using the boost version we just installed
 ```
 spack install cgal@4.13 ^boost@1.69.0+mpi ^openmpi@3.1.4 target=x86_64
 ```
 Sometimes cgal failed to get build because the mprf failed to apply a patch. One way to bypass this is by using a different version of mpfr
 ```
 spack install cgal@4.13 ^boost@1.69.0+mpi ^openmpi@3.1.4 ^mpfr@4.0.0
 ```
 
 ### Nanoflann (Skip this nanoflann is no longer needed)
 The last required library is the nanoflann. THis is just a header that has no dependancies so all is needed is
 ```
 spack install nanoflann target=x86_64
 ```
 ### Building all
 Assuming that the current directory is `../ichnos_code/` where the CMakeLists.txt file is, create a folder, get into the folder and run cmake as 
 ```
 mkdir build
 cd build
 cmake ..
 ```

 If cmake finds the wrong boost version e.g the one under /usr/ for example then passing the following during cmake might help 
 ```
cmake -DBoost_NO_BOOST_CMAKE=TRUE \
      -DBoost_NO_SYSTEM_PATHS=TRUE \
      -DBOOST_ROOT:PATHNAME=/path/to/spack/var/spack/environments/ichnos/.spack-env/view/ ..
 ```
 see [more](https://stackoverflow.com/questions/3016448/how-can-i-get-cmake-to-find-my-alternative-boost-installation).


 ### Building in aqua cluster
 This guide has not been updated but contains usefull information
 First clean the modules and load the openmpi
 ```
module purge
module load system-gcc/openmpi-2.1.5
```
## Dependencies 
The code depends on three libraries boost, nanoflann and CGAL.
### Boost 
At the moment of writing aqua cluster has boost 1.58 available. However the default CMakeLists.txt file requests for 1.69. Therefore we have to edit the file adn downdgrade the requested version
### nanoflann
Although this is a header only library its a good practice to built it so that we can use the find_package function of cmake. To get build and install nanoflann do the following:
```
git clone https://github.com/jlblancoc/nanoflann.git
mkdir install
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/nanoflann/install ..
make install
```

Inside the ichnos_code folder  create a build directory
```
cd ichnos_code
mkdir build
cd build
```
Next run cmake to configure the project. Since this is for final runs we should use the optimize version.
```
cmake -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR:PATH=/path/to/cgal-releases-CGAL-4.11.3/build/release ..
```

At the time of writing, aqua cluster has boost verion 1.58 while the cmake requires the 1.69. to bypass this just downgrade the requested version.


# Running Ichnos
Ichnos reads the inputs from ascii files and writes the outputs to ascii files as well. </br>
To run Ichnos one need to prepare two configuration files as well as a number of other input files that descibe the velocity fields and the domain.
You can run the code with the following options:

* Read the version of the executable
```
ichnos.exe -v
```
* Get help options
```
ichnos.exe -h
```

* Run the code
```
ichnos.exe -c config.ini
```

