# Ichnos overview
Ichnos is a general purpose particle tracking code. As such, Ichnos can trace virtual particles in velocity fields. The development of the code is heavily influenced by the projects that Ichnos was developed for, which is primarily particle tracking in groundwater flow fields. In particular, in groundwater flow fields calculated from finite element groundwater models. Yet, groundwater velocity fields from finite difference and finite volume codes can also be used. A key consideration during the development of the code is to handle simulations with multi-million velocity nodes. Therefore the code can be used with velocitiy fields that can be split across many processors.

A brief demonstration of the capabilities of this code can be found at our [AGU 2020 poster](https://agu2020fallmeeting-agu.ipostersessions.com/?s=0C-96-C6-05-F8-AB-22-9D-63-C2-37-9D-96-2E-CD-B8).

Note that this is a continuation of a previous project [IWFMtrack](https://gwt.ucdavis.edu/research-tools-and-applications/iwfm-track).

In particle tracking codes there are essentially two main functionalities. i) Tracing the particles in a velocity field. ii) Interpolating a velocity field. The code has been designed so that tracing functions are unaware of the details of the velocity field. To do so the code contains a velocity base class that communicates with the tracing class. Then each different type of field has to derive from the base velocity class and overwrite a few functions. The velocity field can be split into 2 components. The positional information of the points where the velocity is known and the velocity itself. Each of those components are defined in based classes. The positional base class `class XYZ` is responsible to calculate the weight contribution of each velocity point in the field for a given particle position. Then it passes the weights to the `Velocity` class which calculates the velocity value. Using this separation it is possible to combine different data structures in the code.

At the moment we provide two types of positional classes `XYZ` 
  1. **CLOUD** </br>
  The CLOUD type assumes that the velocity information is defined into a cloud of points without any particular relation other that proximity between the nodes.</br>
  
  2. **MESH2D** </br>
  The MESH2D type assumes that the velocity is organized via a 2D mesh that is extruded in the vertical direction. Within the extruded 2D mesh the velocity can be defined at the element barycenters, which assumes constant velocity within the elements, or at the mesh nodes, or at the Mesh faces. In the latter case the interpolation is carried out using the [Raviart–Thomas basis functions](https://en.wikipedia.org/wiki/Raviart%E2%80%93Thomas_basis_functions).
 
The velocity can be defined as:
  * Transient or Steady state
  * Stochastic flow fields (This is in an experimental stage)
  * Random Walk (This is in an experimental stage)

------------
# Outline of this repository
Here we breifly outline the Ichnos repository structure and what information each of these folders contain.
* **Bin** is the folder which holds the windows executable of ichnos. 
* **IchnosRwrksp** contains R scripts related to ichnos. (Currently we are using primarily Matlab threfore there is no much R scripting there)
* **WikiExamples** contains the configuration files for the examples presented at the [wiki](https://github.com/giorgk/ichnos/wiki)
* **docs** Was supposed to hold the doxygen documentation for the code. Currently is empty as we use codeblocks 
* **ichnos_code** contains the c++ code files.
* **ichnos_hou** contains project files and data for various [Houdini visualizations](https://www.sidefx.com/)

--------------
# Getting the code
The code relies on libraries that are available to all major operating OS.
## Windons
Windows executables are available under the [Bin](Bin/) folder. However, those executables may not contain the latest updates, fixes and new bugs.
Under the Bin folder there is a [Dev](Bin/Dev) folder which contains the latest build. 
<!-- You can check [here](https://gist.github.com/giorgk/f727351f5efd8c58f8bb885a87f41978#file-ichnos_test_list-md) which examples have been tested so far with the current development version. -->

Put under the same folder the `ichnos.exe` and all the `*.dlls`. Make sure the antivirus is not blocking the execution. To verify that the code can run try the following under a powershell or command prompt.
```
ichnos.exe -v
```
## Linux
For linux you have to compile the code. See the section [Builinding Ichnos](https://github.com/giorgk/ichnos?tab=readme-ov-file#building-ichnos) below.    

------------
# Get started
The [wiki](https://github.com/giorgk/ichnos/wiki) is the best source to get you started.

## Suggested Roadmap
Because the wiki contains alot of information and examples it can be a bit chaotic. Here we provide a roadmap on how to go through the wiki documentation

* The main page of [wiki](https://github.com/giorgk/ichnos/wiki) is very important as it provides a list of the required information of the two  required configuration files. It lists also all possible options. Therefore it is strongly suggested to go through this page to get an idea of the available options and requirements for a simulation run. This is also the page you would need to consult when setting up the two configuration files  

* Next it is recommended to go through, first the [simple example](https://github.com/giorgk/ichnos/wiki/Simple-Example) and seconldy the [sample face example](https://github.com/giorgk/ichnos/wiki/Simple-Face-Example) excersizes as they go through all required file preparation that one would need to do for an Ichnos simulation. 

* Then it is recommended to go thought the tutorials that run on a single core. The single core examples can also be run as multicore by tracing particles concurrently.

* Multiple core simulations where the domain is split into multiple subdomains require an extra effort and therefore one need to have a good understanding of the file requirements.

* Besides the two configuration files ichnos requires a specific format about how to write any information. The [File format](https://github.com/giorgk/ichnos/wiki/File-formats) should be used as a reference to the format of the files.  


------------

# Building Ichnos
To build ichnos we use [cmake](https://cmake.org/). The code contains the required [CMakeList.txt](https://github.com/giorgk/ichnos/blob/master/ichnos_code/CMakeLists.txt). However it is possible that it may need adjustments for each system.

Ichnos depends on a few libraries

* MPI. All mpi frameworks should work. 
* [Boost](https://www.boost.org/) with mpi support. Any relatively new version of boost is good. In the CMakeList.txt we require 1.58 but the code works with lower versions.
* [CGAL](https://www.cgal.org/) Similarly any relatively recent version should work.
* [HighFive](https://github.com/BlueBrain/HighFive). This is optional but speeds up the reading of large files and it is highly recommended. 

## Bulding the dependencies
Getting the dependencies right sometimes can be quite tricky. To build Ichnos we have used three options for different systems.

* [vcpkg](https://github.com/microsoft/vcpkg) is probably the easiest method to build the dependencies. It has always worked under windows without issues. Under linux can work but not without issues.
* [Spack](https://spack.io/) Using spack is probably the best option for linux, yet is a bit more involved compared to the vcpkg. There is more info about spack in the next section
* **Build everything** That requires some experience and it is recommended on systems without administration privileges, e.g clusters.

## Building with vcpkg
First install all the libraries that Ichnos depend on. Make sure you install the boost with mpi</br> 
Then pass the following argument to cmake
```
-DCMAKE_TOOLCHAIN_FILE=path\to\vcpkg\scripts\buildsystems\vcpkg.cmake -DUSEHF=True
```
If [HighFive](https://github.com/BlueBrain/HighFive) is not available then set `-DUSEHF=False`

## Building with Spack
The following is just a guide which at times is not working and requires workarounds.

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
Next we are going to install cmake. Any version should work however I'm specifying here 3.9.4 because this has already been built in other spack environments.
 ```
 spack install cmake@3.9.4 target=x86_64
 ```
 Now the command `spack find` should find that the cmake is installed and the `cmake --version` should return the version information.
 
 ### MPI
 Next we will install openmpi. The `spack info openmpi` will return the available openmpi versions, while at the top specifies which one is the prefered. At the time of writing this that was the 3.1.4 so we build it as
 ```
 spack install openmpi@3.1.4 target=x86_64
 ```
 Usually clusters provide mpi frameworks and it is always a good practice to use these. To instruct Spack to use an external package one has to create a `packages.yaml` file with content similar to the following:
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
 
 ### Highfive
 The last library Highfive is optional. 
 ```
 spack install highfive target=x86_64
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
cmake -DUSEHF=True \
      -DBoost_NO_BOOST_CMAKE=TRUE \
      -DBoost_NO_SYSTEM_PATHS=TRUE \
      -DBOOST_ROOT:PATHNAME=/path/to/spack/var/spack/environments/ichnos/.spack-env/view/ ..
 ```
 see [more](https://stackoverflow.com/questions/3016448/how-can-i-get-cmake-to-find-my-alternative-boost-installation).

------------


# Running Ichnos
More information about running Ichnos can he found in the [wiki section](https://github.com/giorgk/ichnos/wiki#run-ichnos). 

------------
# Who do I talk to
It depends on what you want to say</br>
The [Discussion](https://github.com/giorgk/ichnos/discussions) section is the best place to ask questions/clarifications regarding the use of Ichnos. </br>
Any issues with the code should be reported at [Issues](https://github.com/giorgk/ichnos/issues) section.

For other information you can reach us via email 
- gkourakos@ucdavis.edu
- hdahlke@ucdavis.edu 
- thharter@ucdavis.edu
