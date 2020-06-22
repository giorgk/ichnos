# ichnos

# Dependencies using Spack
## Prepare the environment
First create an environment that we nname it ichnos
```
spack env create ichnos
```
and make it active
```
spack env activate ichnos
```
## CMake
First we are going to install cmake. Any version should work however I'm specifing here 3.9.4 because this has already been built in other spack environments
 ```
 spack install cmake@3.9.4 target=x86_64
 ```
 Now the command `spack find` should find that the cmake is installed anc the `cmake --version` should return the version information.
 
 ## MPI
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
  
## Boost
 The next piece ichnos needs is boost. I'll install 1.69 as this is the one already exists from previous projects. 
 ```
 spack install boost@1.69.0+mpi ^openmpi@3.1.4 target=x86_64
 ```
 ## CGAL
 CGAL is the next library that is required by ichnos. CGAL depends on boost so we will specify that we want to build cgal using the boost version we just installed
 ```
 spack install cgal@4.13 ^boost@1.69.0+mpi ^openmpi@3.1.4 target=x86_64
 ```
 Sometimes cgal failed to get build because the mprf failed to apply a patch. One way to bypass this is by using a different version of mpfr
 ```
 spack install cgal@4.13 ^boost@1.69.0+mpi ^openmpi@3.1.4 ^mpfr@4.0.0
 ```
 
 ## Nanoflann
 The last required library is the nanoflann. THis is just a header that has no dependancies so all is needed is
 ```
 spack install nanoflann target=x86_64
 ```
 # Building
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

The following is not so relevant to be deleted soon...

 # Building in aqua cluster
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



