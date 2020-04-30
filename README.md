# ichnos

# Dependencies using Spack
First create an environment that we nname it ichnos
```
spack env create ichnos
```
and make it active
```
spack env activate ichnos
```
First we are going to install cmake. Any version should work however I'm specifing here 3.9.4 because this has already been built in other spack environments
 ```
 spack install cmake@3.9.4 target=x86_64
 ```
 Now the command `spack find` should find that the cmake is installed anc the `cmake --version` should return the version information.
 
 Next we will install openmpi. The `spack info openmpi` will return the available openmpi versions, while at the top specifies which one is the preffered. At the time of writing this that was the 3.1.4 so we build it as
 ```
 spack install openmpi@3.1.4 target=x86_64
 ```
 The next piece ichnos needs is boost. I'll install 1.69 as this is the one already exists from previous projects. 
 ```
 spack install boost@1.69.0+mpi ^openmpi@3.1.4 target=x86_64
 ```
 CGAL is the next library that is required by ichnos. CGAL depends on boost so we will specify that we want to build cgal using the boost version we just installed
 ```
 spack install cgal@4.13 ^boost@1.69.0+mpi ^openmpi@3.1.4 target=x86_64
 ```
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



