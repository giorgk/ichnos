cmake_minimum_required(VERSION 3.15)

set(CMAKE_BUILD_TYPE Debug)
if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
endif(UNIX)
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(ichnos)

#OPTION (USE_OpenMP "Use OpenMP" ON)
#IF(USE_OpenMP)
#    FIND_PACKAGE(OpenMP)
#    IF(OPENMP_FOUND)
#        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    ENDIF()
#ENDIF()
find_package(MPI REQUIRED)

if(USEHF)
  message(STATUS "Compile with HDF5")
  find_package(HighFive REQUIRED)
  find_package(HDF5 REQUIRED)
  add_definitions(-D_USEHF=1)
else()
  message(STATUS "Compile without HDF5")
  add_definitions(-D_USEHF=0)
endif(USEHF)

find_package( CGAL QUIET COMPONENTS  )
include( ${CGAL_USE_FILE} )

find_package(Boost 1.58 REQUIRED COMPONENTS mpi serialization program_options system)
if(USEHF)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH} ${Boost_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
else()
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH} ${Boost_INCLUDE_DIRS})
endif(USEHF)

#INCLUDE_DIRECTORIES(${GRID_INTERP})

add_executable(${PROJECT_NAME} "ichnos.cpp")

if(USEHF)
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Boost_LIBRARIES} ${CGAL_LIBRARIES} ${MPI_LIBRARIES} ${HDF5_C_LIBRARIES})
else()
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Boost_LIBRARIES} ${CGAL_LIBRARIES} ${MPI_LIBRARIES})
endif(USEHF)

if(MPI_COMPILE_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
