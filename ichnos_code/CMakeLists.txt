cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

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

find_package( CGAL QUIET COMPONENTS  )
include( ${CGAL_USE_FILE} )

find_package(Boost 1.69 REQUIRED COMPONENTS mpi serialization program_options system)
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH} ${Boost_INCLUDE_DIRS})


add_executable(${PROJECT_NAME} "ichnos.cpp")

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Boost_LIBRARIES} ${CGAL_LIBRARIES} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
