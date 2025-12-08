# (Required) Minimum CMake version
cmake_minimum_required(VERSION 3.5)

# (Optional) Project/program/software name
project(phipsi)

# (Required) Enable Fortran support
enable_language(Fortran)

# Enable C support
enable_language(C)

# (Optional) Set build type to Release, which enables -O3 optimization
set(CMAKE_BUILD_TYPE Release)

# Set an option named ENABLE_OPENMP with default value ON
option(ENABLE_OPENMP "Enable OpenMP" ON)
# Enable OpenMP if ENABLE_OPENMP is ON
if(ENABLE_OPENMP)
    find_package(OpenMP) 
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}") 
endif()

# include_directories(${CMAKE_SOURCE_DIR}/inc)
link_directories(${CMAKE_SOURCE_DIR}/lib)

# Add all Fortran source files under src folder (including subfolders)
file(GLOB_RECURSE SRCS RELATIVE ${PROJECT_SOURCE_DIR} src/*)

# Static compilation.
# set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# set(BUILD_SHARED_LIBS OFF)
# set(CMAKE_EXE_LINKER_FLAGS "-static")
# Static compilation.
# set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# set(BUILD_SHARED_LIBS OFF)
# set(CMAKE_EXE_LINKER_FLAGS "-static")

# Only add -z,noexecstack on Linux (or non-Apple systems)
if(UNIX AND NOT APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-z,noexecstack")
endif()
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod")
add_executable(phipsi ${SRCS})

# Add conditional compilation macros. 
add_definitions(-Dgfortran -Dlinux -Dgithub)

# Add include file
target_include_directories(phipsi PUBLIC "include/")

# Add compiled libraries
target_link_libraries(phipsi libsmumps.a)
target_link_libraries(phipsi libdmumps.a)
target_link_libraries(phipsi libmumps_common.a)
target_link_libraries(phipsi libpord.a)
target_link_libraries(phipsi liblapack.a)
target_link_libraries(phipsi librefblas.a)
target_link_libraries(phipsi libmpiseq.a)

# Add libraries installed in Linux system
# target_link_libraries(phipsi -lsuperlu)