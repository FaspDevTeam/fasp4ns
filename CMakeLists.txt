# CMakeLists for FASP4NS

###############################################################################
## General environment setting
###############################################################################

# Minimum cmake version needed
cmake_minimum_required(VERSION 2.8.12)

# Name of the project
project(fasp4ns C CXX)

# Set module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/modules)

###############################################################################
## Set configuration types and their default compiler flags
###############################################################################

set(CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE STRING "Configs" FORCE)
if (MSVC)
    # Compiler flags for Windows MSVC, skip Fortran
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O3")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /W4")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O3")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W4")
else()
    # Compiler flags for Linux or Mac OS X
    set(ADD_STDLIBS "${ADD_STDLIBS} m")
    string(STRIP ${ADD_STDLIBS} ADD_STDLIBS)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
endif()

set(build_type "Release")
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "${build_type}" CACHE STRING "the build type" FORCE)
endif()
message(STATUS "Set build type to ${CMAKE_BUILD_TYPE} configuration")

###############################################################################
# OpenMP : defined on command line in the top Makefile
# For cmake since 3.9, we can use find_package(OpenMP REQUIRED)
# We keep this for older cmake versions!!!
###############################################################################

if(USE_OPENMP)
    if(CMAKE_VERSION VERSION_GREATER 3)
        cmake_policy(SET "CMP0054" NEW)
    endif()

    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")	 
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
    else(OPENMP_FOUND)
        message("   WARNING: OpenMP was requested but disabled!")
    endif(OPENMP_FOUND)
endif(USE_OPENMP)

###############################################################################
## Project specific parameters
###############################################################################

# Find source files required for PennSim library
file(GLOB S1 ${PROJECT_SOURCE_DIR}/src/*.c)
add_library(fasp4ns STATIC ${S1})

# Add include directory
include_directories(${PROJECT_SOURCE_DIR}/include/)

# Find FASP, requires faspsolver and fasp4blkoil     
find_package(FASP REQUIRED)

if (FASP_FOUND)
    include_directories(${FASP_INCLUDE_DIRS})
    include_directories(${FASP_SRC_DIRS})
    set(CMAKE_REQUIRED_INCLUDES ${FASP_INCLUDE_DIRS} ${FASP_DIR}/base/src)
else(FASP_FOUND)
    message("   ERROR: FASP was requested but not supported!")
endif(FASP_FOUND)

###############################################################################
## External libraries
###############################################################################

include(external)

###############################################################################
## Targets for examples and installation
###############################################################################

# Set install location for fasp4ns library
set(CMAKE_INSTALL_LIBDIR "${PROJECT_SOURCE_DIR}/lib/")
install(TARGETS fasp4ns
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Compile and link NavierStokes.ex
add_executable(NavierStokes.ex main/NavierStokes.c)
target_link_libraries(NavierStokes.ex fasp4ns 
                      ${FASP_LIBRARIES}
                      ${BLAS_LIBRARIES} 
                      ${UMFPACK_LIBRARIES} 
                      ${SUPERLU_LIBRARIES} 
                      ${MUMPS_LIBRARIES}
                      ${MKL_LIBRARIES}
                      ${ADD_STDLIBS})

install(TARGETS NavierStokes.ex 
        DESTINATION ${PROJECT_SOURCE_DIR})

# Auto generating headers if needed
add_custom_target(headers ${PROJECT_SOURCE_DIR}/util/fhead.sh "${PROJECT_SOURCE_DIR}"
                  DEPENDS ${PROJECT_SOURCE_DIR}
                  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/util"
                  COMMENT "Auto-generating header file for C functions..."
                  VERBATIM)

# CMake uninstall target
if(NOT TARGET uninstall)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/modules/cmake_uninstall.cmake.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
                   IMMEDIATE @ONLY)
    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
