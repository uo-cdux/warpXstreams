##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##============================================================================
cmake_minimum_required(VERSION 3.8...3.15 FATAL_ERROR)
project(advection CXX)

list(APPEND CMAKE_MODULE_PATH "~/repositories/CMakeFindModules")

# Find the VTK-m package
find_package(VTKm REQUIRED QUIET)

# Boost
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

# PAPI
find_package(PAPI REQUIRED)
include_directories(${PAPI_INCLUDE_DIRS})
link_directories(${PAPI_LIBRARIES})

# LIKWID
find_package(LIKWID REQUIRED QUIET)
include_directories(${LIKWID_INCLUDE_DIRS})
link_directories(${LIKWID_LIBRARIES})

# Intel Advisor
# find_package(ADVISOR REQUIRED QUIET)
# include_directories(${ADVISOR_INCLUDE_DIR})
# link_directories(${ADVISOR_LIBRARIES})

# Variorum
find_package(VARIORUM REQUIRED QUIET)
include_directories(${VARIORUM_INCLUDE_DIRS})
link_directories(${VARIORUM_LIBRARIES})

add_executable(advection advection.cxx SeedingConfig.h SeedGenerator.hxx ValidateOptions.hxx)
#target_link_libraries(advection PRIVATE vtkm_cont vtkm_filter papi pfm ${Boost_LIBRARIES} ${ADVISOR_LIBRARIES} ${CMAKE_DL_LIBS})
target_link_libraries(advection PRIVATE vtkm_cont vtkm_filter ${PAPI_LIBRARIES} ${Boost_LIBRARIES})
#target_link_libraries(advection PRIVATE vtkm_cont vtkm_filter ${VARIORUM_LIBRARIES} ${Boost_LIBRARIES})
#target_link_libraries(advection PRIVATE vtkm_cont vtkm_filter ${LIKWID_LIBRARIES} ${Boost_LIBRARIES})
vtkm_add_target_information(advection
                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
                            DEVICE_SOURCES advection.cxx)
if(TARGET vtkm::tbb)
  target_compile_definitions(advection PRIVATE BUILDING_TBB_VERSION)
endif()

add_executable(simpledata simpledata.cxx)
target_link_libraries(simpledata PRIVATE vtkm_cont ${Boost_LIBRARIES})
vtkm_add_target_information(simpledata
                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
                            DEVICE_SOURCES simple.cxx)

#add_executable(roofexperiments roofexperiments.cxx SeedingOptions.hxx)
#target_link_libraries(roofexperiments PRIVATE vtkm_cont vtkm_filter ${Boost_LIBRARIES})
#vtkm_add_target_information(roofexperiments
#                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
#                            DEVICE_SOURCES roofexperiments.cxx)
#if(TARGET vtkm::tbb)
#  target_compile_definitions(roofexperiments PRIVATE BUILDING_TBB_VERSION)
#endif()
