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

#Find the VTK-m package
find_package(VTKm REQUIRED QUIET)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

add_executable(advection advection.cxx)
target_link_libraries(advection PRIVATE vtkm_cont vtkm_filter ${Boost_LIBRARIES})
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

add_executable(roofexperiments roofexperiments.cxx)
target_link_libraries(roofexperiments PRIVATE vtkm_cont vtkm_filter ${Boost_LIBRARIES})
vtkm_add_target_information(roofexperiments
                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
                            DEVICE_SOURCES roofexperiments.cxx)
if(TARGET vtkm::tbb)
  target_compile_definitions(roofexperiments PRIVATE BUILDING_TBB_VERSION)
endif()
