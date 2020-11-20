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

add_executable(advection advection.cxx Config.h SeedGenerator.hxx ValidateOptions.hxx)
add_executable(electrons MakeElectronDataset.cxx)
target_link_libraries(advection PRIVATE vtkm_cont vtkm_io vtkm_worklet vtkm_filter ${Boost_LIBRARIES})
target_link_libraries(electrons PRIVATE vtkm_cont vtkm_io vtkm_worklet vtkm_filter)
vtkm_add_target_information(advection
                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
                            DEVICE_SOURCES advection.cxx)
if(TARGET vtkm::tbb)
  target_compile_definitions(advection PRIVATE BUILDING_TBB_VERSION)
endif()

#add_executable(simpledata simpledata.cxx)
#target_link_libraries(simpledata PRIVATE vtkm_cont ${Boost_LIBRARIES})
#vtkm_add_target_information(simpledata
#                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
#                            DEVICE_SOURCES simple.cxx)
