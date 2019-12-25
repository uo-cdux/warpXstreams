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

add_executable(advection advection.cxx)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
target_link_libraries(advection PRIVATE vtkm_filter ${Boost_LIBRARIES})
vtkm_add_target_information(advection
                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
                            DEVICE_SOURCES advection.cxx)
if(TARGET vtkm::tbb)
  target_compile_definitions(advection PRIVATE BUILDING_TBB_VERSION)
endif()
