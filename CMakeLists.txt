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

set(VTKm_DIR "/home/local/KHQ/abhi.yenpure/repositories/mine/build/lib/cmake/vtkm-1.8/")
set(VTK_DIR "/home/local/KHQ/abhi.yenpure/repositories/vtk/build/lib/cmake/vtk-9.2")
set(CMAKE_CUDA_ARCHITECTURES 75)

# Find the VTK-m package
find_package(VTK REQUIRED QUIET)
find_package(VTKm REQUIRED QUIET)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

add_executable(savedata savedata.cxx Config.h SeedGenerator.hxx ValidateOptions.hxx FilterStreamlines.h)
target_link_libraries(savedata PRIVATE vtkm_cont vtkm_io vtkm_worklet vtkm_filter ${Boost_LIBRARIES} ${VTK_LIBRARIES})
