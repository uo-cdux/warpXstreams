#include <stdio.h>

#include <iostream>
#include <string>

#include "boost/program_options.hpp"

#include <vtkm/Types.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/io/reader/VTKDataSetReader.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/WorkletMapField.h>

#include <vtkm/worklet/particleadvection/Field.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
#include <vtkm/worklet/particleadvection/Integrators.h>
#include <vtkm/worklet/particleadvection/ParticleAdvectionWorklets.h>

#include "SeedingConfig.h"
#include "SeedGenerator.hxx"
#include "ValidateOptions.hxx"

int main(int argc, char **argv) {
  vtkm::cont::SetStderrLogLevel(vtkm::cont::LogLevel::Off);

  namespace options = boost::program_options;
  options::options_description desc("Options");
  desc.add_options()("data", options::value<std::string>()->required(), "Path to dataset")
                    ("field", options::value<std::string>()->required(), "Name of vector field")
                    ("steps", options::value<vtkm::Id>()->required(), "Number of Steps")
                    ("length", options::value<vtkm::FloatDefault>()->required(), "Length of a single step")
                    ("seeds", options::value<std::string>()->required(), "VTK file to read electrons from");
  options::variables_map vm;
  options::store(options::parse_command_line(argc, argv, desc), vm); // can throw
  options::notify(vm);

  //seeding::SeedingConfig config;
  if (!(vm.count("data")
      && vm.count("steps")
      && vm.count("length")
      && vm.count("field")
      && vm.count("seeds")))
  {
    std::cout << "Advection Benchmark" << std::endl << desc << std::endl;
  }

  std::string data = vm["data"].as<std::string>();
  std::string fieldname = vm["field"].as<std::string>();
  vtkm::Id steps = vm["steps"].as<vtkm::Id>();
  vtkm::FloatDefault length = vm["length"].as<vtkm::FloatDefault>();
  std::string seeddata = vm["seeds"].as<std::string>();

  using ArrayType = vtkm::cont::ArrayHandle<vtkm::Vec3f>;
  using FieldType = vtkm::worklet::particleadvection::ElectroMagneticField<ArrayType>;
  using SeedsType = vtkm::cont::ArrayHandle<vtkm::Electron>;
  using EvaluatorType = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
  using IntegratorType = vtkm::worklet::particleadvection::RK4Integrator<EvaluatorType>;
  using ParticleType = vtkm::worklet::particleadvection::Particles<vtkm::Electron>;
  using AdvectionWorklet = vtkm::worklet::particleadvection::ParticleAdvectWorklet;

  vtkm::io::reader::VTKDataSetReader dataReader(data);
  vtkm::cont::DataSet dataset = dataReader.ReadDataSet();
  vtkm::cont::DynamicCellSet cells = dataset.GetCellSet();
  vtkm::cont::CoordinateSystem coords = dataset.GetCoordinateSystem();

  vtkm::cont::Timer timer;

  timer.Start();

  ArrayType electric, magnetic;
  dataset.GetField("E").GetData().CopyTo(electric);
  dataset.GetField("B").GetData().CopyTo(magnetic);
  FieldType electromagnetic(electric, magnetic);

  EvaluatorType evaluator(coords, cells, electromagnetic);
  IntegratorType integrator(evaluator, length);

  /*
   * Make seeds based on the seeding option.
   */
  SeedsType seeds;
  vtkm::io::reader::VTKDataSetReader seedsReader(seeddata);
  vtkm::cont::DataSet seedsData = seedsReader.ReadDataSet();
  //seeding::GenerateSeeds(config, dataset, seeds);
  seeding::GenerateElectrons(seedsData, seeds);

  vtkm::cont::ArrayHandleConstant<vtkm::Id> particleSteps(
      steps, seeds.GetNumberOfValues());
  ParticleType particles(seeds, steps);
  vtkm::cont::ArrayHandleIndex indices(seeds.GetNumberOfValues());

  timer.Stop();
  std::cout << "Pre-requisite : " << timer.GetElapsedTime() << std::endl;
  timer.Reset();

  timer.Start();
  vtkm::cont::Invoker invoker;
  invoker(AdvectionWorklet{}, indices, integrator, particles, particleSteps);
  timer.Stop();

  return 1;
}
