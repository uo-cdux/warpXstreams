#include <stdio.h>

#include <iostream>
#include <string>

#include "boost/program_options.hpp"

#include <vtkm/Types.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/WorkletMapField.h>

#include <vtkm/worklet/particleadvection/Field.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
#include <vtkm/worklet/particleadvection/EulerIntegrator.h>
#include <vtkm/worklet/particleadvection/ParticleAdvectionWorklets.h>

#include "SeedingConfig.h"
#include "SeedGenerator.hxx"
#include "ValidateOptions.hxx"

void GenerateRandomIndices(std::vector<vtkm::Id>& randoms, vtkm::Id numberOfSeeds, vtkm::Id total)
{
  srand(314);
  for (int index = 0; index < numberOfSeeds; index++)
  {
    randoms.push_back(rand() % total);
  }
}

int main(int argc, char **argv) {
  vtkm::cont::SetStderrLogLevel(vtkm::cont::LogLevel::Off);

  namespace options = boost::program_options;
  options::options_description desc("Options");
  desc.add_options()("data", options::value<std::string>()->required(), "Path to dataset")
                    ("field", options::value<std::string>()->required(), "Name of vector field")
                    ("steps", options::value<vtkm::Id>()->required(), "Number of Steps")
                    ("length", options::value<vtkm::FloatDefault>()->required(), "Length of a single step")
                    ("seeds", options::value<vtkm::Id>()->required(), "Number of Seeds")
                    ("file", options::value<std::string>()->required(), "VTK file to read electrons from");
  options::variables_map vm;
  options::store(options::parse_command_line(argc, argv, desc), vm); // can throw
  options::notify(vm);

  //seeding::SeedingConfig config;
  if (!(vm.count("data")
      && vm.count("steps")
      && vm.count("length")
      && vm.count("field")
      && vm.count("seeds")
      && vm.count("file")))
  {
    std::cout << "Advection Benchmark" << std::endl << desc << std::endl;
  }

  std::string data = vm["data"].as<std::string>();
  std::string fieldname = vm["field"].as<std::string>();
  vtkm::Id steps = vm["steps"].as<vtkm::Id>();
  vtkm::FloatDefault length = vm["length"].as<vtkm::FloatDefault>();
  vtkm::Id numSeeds = vm["seeds"].as<vtkm::Id>();
  std::string seeddata = vm["file"].as<std::string>();

  using ArrayType = vtkm::cont::ArrayHandle<vtkm::Vec3f>;
  using FieldType = vtkm::worklet::particleadvection::ElectroMagneticField<ArrayType>;
  using IndexType = vtkm::cont::ArrayHandle<vtkm::Id>;
  using SeedsType = vtkm::cont::ArrayHandle<vtkm::Electron>;
  using EvaluatorType = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
  using IntegratorType = vtkm::worklet::particleadvection::EulerIntegrator<EvaluatorType>;
  using ParticleType = vtkm::worklet::particleadvection::Particles<vtkm::Electron>;
  using AdvectionWorklet = vtkm::worklet::particleadvection::ParticleAdvectWorklet;

  vtkm::io::VTKDataSetReader dataReader(data);
  vtkm::cont::DataSet dataset = dataReader.ReadDataSet();
  vtkm::cont::DynamicCellSet cells = dataset.GetCellSet();
  vtkm::cont::CoordinateSystem coords = dataset.GetCoordinateSystem();

  auto bounds = coords.GetBounds();
  std::cout << "Bounds : " << bounds << std::endl;

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

  {
    SeedsType allSeeds;
    vtkm::io::VTKDataSetReader seedsReader(seeddata);
    vtkm::cont::DataSet seedsData = seedsReader.ReadDataSet();
    seeding::GenerateElectrons(seedsData, allSeeds);

    std::vector<vtkm::Id> randoms;
    GenerateRandomIndices(randoms, numSeeds, allSeeds.GetNumberOfValues());
    IndexType toKeep = vtkm::cont::make_ArrayHandle(randoms);

    std::cout << "To keep : " << toKeep.GetNumberOfValues() <<  std::endl;

    vtkm::cont::ArrayHandlePermutation<IndexType, SeedsType> temp(toKeep, allSeeds);
    vtkm::cont::Algorithm::Copy(temp, seeds);
    auto portal = seeds.ReadPortal();
    for (vtkm::Id i = 0; i < portal.GetNumberOfValues(); i++)
    {
      auto electron = portal.Get(i);
      std::cout << electron.Pos 
                << ", m : " << electron.Mass 
                << ", u : " << electron.Momentum 
                << std::endl;
    }
    std::cout << "Permutation : " << seeds.GetNumberOfValues() <<  std::endl;
  }

  std::cout << "Advecting " << seeds.GetNumberOfValues() << " particles" << std::endl;

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

  std::cout << "Advection : " << timer.GetElapsedTime() << std::endl;

  return 1;
}
