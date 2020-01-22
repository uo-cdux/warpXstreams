#include <iostream>
#include <string>

#include "boost/program_options.hpp"

#include <vtkm/Types.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/io/reader/VTKDataSetReader.h>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
#include <vtkm/worklet/particleadvection/Integrators.h>
#include <vtkm/worklet/particleadvection/ParticleAdvectionWorklets.h>

#include "SeedingOptions.hxx"

int main(int argc, char **argv) {
  namespace options = boost::program_options;
  options::options_description desc("Options");
  desc.add_options()("data", options::value<std::string>()->required(), "Path to dataset")
                    ("field", options::value<std::string>()->required(), "Name of vector field")
                    ("steps", options::value<vtkm::Id>()->required(), "Number of Steps")
                    ("length", options::value<vtkm::FloatDefault>()->required(), "Length of a single step")
                    ("seeding", options::value<int>()->required(), "Seeding options : UNIFORM, RANDOM");
                    ("density", options::value<std::vector<vtkm::Id>>()->multitoken(), "Reduction factor in density");
                    ("point", options::value<std::vector<vtkm::Id>>()->multitoken(), "Reduction factor in density");
                    ("seeds", options::value<std::vector<vtkm::Id>>()->multitoken(), "Number of Seeds");
  options::variables_map vm;
  options::store(options::parse_command_line(argc, argv, desc), vm); // can throw
  options::notify(vm);

  seeding::SeedingConfig config;
  if (!(vm.count("data")
      && vm.count("steps")
      && vm.count("length")
      && vm.count("field"))
      && seeding::ValidateSeedingOptions(vm, config))
  {
    std::cout << "Advection Benchmark" << std::endl << desc << std::endl;
  }
  // Assuming uniform seeding for now.

  std::string data = vm["data"].as<std::string>();
  std::string fieldname = vm["field"].as<std::string>();
  vtkm::Id steps = vm["steps"].as<vtkm::Id>();
  vtkm::FloatDefault length = vm["length"].as<vtkm::FloatDefault>();

  using FieldType = vtkm::cont::ArrayHandle<vtkm::Vec3f>;
  using SeedsType = vtkm::cont::ArrayHandle<vtkm::Particle>;
  using EvaluatorType =
      vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
  using IntegratorType =
      vtkm::worklet::particleadvection::RK4Integrator<EvaluatorType>;
  using ParticleType = vtkm::worklet::particleadvection::Particles;

  vtkm::io::reader::VTKDataSetReader dataReader(data);
  vtkm::cont::DataSet dataset = dataReader.ReadDataSet();
  vtkm::cont::DynamicCellSet cells = dataset.GetCellSet();
  vtkm::cont::CoordinateSystem coords = dataset.GetCoordinateSystem();

  vtkm::cont::Timer timer;

  timer.Start();

  FieldType field;
  dataset.GetField(fieldname).GetData().CopyTo(field);

  EvaluatorType evaluator(coords, cells, field);
  IntegratorType integrator(evaluator, length);

  /*
   * Make seeds based on the seeding option.
   */
  SeedsType seeds;
  seeding::GenerateSeeds(config, dataset, seeds);
  vtkm::cont::ArrayHandleConstant<vtkm::Id> particleSteps(
      steps, seeds.GetNumberOfValues());
  ParticleType particles(seeds, steps);

  vtkm::cont::ArrayHandleIndex indices(seeds.GetNumberOfValues());

  timer.Stop();

  std::cout << "Pre-requisite : " << timer.GetElapsedTime() << std::endl;

  timer.Reset();
  timer.Start();

  using AdvectionWorklet =
      vtkm::worklet::particleadvection::ParticleAdvectWorklet;
  vtkm::cont::Invoker invoker;
  invoker(AdvectionWorklet{}, indices, integrator, particles, particleSteps);

  timer.Stop();

  std::cout << "Advection : " << timer.GetElapsedTime() << std::endl;

  return 1;
}
