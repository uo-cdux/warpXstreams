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

#include "SeedingConfig.h"
#include "SeedGenerator.hxx"
#include "ValidateOptions.hxx"

//#include <papi.h>
//#include <likwid.h>
#include <variorum.h>

//#include "advisor-annotate.h"

#define NUM_EVENTS 4
#define THRESHOLD 10000
#define ERROR_RETURN(retval) { fprintf(stderr, "Error %d %s:line %d: \n", retval,__FILE__,__LINE__);  exit(retval); } 

int main(int argc, char **argv) {
  namespace options = boost::program_options;
  options::options_description desc("Options");
  desc.add_options()("data", options::value<std::string>()->required(), "Path to dataset")
                    ("field", options::value<std::string>()->required(), "Name of vector field")
                    ("steps", options::value<vtkm::Id>()->required(), "Number of Steps")
                    ("length", options::value<vtkm::FloatDefault>()->required(), "Length of a single step")
                    ("seeding", options::value<int>()->required(), "Seeding options : UNIFORM, RANDOM")
                    ("density", options::value<std::vector<vtkm::Id>>()->multitoken(), "Reduction factor in density")
                    ("point", options::value<std::vector<vtkm::FloatDefault>>()->multitoken(), "Point for single seed")
                    ("seeds", options::value<vtkm::Id>()->multitoken(), "Number of Seeds");
  options::variables_map vm;
  options::store(options::parse_command_line(argc, argv, desc), vm); // can throw
  options::notify(vm);

  seeding::SeedingConfig config;
  if (!(vm.count("data")
      && vm.count("steps")
      && vm.count("length")
      && vm.count("field")
      && validate::ValidateSeedingOptions(vm, config)))
  {
    std::cout << "Advection Benchmark" << std::endl << desc << std::endl;
  }

  std::string data = vm["data"].as<std::string>();
  std::string fieldname = vm["field"].as<std::string>();
  vtkm::Id steps = vm["steps"].as<vtkm::Id>();
  vtkm::FloatDefault length = vm["length"].as<vtkm::FloatDefault>();

  using FieldType = vtkm::cont::ArrayHandle<vtkm::Vec3f>;
  using SeedsType = vtkm::cont::ArrayHandle<vtkm::Particle>;
  using EvaluatorType = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
  using IntegratorType = vtkm::worklet::particleadvection::RK4Integrator<EvaluatorType>;
  using ParticleType = vtkm::worklet::particleadvection::Particles;
  using AdvectionWorklet = vtkm::worklet::particleadvection::ParticleAdvectWorklet;

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

  //LIKWID_MARKER_INIT;
  //LIKWID_MARKER_START("advection");
  //ANNOTATE_SITE_BEGIN(advect);
  variorum_monitoring(stdout);

  vtkm::cont::Invoker invoker;
  invoker(AdvectionWorklet{}, indices, integrator, particles, particleSteps);

  //ANNOTATE_SITE_END(advect);
  //LIKWID_MARKER_STOP("advection");
  //LIKWID_MARKER_CLOSE;

  timer.Stop();
  std::cout << "Advection : " << timer.GetElapsedTime() << std::endl;
  //VerifySeeds(seeds);

  return 1;
}
