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
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
#include <vtkm/worklet/particleadvection/Integrators.h>
#include <vtkm/worklet/particleadvection/ParticleAdvectionWorklets.h>

#include "SeedingConfig.h"
#include "SeedGenerator.hxx"
#include "ValidateOptions.hxx"

#include <papi.h>
#include <likwid.h>
//extern "C" {
//#include <variorum.h>
//}
//#include "advisor-annotate.h"

#define NUM_EVENT 3
#define THRESHOLD 10000
#define ERROR_RETURN(retval) { fprintf(stderr, "Error %d %s:line %d: \n", retval,__FILE__,__LINE__);  exit(retval); }

int main(int argc, char **argv) {
  vtkm::cont::SetStderrLogLevel(vtkm::cont::LogLevel::Off);

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
  using ParticleType = vtkm::worklet::particleadvection::StateRecordingParticles;
  //using AdvectionWorklet = vtkm::worklet::particleadvection::ParticleAdvectWorklet;
  using AdvectionWorklet = vtkm::worklet::StreamlineWorklet;

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
  //std::cout << "Pre-requisite : " << timer.GetElapsedTime() << std::endl;
  timer.Reset();

  //LIKWID_MARKER_INIT;
  //LIKWID_MARKER_START("advection");
  //ANNOTATE_SITE_BEGIN(advect);
  //FILE* fp;
  //fp = fopen("prof.txt", "w+");
  //variorum_monitoring(stdout);

  int retval;
  int EventSet = PAPI_NULL;
  /*must be initialized to PAPI_NULL before calling PAPI_create_event*/
  int event_codes[NUM_EVENT]=  {PAPI_FP_OPS, PAPI_L3_TCM, PAPI_SR_INS};
  long long values[NUM_EVENT];
  if((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT )
    ERROR_RETURN(retval);
  /* Creating event set   */
  if ((retval=PAPI_create_eventset(&EventSet)) != PAPI_OK)
  ERROR_RETURN(retval);
  /* Add the array of events PAPI_TOT_INS and PAPI_TOT_CYC to the eventset*/
  if ((retval=PAPI_add_events(EventSet, event_codes, NUM_EVENT)) != PAPI_OK)
    ERROR_RETURN(retval);
  /* Start counting */
  if ( (retval=PAPI_start(EventSet)) != PAPI_OK)
    ERROR_RETURN(retval);

  timer.Start();
  //vtkm::cont::Invoker invoker;
  //invoker(AdvectionWorklet{}, indices, integrator, particles, particleSteps);
  AdvectionWorklet worklet;
  vtkm::worklet::StreamlineResult result;
  result = worklet.Run(integrator, particles, steps);
  timer.Stop();

  vtkm::cont::DataSet output;
  vtkm::cont::CoordinateSystem outputCoords("coordinates", res.Positions);
  output.SetCellSet(res.PolyLines);
  output.AddCoordinateSystem(outputCoords);
  vtkm::io::writer::VTKDataSetWriter writer("output.vtk");
  writer.WriteDataSet(output);

  /* Stop counting, this reads from the counter as well as stop it. */
  if ( (retval=PAPI_stop(EventSet,values)) != PAPI_OK)
    ERROR_RETURN(retval);

  long long flop = values[0];
  long long memory = (values[1] + values[2]) * 8;
  double time = timer.GetElapsedTime();
  // Convert to GFlop/s
  double flop_per_sec = static_cast<double>(flop) / (1000*1000*1000) / time;
  double flop_per_byte = static_cast<double>(flop) / static_cast<double>(memory);

  //variorum_monitoring(stdout);
  //ANNOTATE_SITE_END(advect);
  //LIKWID_MARKER_STOP("advection");
  //LIKWID_MARKER_CLOSE;

  std::cout << "MFLOP/sec : " << flop_per_sec << ", FLOP/Byte : " << flop_per_byte << std::endl;

  //VerifySeeds(seeds);
  //fclose(fp);
  return 1;
}
