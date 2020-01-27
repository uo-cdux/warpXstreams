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

#include <papi.h>

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

  int retval, num_hwcntrs = 0;
  int Events[NUM_EVENTS] = {PAPI_FP_OPS, PAPI_SP_OPS, PAPI_L3_TCA, PAPI_L3_TCM};
  char errstring[PAPI_MAX_STR_LEN];
  /*This is going to store our list of results*/
  long long values[NUM_EVENTS];

  /***************************************************************************
  *  This part initializes the library and compares the version number of the*
  * header file, to the version of the library, if these don't match then it *
  * is likely that PAPI won't work correctly.If there is an error, retval    *
  * keeps track of the version number.                                       *
  ***************************************************************************/
  if((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT )
  {
     fprintf(stderr, "Error: %d %s\n",retval, errstring);
     exit(1);
  }
  /**************************************************************************
   * PAPI_num_counters returns the number of hardware counters the platform *
   * has or a negative number if there is an error                          *
   **************************************************************************/
  if ((num_hwcntrs = PAPI_num_counters()) < PAPI_OK)
  {
     printf("There are no counters available. \n");
     exit(1);
  }
  printf("There are %d counters in this system\n",num_hwcntrs);


  timer.Start();

  if ( (retval = PAPI_start_counters(Events, NUM_EVENTS)) != PAPI_OK)
     ERROR_RETURN(retval);

  printf("\nCounter Started: \n");

  vtkm::cont::Invoker invoker;
  invoker(AdvectionWorklet{}, indices, integrator, particles, particleSteps);

  if ( (retval=PAPI_read_counters(values, NUM_EVENTS)) != PAPI_OK)
     ERROR_RETURN(retval);

  printf("Read successfully\n");
  printf("The floating point operations : %lld \n",values[0]);
  printf("The single precision floating point operations : %lld \n", values[1] );
  printf("Last level data cache accesses : %lld \n", values[2] );
  printf("Last level data cache misses : %lld \n", values[3] );

  if ((retval=PAPI_stop_counters(values, NUM_EVENTS)) != PAPI_OK)
    ERROR_RETURN(retval);

  timer.Stop();

  std::cout << "Advection : " << timer.GetElapsedTime() << std::endl;
  VerifySeeds(seeds);

  return 1;
}
