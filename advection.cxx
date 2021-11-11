#include <stdio.h>

#include <iostream>
#include <string>

#include "boost/program_options.hpp"

#include <vtkm/Types.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandlePermutation.h>
#include <vtkm/cont/Timer.h>

#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/particleadvection/Field.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
#include <vtkm/worklet/particleadvection/EulerIntegrator.h>
#include <vtkm/worklet/particleadvection/RK4Integrator.h>
#include <vtkm/worklet/particleadvection/ParticleAdvectionWorklets.h>

#include "Config.h"
#include "SeedGenerator.hxx"
#include "ValidateOptions.hxx"

namespace detail
{
class GetSteps : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  GetSteps() {}
  using ControlSignature = void(FieldIn, FieldOut);
  using ExecutionSignature = void(_1, _2);
  VTKM_EXEC void operator()(const vtkm::ChargedParticle& p, vtkm::Id& numSteps) const
  {
    numSteps = p.NumSteps;
  }
};

class ComputeNumPoints : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  ComputeNumPoints() {}
  using ControlSignature = void(FieldIn, FieldIn, FieldOut);
  using ExecutionSignature = void(_1, _2, _3);

  // Offset is number of points in streamline.
  // 1 (inital point) + number of steps taken (p.NumSteps - initalNumSteps)
  VTKM_EXEC void operator()(const vtkm::ChargedParticle& p,
                            const vtkm::Id& initialNumSteps,
                            vtkm::Id& diff) const
  {
      diff = 1 + p.NumSteps - initialNumSteps;
  }
};

} // namespace detail

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
  desc.add_options()("data",    options::value<std::string>()->required(),        "Path to dataset")
                    ("steps",   options::value<vtkm::Id>()->required(),           "Number of Steps")
                    ("length",  options::value<vtkm::FloatDefault>()->required(), "Length of a single step")
                    ("seeds",   options::value<vtkm::Id>(),        "Number of seeds for random/single seeding")
                    ("seeddata",  options::value<std::string>()->required(), "VTK file to read electrons from")
                    ("threshold", options::value<vtkm::FloatDefault>()->required(), "Foltering threshold")
                    ("sampleX",  options::value<std::string>(), "Seed sampling range X")
                    ("sampleY",  options::value<std::string>(), "Seed sampling range Y")
                    ("sampleZ",  options::value<std::string>(), "Seed sampling range Z");

  options::variables_map vm;
  std::ifstream settings_file(std::string(argv[1]), std::ifstream::in);
  options::store(options::parse_config_file(settings_file, desc), vm);
  settings_file.close();
  options::notify(vm);

  config::Config config;
  int res = validate::ValidateOptions(vm, config);
  if(res < 0)
  {
    std::cout << "Advection Benchmark" << std::endl << desc << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string data = config.GetDataSetName();
  vtkm::Id steps = config.GetNumSteps();
  vtkm::FloatDefault length = config.GetStepLength();
  vtkm::Id numSeeds = config.GetNumSeeds();
  std::string seeddata = config.GetSeedData();
  vtkm::FloatDefault threshold = config.GetThreshold();

  using ArrayType = vtkm::cont::ArrayHandle<vtkm::Vec3f>;
  using FieldType = vtkm::worklet::particleadvection::ElectroMagneticField<ArrayType>;
  using IndexType = vtkm::cont::ArrayHandle<vtkm::Id>;
  using SeedsType = vtkm::cont::ArrayHandle<vtkm::ChargedParticle>;
  using EvaluatorType = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
  using IntegratorType = vtkm::worklet::particleadvection::RK4Integrator<EvaluatorType>;
  using Stepper = vtkm::worklet::particleadvection::Stepper<IntegratorType, EvaluatorType>;
  using ParticleType = vtkm::worklet::particleadvection::StateRecordingParticles<vtkm::ChargedParticle>;
  using AdvectionWorklet = vtkm::worklet::particleadvection::ParticleAdvectWorklet;

  vtkm::io::VTKDataSetReader dataReader(data);
  vtkm::cont::DataSet dataset = dataReader.ReadDataSet();
  vtkm::cont::DynamicCellSet cells = dataset.GetCellSet();
  vtkm::cont::CoordinateSystem coords = dataset.GetCoordinateSystem();

  auto bounds = coords.GetBounds();
  std::cout << "Bounds : " << bounds << std::endl;
  using Structured3DType = vtkm::cont::CellSetStructured<3>;
  Structured3DType castedCells = cells.Cast<Structured3DType>();
  auto dims = castedCells.GetSchedulingRange(vtkm::TopologyElementTagPoint());
  vtkm::Vec3f spacing = {bounds.X.Length() / (dims[0] - 1),
                         bounds.Y.Length() / (dims[1] - 1),
                         bounds.Z.Length() / (dims[2] - 1)};
  std::cout << spacing << std::endl;
  constexpr static vtkm::FloatDefault SPEED_OF_LIGHT =
    static_cast<vtkm::FloatDefault>(2.99792458e8);
  spacing = spacing * spacing;
  length = 1.0 / (SPEED_OF_LIGHT * vtkm::Sqrt(1./spacing[0] + 1./spacing[1] + 1./spacing[2]));
  std::cout << "CFL length : " << length << std::endl;

  vtkm::cont::Timer timer;
  timer.Start();

  ArrayType electric, magnetic;
  dataset.GetField("E").GetData().AsArrayHandle(electric);
  dataset.GetField("B").GetData().AsArrayHandle(magnetic);
  FieldType electromagnetic(electric, magnetic);

  EvaluatorType evaluator(coords, cells, electromagnetic);
  Stepper stepper(evaluator, length);

  /*
   * Make seeds based on the seeding option.
   */
  SeedsType seeds;

  {
    SeedsType allSeeds;
    vtkm::io::VTKDataSetReader seedsReader(seeddata);
    vtkm::cont::DataSet seedsData = seedsReader.ReadDataSet();
    vtkm::cont::ArrayHandle<vtkm::Id> filter;
    seeding::GenerateChargedParticles(config, seedsData, allSeeds, filter);
    SeedsType _allSeeds;
    vtkm::cont::Algorithm::CopyIf(allSeeds, filter, _allSeeds);

    auto count = vtkm::cont::Algorithm::Reduce(filter, static_cast<vtkm::Id>(0));
    std::cout << "Sampled " << count << " electrons" << std::endl;

    std::vector<vtkm::Id> randoms;
    GenerateRandomIndices(randoms, numSeeds, _allSeeds.GetNumberOfValues());
    IndexType toKeep = vtkm::cont::make_ArrayHandle(randoms, vtkm::CopyFlag::On);

    vtkm::cont::ArrayHandlePermutation<IndexType, SeedsType> temp(toKeep, _allSeeds);
    vtkm::cont::Algorithm::Copy(temp, seeds);
  }

  vtkm::cont::Invoker invoker;
  std::cout << "Advecting " << seeds.GetNumberOfValues() << " particles" << std::endl;

  vtkm::cont::ArrayHandle<vtkm::Id> initSteps;
  invoker(detail::GetSteps{}, seeds, initSteps);

  vtkm::cont::ArrayHandleConstant<vtkm::Id> particleSteps(steps, seeds.GetNumberOfValues());
  ParticleType particles(seeds, steps);
  vtkm::cont::ArrayHandleIndex indices(seeds.GetNumberOfValues());

  timer.Stop();
  std::cout << "Pre-requisite : " << timer.GetElapsedTime() << std::endl;
  timer.Reset();

  timer.Start();
  invoker(AdvectionWorklet{}, indices, stepper, particles, particleSteps);
  timer.Stop();

  std::cout << "Advection : " << timer.GetElapsedTime() << std::endl;

  // Has the count of points in a streamline
  vtkm::cont::ArrayHandle<vtkm::Id> numPoints;
  invoker(detail::ComputeNumPoints{}, seeds, initSteps, numPoints);
  // Has all points for the streamline
  vtkm::cont::ArrayHandle<vtkm::Vec3f> streams;
  particles.GetCompactedHistory(streams);
  // Calculate the indices for streamlines cell set
  vtkm::cont::ArrayHandle<vtkm::Id> cellIndex;
  vtkm::Id connectivityLen = vtkm::cont::Algorithm::ScanExclusive(numPoints, cellIndex);
  // Connectivity for the cells
  vtkm::cont::ArrayHandleCounting<vtkm::Id> connCount(0, 1, connectivityLen);
  vtkm::cont::ArrayHandle<vtkm::Id> connectivity;
  vtkm::cont::ArrayCopy(connCount, connectivity);
  vtkm::cont::ArrayHandle<vtkm::UInt8> cellTypes;
  auto polyLineShape =
    vtkm::cont::make_ArrayHandleConstant<vtkm::UInt8>(vtkm::CELL_SHAPE_POLY_LINE, numSeeds);
  vtkm::cont::ArrayCopy(polyLineShape, cellTypes);
  auto numIndices = vtkm::cont::make_ArrayHandleCast(numPoints, vtkm::IdComponent());
  auto offsets = vtkm::cont::ConvertNumComponentsToOffsets(numIndices);
  vtkm::cont::CellSetExplicit<> polylines;
  polylines.Fill(streams.GetNumberOfValues(), cellTypes, connectivity, offsets);

  vtkm::cont::DataSet output;
  output.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coords", streams));
  output.SetCellSet(polylines);

  vtkm::io::VTKDataSetWriter writer1("streams.vtk");
  writer1.WriteDataSet(output);

  return 1;
}
