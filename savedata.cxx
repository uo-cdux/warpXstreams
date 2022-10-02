#include <stdio.h>

#include <iostream>
#include <string>


#include <vtkm/Types.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandlePermutation.h>
#include <vtkm/cont/Timer.h>

#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/worklet/WorkletMapField.h>

#include <vtkm/filter/flow/Streamline.h>

#include "Config.h"
#include "SeedGenerator.hxx"
#include "ValidateOptions.hxx"

#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>

#include <stdio.h>

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

class ExtractParticleData : public vtkm::worklet::WorkletMapField
{
public:
  ExtractParticleData() {};
  using ControlSignature = void (FieldIn, FieldOut, FieldOut, FieldOut, FieldOut, FieldOut);

  VTKM_EXEC void operator()(const vtkm::ChargedParticle& particle,
                            vtkm::Vec3f& position,
                            vtkm::Vec3f& momentum,
                            vtkm::FloatDefault& mass,
                            vtkm::FloatDefault& charge,
                            vtkm::FloatDefault& weighting) const
  {
    position  = particle.Pos;
    momentum  = particle.Momentum;
    mass      = particle.Mass;
    charge    = particle.Charge;
    weighting = particle.Weighting;
  }
};

void ExtractDataSetFromSeeds(const vtkm::cont::ArrayHandle<vtkm::ChargedParticle>& seeds)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> pos, mom;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> mass, charge, weighting;
  vtkm::cont::Invoker invoker;
  invoker(ExtractParticleData{}, seeds, pos, mom, mass, charge, weighting);
  auto numvalues = seeds.GetNumberOfValues();

/*  vtkm::cont::Token token;
  const double* p_pos  = pos.GetBuffers()[0].ReadPointerHost(token);
  const double* p_mom  = mom.GetBuffers()[0].ReadPointerHost(token);
  const double* p_mass = mass.GetBuffers()[0].ReadPointerHost(token);
  const double* p_c    = charge.GetBuffers()[0].ReadPointerHost(token);
  const double* p_w    = weighting.GetBuffers()[0].ReadPointerHost(token);
*/
  vtkNew<vtkPoints> points;
  points->SetDataTypeToDouble();
  points->SetNumberOfPoints(numvalues);
  auto p_pos = pos.ReadPortal();
  for(int i = 0; i < numvalues; i++)
  {
    auto _point = p_pos.Get(i);
    double point[3] = {_point[0], _point[1], _point[2]};
    points->SetPoint(i, point);
  }
  vtkNew<vtkPolyData> polyData;
  polyData->SetPoints(points);

  vtkNew<vtkDoubleArray> momentum;
  momentum->SetNumberOfComponents(3);
  momentum->SetName("Momentum");
  // Create the data to store (here we just use (0,0,0))
  auto p_mom = mom.ReadPortal();
  for(int i = 0; i < numvalues; i++)
  {
    auto _mom = p_mom.Get(i);
    double momx[3] = {_mom[0], _mom[1], _mom[2]};
    momentum->InsertNextTuple(momx);
  }
  // The data is added to FIELD data (rather than POINT data as usual)
  polyData->GetPointData()->AddArray(momentum);

  vtkNew<vtkDoubleArray> mass1;
  mass1->SetNumberOfComponents(1);
  mass1->SetName("Mass");
  // Create the data to store (here we just use (0,0,0))
  auto p_mass = mass.ReadPortal();
  for(int i = 0; i < numvalues; i++)
  {
    double _mass = p_mass.Get(i);
    mass1->InsertNextTuple1(_mass);
  }
  // The data is added to FIELD data (rather than POINT data as usual)
  polyData->GetPointData()->AddArray(mass1);

  vtkNew<vtkDoubleArray> char1;
  char1->SetNumberOfComponents(1);
  char1->SetName("Charge");
  // Create the data to store (here we just use (0,0,0))
  auto p_char = charge.ReadPortal();
  for(int i = 0; i < numvalues; i++)
  {
    double _char = p_char.Get(i);
    char1->InsertNextTuple1(_char);
  }
  // The data is added to FIELD data (rather than POINT data as usual)
  polyData->GetPointData()->AddArray(char1);

  vtkNew<vtkDoubleArray> we1;
  we1->SetNumberOfComponents(1);
  we1->SetName("Weighting");
  // Create the data to store (here we just use (0,0,0))
  auto p_we = weighting.ReadPortal();
  for(int i = 0; i < numvalues; i++)
  {
    double _we = p_we.Get(i);
    we1->InsertNextTuple1(_we);
  }
  // The data is added to FIELD data (rather than POINT data as usual)
  polyData->GetPointData()->AddArray(we1);

  vtkNew<vtkPolyDataWriter> writer;
  writer->SetFileName("output.vtk");
  writer->SetInputData(polyData);
  writer->SetFileTypeToBinary();
  writer->Write();
  std::cerr << writer->GetOutputString() << std::endl;
/*
  FILE* f_out;
  f_out = fopen("temp.vtk", "w");
  fprintf(f_out, "# vtk DataFile Version 2.0\n");
  fprintf(f_out, "WarpX particle data\n");
  fprintf(f_out, "BINARY\n");
  fprintf(f_out, "DATASET POLYDATA\n");
  fprintf(f_out, "POINTS %ld double\n", numvalues);
  fwrite(p_pos, sizeof(vtkm::Vec3f), numvalues, f_out);
  fprintf(f_out, "\n");
  fprintf(f_out, "POINT_DATA %ld\n", numvalues);
  fprintf(f_out, "VECTORS mom double\n");
  fwrite(p_mom, sizeof(vtkm::Vec3f), numvalues, f_out);
  fprintf(f_out, "\n");
  fprintf(f_out, "SCALARS mass double 1\n");
  fprintf(f_out, "LOOKUP_TABLE default\n");
  fwrite(p_mass, sizeof(vtkm::FloatDefault), numvalues, f_out);
  fprintf(f_out, "\n");
  fprintf(f_out, "SCALARS charge double 1\n");
  fprintf(f_out, "LOOKUP_TABLE default\n");
  fwrite(p_c, sizeof(vtkm::FloatDefault), numvalues, f_out);
  fprintf(f_out, "\n");
  fprintf(f_out, "SCALARS w double 1\n");
  fprintf(f_out, "LOOKUP_TABLE default\n");
  fwrite(p_w, sizeof(vtkm::FloatDefault), numvalues, f_out);
  fprintf(f_out, "\n");
  fclose(f_out);
*/
}

void PrintSeeds(const vtkm::cont::ArrayHandle<vtkm::ChargedParticle>& seeds)
{
  auto portal = seeds.ReadPortal();
  for(int i = 0; i < portal.GetNumberOfValues(); i++)
  {
    auto s = portal.Get(i);
    std::cout << s.Pos << " " << s.Momentum  << " " << s.Mass  << " " << s.Charge  << " " << s.Weighting << std::endl;
  }
}

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
  using SeedsType = vtkm::cont::ArrayHandle<vtkm::ChargedParticle>;
  using IndexType = vtkm::cont::ArrayHandle<vtkm::Id>;
/*
  using FieldType = vtkm::worklet::flow::ElectroMagneticField<ArrayType>;
  using IndexType = vtkm::cont::ArrayHandle<vtkm::Id>;
  using EvaluatorType = vtkm::worklet::flow::GridEvaluator<FieldType>;
  using IntegratorType = vtkm::worklet::flow::RK4Integrator<EvaluatorType>;
  using Stepper = vtkm::worklet::flow::Stepper<IntegratorType, EvaluatorType>;
  using ParticleType = vtkm::worklet::flow::StateRecordingParticles<vtkm::ChargedParticle>;
  using AdvectionWorklet = vtkm::worklet::flow::ParticleAdvectWorklet;
*/
  SeedsType seeds;
/*
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
*/
  //std::cout << "Original data" << std::endl;
  //detail::PrintSeeds(seeds);

  //detail::ExtractDataSetFromSeeds(seeds);

  {
    SeedsType allSeeds;
    vtkm::io::VTKDataSetReader seedsReader("output.vtk");
    vtkm::cont::DataSet seedsData = seedsReader.ReadDataSet();
    vtkm::cont::ArrayHandle<vtkm::Vec3f> pos, mom;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> mass, charge, w;

    seedsData.GetCoordinateSystem().GetData().AsArrayHandle(pos);
    seedsData.GetField("Momentum").GetData().AsArrayHandle(mom);
    seedsData.GetField("Mass").GetData().AsArrayHandle(mass);
    seedsData.GetField("Charge").GetData().AsArrayHandle(charge);
    seedsData.GetField("Weighting").GetData().AsArrayHandle(w);
    seeding::GenerateChargedParticles(pos, mom, mass, charge, w, allSeeds);
    vtkm::cont::Algorithm::Copy(allSeeds, seeds);
  }

  std::cout << "Reconstructed data" << std::endl;
  detail::PrintSeeds(seeds);

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

  vtkm::filter::flow::Streamline streamline;

  streamline.SetStepSize(length);
  streamline.SetNumberOfSteps(50);
  streamline.SetSeeds(seeds);
  streamline.SetVectorFieldType(vtkm::filter::flow::VectorFieldType::ELECTRO_MAGNETIC_FIELD_TYPE);
  streamline.SetEField("E");
  streamline.SetBField("B");
  auto output = streamline.Execute(dataset);

  vtkm::io::VTKDataSetWriter writer1("streams.vtk");
  writer1.WriteDataSet(output);

  /*vtkm::cont::ArrayHandle<vtkm::Id> initSteps;
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
*/
}
