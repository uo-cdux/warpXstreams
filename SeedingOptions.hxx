#ifndef seeding_options_hxx
#define seeding_options_hxx

#include <vector>

#include "boost/program_options.hpp"

#include <vtkm/Types.h>
#include <vtkm/Particle.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/Invoker.h>
#include <vtkm/worklet/WorkletMapField.h>

namespace seeding
{

namespace options = boost::program_options;

enum class Options
{
  UNIFORM = 0,
  UNIFORM_SPARSE,
  SINGLE,
  SINGLE_COPIES,
  RANDOM,
};

class SeedingConfig
{
public:
  SeedingConfig() = default;

  void SetOption(Options option) {this->Option = option;}
  Options GetOption() {return this->Option;}

  void SetPoint(vtkm::Vec3f point) {this->Point = point;}
  vtkm::Vec3f GetPoint() {return this->Point;}

  void SetDensity(vtkm::Id3 density) {this->Density = density;}
  vtkm::Id3 GetDensity() {return this->Density;}

  void SetSeedCount(vtkm::Id seedCount) {this->SeedCount = seedCount; }
  vtkm::Id GetSeedCount() {return this->SeedCount;}

private:
  Options Option;
  vtkm::Vec3f Point;
  vtkm::Id3 Density;
  vtkm::Id SeedCount;
};

bool ValidateSeedingOptions(options::variables_map& vm,
                            SeedingConfig& config)
{
  if(!vm.count("seeding"))
    return false;
  int param = vm["seeding"].as<int>();
  Options option = static_cast<Options>(param);
  std::cout << "Option : " << param << std::endl;
  config.SetOption(option);
  switch(option)
  {
    case Options::UNIFORM:
    {
      std::cout << "uniform seeding" << std::endl;
      return true;
    }
    break;

    case Options::UNIFORM_SPARSE:
    {
      std::cout << "uniform sparse seeding" << std::endl;
      std::vector<vtkm::Id> vec = vm["density"].as<std::vector<vtkm::Id>>();
      vtkm::Id3 density{vec[0], vec[1], vec[2]};
      config.SetDensity(density);
      return true;
    }
    break;

    case Options::SINGLE:
    // Needs a point
    {
      std::cout << "single seeding" << std::endl;
      std::vector<vtkm::FloatDefault> vec = vm["point"].as<std::vector<vtkm::FloatDefault>>();
      vtkm::Vec3f point{vec[0], vec[1], vec[2]};
      config.SetPoint(point);
      return true;
    }
    break;

    case Options::SINGLE_COPIES:
    // Needs a point
    {
      std::cout << "single seeding" << std::endl;
      std::vector<vtkm::FloatDefault> vec = vm["point"].as<std::vector<vtkm::FloatDefault>>();
      vtkm::Vec3f point{vec[0], vec[1], vec[2]};
      config.SetPoint(point);
      vtkm::Id seedCount = vm["seeds"].as<vtkm::Id>();
      config.SetSeedCount(seedCount);
      return true;
    }
    break;

    case Options::RANDOM:
    // Needs a random number generator
    // Needs number of seeds
    {
      vtkm::Id seedCount = vm["seeds"].as<vtkm::Id>();
      config.SetSeedCount(seedCount);
      return true;
    }
    break;

    default:
      return false;
    break;
  }
}

class SeedsFromCoordinates : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn, FieldOut);
  using ExecutionSignature = void(WorkIndex, _1, _2);

  template <typename PointType>
  VTKM_EXEC void operator()(const vtkm::Id &index,
                            const PointType &point,
                            vtkm::Particle &particle) const
  {
    particle.ID = index;
    particle.Pos = point;
  }
};

void MakeSeedsFromCoordinates(vtkm::cont::CoordinateSystem& coords,
                              vtkm::cont::ArrayHandle<vtkm::Particle>& seeds)
{
  vtkm::cont::Invoker invoker;
  invoker(SeedsFromCoordinates{}, coords.GetData(), seeds);
}

void MakeSingleSeed(vtkm::Id seedCount,
                    vtkm::Vec3f& point,
                    vtkm::cont::ArrayHandle<vtkm::Particle>& seeds) {
  // Not implemented yet.
}

void MakeRandomSeeds(vtkm::cont::ArrayHandle<vtkm::Particle>& seeds) {
  // Not implemented yet.
}

void GenerateSeeds(SeedingConfig& config,
                   vtkm::cont::DataSet& dataset,
                   vtkm::cont::ArrayHandle<vtkm::Particle>& seeds)
{
  Options option = config.GetOption();
  switch(option)
  {
    case Options::UNIFORM:
    {
      std::cout << "Generating Uniform Seeds" << std::endl;
      vtkm::cont::CoordinateSystem coords = dataset.GetCoordinateSystem();
      MakeSeedsFromCoordinates(coords, seeds);
    }
    break;

    case Options::UNIFORM_SPARSE:
    // Needs coordinate system
    break;

    case Options::SINGLE:
    // Needs a point
    break;

    case Options::SINGLE_COPIES:
    // Needs a point
    break;

    case Options::RANDOM:
    // Needs a random number generator
    break;
  }
}

} //namespace seeding

#endif
