#ifndef seeding_options_hxx
#define seeding_options_hxx

#include <vector>

#include "boost/program_options.hpp"

#include <vtkm/Types.h>

#include "SeedingConfig.h"

namespace validate
{

namespace options = boost::program_options;

bool ValidateSeedingOptions(options::variables_map& vm,
                            seeding::SeedingConfig& config)
{
  if(!vm.count("seeding"))
    return false;
  int param = vm["seeding"].as<int>();
  seeding::Options option = static_cast<seeding::Options>(param);
//  std::cout << "Option : " << param << std::endl;
  config.SetOption(option);
  switch(option)
  {
    case seeding::Options::UNIFORM:
    {
//      std::cout << "uniform seeding" << std::endl;
      return true;
    }
    break;

    case seeding::Options::UNIFORM_SPARSE:
    {
//      std::cout << "uniform sparse seeding" << std::endl;
      std::vector<vtkm::Id> vec = vm["density"].as<std::vector<vtkm::Id>>();
      vtkm::Id3 density{vec[0], vec[1], vec[2]};
      config.SetDensity(density);
      return true;
    }
    break;

    case seeding::Options::SINGLE:
    // Needs a point
    {
 //     std::cout << "single seeding" << std::endl;
      std::vector<vtkm::FloatDefault> vec = vm["point"].as<std::vector<vtkm::FloatDefault>>();
      vtkm::Vec3f point{vec[0], vec[1], vec[2]};
      config.SetPoint(point);
      return true;
    }
    break;

    case seeding::Options::SINGLE_COPIES:
    // Needs a point
    {
//      std::cout << "single seeding" << std::endl;
      std::vector<vtkm::FloatDefault> vec = vm["point"].as<std::vector<vtkm::FloatDefault>>();
      vtkm::Vec3f point{vec[0], vec[1], vec[2]};
      config.SetPoint(point);
      vtkm::Id seedCount = vm["seeds"].as<vtkm::Id>();
      config.SetSeedCount(seedCount);
      return true;
    }
    break;

    case seeding::Options::RANDOM:
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

} // namespace validate



#endif
