#ifndef validate_options_hxx
#define validate_options_hxx

#include <sstream>
#include <vector>

#include "boost/program_options.hpp"

#include <vtkm/Types.h>

#include "Config.h"

namespace validate
{

namespace options = boost::program_options;

template <typename T>
std::vector<T> Tokenize(std::string toTokenize)
{
  std::vector<T> values;
  std::string temp;
  std::stringstream tokenizer(toTokenize);
  const char delim = ':';
  while(std::getline(tokenizer, temp, delim))
  {
     std::istringstream stream(temp);
     T val;
     stream >> val;
     values.push_back(val);
  }
  return values;
}

int ValidateOptions(options::variables_map& vm,
                    config::Config& config)
{
  if(!vm.count("data"))
    return -1;
  config.SetDataSetName(vm["data"].as<std::string>());
//  if(!vm.count("field"))
//    return -1;
//  config.SetFieldName(vm["field"].as<std::string>());
  if(!vm.count("steps"))
    return -1;
  config.SetNumSteps(vm["steps"].as<vtkm::Id>());
  if(!vm.count("length"))
    return -1;
  config.SetStepLength(vm["length"].as<vtkm::FloatDefault>());
  if(!vm.count("seeds"))
    return -1;
  vtkm::Id numSeeds = vm["seeds"].as<vtkm::Id>();
  config.SetNumSeeds(numSeeds);
  if(!vm.count("seeddata"))
    return -1;
  config.SetSeedData(vm["seeddata"].as<std::string>());

/*  config::SeedingOption seeding  = static_cast<config::SeedingOption>(vm["seeding"].as<int>());
  config.SetSeeding(seeding);
  // Get seeding rake
  // If not provided / use dataset extents.
  vtkm::Bounds bounds;
  vtkm::Id3 userExtents(0,0,0);
  if(vm.count("xextent"))
  {
    userExtents[0] = 1;
    std::vector<vtkm::FloatDefault> xextent = Tokenize<vtkm::FloatDefault>(vm["xextent"].as<std::string>());
    bounds.X = vtkm::Range(xextent.at(0), xextent.at(1));
  }
  if(vm.count("yextent"))
  {
    userExtents[1] = 1;
    std::vector<vtkm::FloatDefault> yextent = Tokenize<vtkm::FloatDefault>(vm["yextent"].as<std::string>());
    bounds.Y = vtkm::Range(yextent.at(0), yextent.at(1));
  }
  if(vm.count("zextent"))
  {
    userExtents[2] = 1;
    std::vector<vtkm::FloatDefault> zextent = Tokenize<vtkm::FloatDefault>(vm["zextent"].as<std::string>());
    bounds.Z = vtkm::Range(zextent.at(0), zextent.at(1));
  }
  config.SetBounds(bounds);
  config.SetUserExtents(userExtents);
  if(seeding == config::SeedingOption::UNIFORM)
  {
    if(!vm.count("dims"))
      return -1;
    std::vector<vtkm::Id> _dims = Tokenize<vtkm::Id>(vm["dims"].as<std::string>());
    vtkm::Id3 dims(_dims.at(0), _dims.at(1), _dims.at(2));
    config.SetDimensions(dims);
  }
  else if(seeding == config::SeedingOption::RANDOM)
  {
    if(!vm.count("seeds"))
      return -1;
    vtkm::Id numSeeds = vm["seeds"].as<vtkm::Id>();
    config.SetNumSeeds(numSeeds);
  }
  else if(seeding == config::SeedingOption::SINGLE)
  {
    if(!vm.count("seeds"))
      return -1;
    if(!vm.count("point"))
      return -1;
    vtkm::Id numSeeds = vm["seeds"].as<vtkm::Id>();
    config.SetNumSeeds(numSeeds);
    std::vector<vtkm::FloatDefault> _point = Tokenize<vtkm::FloatDefault>(vm["point"].as<std::string>());
    vtkm::Vec3f point(_point.at(0), _point.at(1), _point.at(2));
    config.SetPoint(point);
  }*/
  return 0;
}

} // namespace validate

#endif

