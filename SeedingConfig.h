#ifndef seeding_config_h
#define seeding_config_h

#include <vtkm/Types.h>

namespace seeding
{
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

} //namespace seeding

#endif
