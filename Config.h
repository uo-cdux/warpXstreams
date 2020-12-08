#ifndef seeding_config_h
#define seeding_config_h

#include <vtkm/Types.h>

namespace config
{
enum class SeedingOption
{
  UNIFORM = 0,
  RANDOM  = 1,
  SINGLE  = 2,
};

class Config
{
public:
  Config()
  : Option(SeedingOption::UNIFORM)
  , Dimensions(-1, -1, -1) // Force native resolution
  {}

  void SetDataSetName(const std::string& dataSetName) {this->DataSetName = dataSetName;}
  std::string GetDataSetName() const {return this->DataSetName;}

  void SetFieldName(const std::string& fieldName) {this->FieldName = fieldName;}
  std::string GetFieldName() const {return this->FieldName;}

  void SetNumSteps(vtkm::Id numSteps) {this->NumSteps = numSteps;}
  vtkm::Id GetNumSteps() const {return this->NumSteps;}

  void SetStepLength(vtkm::FloatDefault stepLength) {this->Length = stepLength;}
  vtkm::FloatDefault GetStepLength() const {return this->Length;}

  void SetSeeding(SeedingOption option) {this->Option = option;}
  SeedingOption GetSeedingOption() const {return this->Option;}

  void SetBounds(vtkm::Bounds& bounds) {this->Bounds = bounds;}
  vtkm::Bounds GetBounds() const {return this->Bounds;}

  void SetDimensions(vtkm::Id3& dims) {this->Dimensions = dims;}
  vtkm::Id3 GetDimensions() const {return this->Dimensions;}

  void SetNumSeeds(vtkm::Id numSeeds) {this->SeedCount = numSeeds;}
  vtkm::Id GetNumSeeds() const {return this->SeedCount;}

  void SetPoint(vtkm::Vec3f point) {this->Point = point;}
  vtkm::Vec3f GetPoint() const {return this->Point;}

  void SetUserExtents(vtkm::Id3& userExtents) {this->UserExtents = userExtents;}
  vtkm::Id3 GetUserExtents() const {return this->UserExtents;}

  void SetSeedData(const std::string& seedData){this->SeedData = seedData;}
  std::string GetSeedData() const {return this->SeedData;}

  void SetThreshold(vtkm::FloatDefault threshold) {this->Threshold = threshold;}
  vtkm::FloatDefault GetThreshold() const {return this->Threshold;}
private:
  std::string DataSetName;
  std::string FieldName;
  vtkm::Id NumSteps;
  vtkm::FloatDefault Length;
  SeedingOption Option;
  vtkm::Id3 UserExtents;
  vtkm::Bounds Bounds;
  vtkm::Vec3f Point;
  vtkm::Id3 Dimensions;
  vtkm::Id SeedCount;
  std::string SeedData;
  vtkm::FloatDefault Threshold;
};

} //namespace seeding

#endif
