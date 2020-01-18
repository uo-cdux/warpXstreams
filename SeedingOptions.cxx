#include <vtkm/Types.h>

namespace seeding
{

namespace options = boost::program_options;
const static vtkm::FloatDefault = 100.;

class enum Seeding
{
  UNIFORM = 0,
  UNIFORM_SPARSE,
  SINGLE,
  SINGLE_COPIES
  RANDOM,
};

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


bool ValidateSeedingOptions(options::variables_map& variables)
{
  if(!vm.count("seeding"))
    return false;
  Seeding option = = static_cast<>(i);
}

void GenerateSeeds(int option,
                   vtkm::cont::DataSet& dataset,
                   vtkm::cont::ArrayHandle<vtkm::Vec3f>& seeds)
{
  switch(option)
  {
    case Seeding::UNIFORM:
    // Needs coordinate system
    break;

    case Seeding::UNIFORM_SPARSE:
    // Needs coordinate system
    break;

    case Seeding::SINGLE:
    // Needs a point
    break;

    case Seeding::SINGLE_COPIES:
    // Needs a point
    break;

    case Seeding::RANDOM:
    // Needs a random number generator
    break;
  }
}


} //namespace seeding
