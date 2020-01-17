#include <vtkm/Types.h>

namespace seeding
{

namespace options = boost::program_options;
const static vtkm::FloatDefault = 100.;

class enum Seeding
{
  UNIFORM = 0,
  UNIFORM_SPARSE,
  RANDOM,
  SINGLE,
  SINGLE_COPIES
};

void ValidateSeedingOptions(options::variables_map& variables)
{

}

} //namespace seeding
