#include <iostream>
#include <string>

#include "boost/program_options.hpp"

#include <vtkm/Types.h>

int main(int argc, char** argv)
{
  namespace options = boost::program_options;
  options::options_description desc("Options");
  desc.add_options()
    ("steps", options::value<vtkm::Id>()->required(), "Number of Steps")
    ("length", options::value<vtkm::FloatDefault>()->required(), "Length of a single step")
    ("seeds", "Seeding options : UNIFORM, RANDOM");
  options::variables_map vm;
  options::store(options::parse_command_line(argc, argv, desc), vm); // can throw
  options::notify(vm);
  if(!(vm.count("steps") && vm.count("length") && vm.count("seeds")))
  {
    std::cout << "Advection Benchmark" << std::endl
              << desc << std::endl;
  }

  vtkm::Id steps = vm["steps"].as<vtkm::Id>();
  vtkm::FloatDefault length = vm["length"].as<vtkm::FloatDefault>();
  return 1;
}
