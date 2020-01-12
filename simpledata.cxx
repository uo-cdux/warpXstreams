#include <vector>

#include <vtkm/Types.h>

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/DataSetFieldAdd.h>

#include <vtkm/io/writer/VTKDataSetWriter.h>

int main()
{
  vtkm::Id3 dimensions(11, 2, 2);
  vtkm::Vec3f origin(0., 0., 0.);
  vtkm::Vec3f spacing(1., 1., 1.);

  vtkm::cont::DataSetBuilderUniform datasetBuilder;
  vtkm::cont::DataSet dataset = datasetBuilder.Create(dimensions, origin, spacing);

  std::vector<vtkm::Vec3f> vectors;
  vtkm::Id numberOfPoints = dimensions[0] * dimensions[1] * dimensions[2];
  for(vtkm::Id index =0; index < numberOfPoints; index++)
  {
    vectors.push_back(vtkm::Vec3f{1., 0., 0.});
  }
  vtkm::cont::DataSetFieldAdd fieldAdder;
  fieldAdder.AddPointField(dataset, "vectors", vectors);

  vtkm::io::writer::VTKDataSetWriter writer("data/simple.vtk");
  writer.WriteDataSet(dataset);
}


