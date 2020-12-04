#include <vtkm/Types.h>
#include <vtkm/Math.h>

#include <vtkm/cont/Invoker.h>

#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletMapTopology.h>

namespace detail
{

class AngularEntropy : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
  using ControlSignature = void(CellSetIn, FieldInPoint, FieldOutCell, FieldOutCell);
  using ExecutionSignature = void(InputIndex, PointCount, _2, _3, _4);

  VTKM_EXEC_CONT
  AngularEntropy(vtkm::FloatDefault threshold)
  : Threshold(threshold)
  {
  }

  template<typename PointVec>
  VTKM_EXEC
  void operator()(const vtkm::IdComponent index,
                  const vtkm::IdComponent numPoints,
                  const PointVec& streamPoints,
                  vtkm::Id& pass,
                  vtkm::FloatDefault& maxCurvature) const
  {
    using Point = typename PointVec::ComponentType;
    vtkm::VecVariable<vtkm::FloatDefault, 1000> angles;
    maxCurvature = 0.;
    for(vtkm::IdComponent i = 0; i < numPoints - 2; i++)
    {
      Point v1 = streamPoints[i+1] - streamPoints[i];
      Point v2 = streamPoints[i+2] - streamPoints[i];
      auto curvature = 2*vtkm::Magnitude(vtkm::Cross(v1,v2)) /
                       (vtkm::Magnitude(streamPoints[i]   - streamPoints[i+1])*
                        vtkm::Magnitude(streamPoints[i+1] - streamPoints[i+2])*
                        vtkm::Magnitude(streamPoints[i+2] - streamPoints[i]  ));
      maxCurvature = vtkm::Max(maxCurvature, curvature);
      /*angles.Append(
           vtkm::ACos(
                 vtkm::Dot(dir1, dir2) /
                 (vtkm::Magnitude(dir1)*vtkm::Magnitude(dir2))
           )
      );*/
      //sum += angles[i];
    }
    /*vtkm::FloatDefault accumulator = 0;
    for(vtkm::IdComponent i = 0; i < numPoints - 2; i++)
    {
      auto temp = (angles[i] / sum);
      accumulator += temp*vtkm::Log2(temp);
    }
    angEntropy = -1.0*accumulator/vtkm::Log2(numPoints-1);*/
    if(maxCurvature > this->Threshold)
      pass = 1;
    else
      pass = 0;
    std::cout << index << "(" << pass << ") : " << maxCurvature << std::endl;
  }
private:
  vtkm::FloatDefault Threshold;
};

class CountAndOffset : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
  using ControlSignature = void(CellSetIn, FieldInCell, FieldInCell, FieldOutCell, FieldOutCell);
  using ExecutionSignature = void(_2, _3, _4, PointCount, _5);

  void operator()(const vtkm::Id filter,
                  const vtkm::Id inoffset,
                  vtkm::Id& outoffset,
                  const vtkm::Id incount,
                  vtkm::Id& outcount) const
  {
    if(filter == 1)
    {
      outoffset = inoffset;
      outcount  = incount;
    }
    else
    {
      outoffset = 0;
      outcount  = 0;
    }
  }
};

} //namespace detail

vtkm::cont::DataSet FilterStreamLines(const vtkm::cont::DataSet& input)
{
  vtkm::cont::Invoker invoker;

  vtkm::cont::ArrayHandle<vtkm::Id> filter;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> angularEntropy;
  vtkm::cont::DynamicCellSet cells = input.GetCellSet();
  vtkm::cont::CoordinateSystem coords = input.GetCoordinateSystem();

  std::cout << "Getting entropies" << std::endl;

  detail::AngularEntropy entropyWorklet(10000.0);
  invoker(entropyWorklet, cells, coords.GetData(), filter, angularEntropy);

  std::cout << "Finished entropies" << std::endl;

  using UnstructuredType = vtkm::cont::CellSetExplicit<>;
  UnstructuredType streams = cells.Cast<UnstructuredType>();

  vtkm::TopologyElementTagCell visitTopo{};
  vtkm::TopologyElementTagPoint incidentTopo{};

  auto inOffsets = streams.GetOffsetsArray(visitTopo, incidentTopo);
  auto inConnectivity = streams.GetConnectivityArray(visitTopo, incidentTopo);

  std::cout << "Getting Offsets and Counts" << std::endl;

  vtkm::cont::ArrayHandle<vtkm::Id> _inOffsets;
  // Truncate the last garbage value
  vtkm::cont::Algorithm::CopySubRange(inOffsets, 0, inOffsets.GetNumberOfValues() - 1, _inOffsets, 0);

  vtkm::cont::ArrayHandle<vtkm::Id> _offsets;
  vtkm::cont::ArrayHandle<vtkm::Id> _counts;
  invoker(detail::CountAndOffset(), cells, filter, _inOffsets, _offsets, _counts);
  std::cout << "Finish Offsets and Counts" << std::endl;

  vtkm::cont::ArrayHandle<vtkm::Id> connoffsets;
  vtkm::cont::ArrayHandle<vtkm::Id> counts;

  std::cout << "Getting Copyif " << std::endl;

  vtkm::cont::Algorithm::CopyIf(_offsets, filter, offsets);
  vtkm::cont::Algorithm::CopyIf(_counts, filter, counts);

  std::cout << "Finished Copyif " << std::endl;

  vtkm::Id totalStreams = vtkm::cont::Algorithm::Reduce(filter, static_cast<vtkm::Id>(0));
  vtkm::Id totalPoints  = vtkm::cont::Algorithm::Reduce(counts, static_cast<vtkm::Id>(0));

  vtkm::cont::ArrayHandle<vtkm::Id> outConnectivity;
  outConnectivity.Allocate(totalPoints);
  std::cout << "Finalizing: " << offsets.GetNumberOfValues() << ", " << counts.GetNumberOfValues() << std::endl;
  auto offsetsPortal = offsets.ReadPortal();
  auto countsPortal  = counts.ReadPortal();
  vtkm::Id runningCount = 0;
  for(vtkm::Id index = 0; index < totalStreams; index++)
  {
    auto copyOffset = offsetsPortal.Get(index);
    auto copyCount = countsPortal.Get(index);
    std::cout << "Copying" << copyOffset << " : " << copyCount << std::endl;
    vtkm::cont::Algorithm::CopySubRange(inConnectivity, copyOffset, copyCount, outConnectivity, runningCount);
    runningCount +=copyCount;
  }

  vtkm::cont::ArrayHandle<vtkm::UInt8> outCellTypes;
  auto polyLineShape =
    vtkm::cont::make_ArrayHandleConstant<vtkm::UInt8>(vtkm::CELL_SHAPE_POLY_LINE, totalStreams);
  vtkm::cont::ArrayCopy(polyLineShape, outCellTypes);

  vtkm::cont::CellSetExplicit<> outStreams;
  {
    auto numVals = offsets.GetNumberOfValues();
    vtkm::cont::ArrayHandle<vtkm::Id> _temp;
    _temp.Allocate(numVals + 1);
    vtkm::cont::Algorithm::CopySubRange(offsets, 0, numVals, _temp, 0);
    auto _tempPortal = _temp.WritePortal();
    _tempPortal.Set(numVals, runningCount);
    vtkm::cont::Algorithm::Copy(_temp, offsets);
  }
  outStreams.Fill(totalPoints, outCellTypes, outConnectivity, offsets);

  vtkm::cont::DataSet output;
  output.AddCoordinateSystem(input.GetCoordinateSystem());
  output.SetCellSet(outStreams);

  return output;
}
