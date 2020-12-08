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
                  vtkm::FloatDefault& output) const
  {
    (void)index;
    using Point = typename PointVec::ComponentType;
    vtkm::VecVariable<vtkm::FloatDefault, 1000> angles;
    vtkm::FloatDefault sumCurvature = 0.;
    vtkm::FloatDefault maxCurvature = 0.;
    for(vtkm::IdComponent i = 0; i < numPoints - 2; i++)
    {
      Point v1 = streamPoints[i+1] - streamPoints[i];
      Point v2 = streamPoints[i+2] - streamPoints[i];
      auto curvature = 2*vtkm::Magnitude(vtkm::Cross(v1,v2)) /
                       (vtkm::Magnitude(streamPoints[i]   - streamPoints[i+1])*
                        vtkm::Magnitude(streamPoints[i+1] - streamPoints[i+2])*
                        vtkm::Magnitude(streamPoints[i+2] - streamPoints[i]  ));
      if(vtkm::IsNan(curvature))
        curvature = 0.0;
      sumCurvature += curvature;
      maxCurvature = vtkm::Max(static_cast<vtkm::FloatDefault>(maxCurvature),
                               static_cast<vtkm::FloatDefault>(curvature));
    }
    //vtkm::FloatDefault meanCurvature = sumCurvature / (vtkm::FloatDefault(numPoints) - 2);
    output = sumCurvature;
    if(output > this->Threshold)
      pass = 1;
    else
      pass = 0;
    //std::cout << index << "(" << pass << ") : " << maxCurvature << std::endl;
  }
private:
  vtkm::FloatDefault Threshold;
};

class CountAndOffset : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
  using ControlSignature = void(CellSetIn, FieldOutCell);
  using ExecutionSignature = void(PointCount, _2);

  void operator()(const vtkm::Id incount,
                  vtkm::Id& outcount) const
  {
      outcount  = incount;
  }
};

class Gather : public vtkm::worklet::WorkletMapField
{
public:
  Gather()
  {
  }

  using ControlSignature = void(FieldIn, WholeArrayIn, FieldOut);

  template <typename PointArrayType>
  VTKM_EXEC
  void operator()(const vtkm::Id index,
                  const PointArrayType points,
                  vtkm::Vec3f& pointAtIndex) const
  {
    auto point = points.Get(index);
    pointAtIndex = vtkm::Vec3f{point[0], point[1], point[2]};
  }
};

} //namespace detail

vtkm::cont::DataSet FilterStreamLines(const vtkm::cont::DataSet& input,
                                      const vtkm::FloatDefault& threshold)
{
  vtkm::cont::Invoker invoker;

  vtkm::cont::ArrayHandle<vtkm::Id> filter;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> maxCurvature;
  vtkm::cont::DynamicCellSet cells = input.GetCellSet();
  vtkm::cont::CoordinateSystem coords = input.GetCoordinateSystem();

  detail::AngularEntropy entropyWorklet(threshold);
  invoker(entropyWorklet, cells, coords.GetData(), filter, maxCurvature);

  {
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> _maxCurvature;
    vtkm::cont::Algorithm::Copy(maxCurvature, _maxCurvature);
    vtkm::cont::Algorithm::Sort(_maxCurvature);
    vtkm::Id values = _maxCurvature.GetNumberOfValues();
    auto portal = _maxCurvature.ReadPortal();
    std::cout << "Curvature (Min/Max) : " << portal.Get(0) << "/" << portal.Get(values-1) << std::endl;
    std::cout << "Curvature 10% : " << portal.Get(90*(vtkm::FloatDefault(values)/100.)) << std::endl;
    std::cout << "Curvature 20% : " << portal.Get(80*(vtkm::FloatDefault(values)/100.)) << std::endl;
    std::cout << "Curvature 50% : " << portal.Get(50*(vtkm::FloatDefault(values)/100.)) << std::endl;
  }

  using UnstructuredType = vtkm::cont::CellSetExplicit<>;
  UnstructuredType streams = cells.Cast<UnstructuredType>();

  vtkm::TopologyElementTagCell visitTopo{};
  vtkm::TopologyElementTagPoint incidentTopo{};

  auto inOffsets = streams.GetOffsetsArray(visitTopo, incidentTopo);
  auto inConnectivity = streams.GetConnectivityArray(visitTopo, incidentTopo);

  vtkm::cont::ArrayHandle<vtkm::Id> inCounts;
  invoker(detail::CountAndOffset(), cells, inCounts);

  vtkm::cont::ArrayHandle<vtkm::Id> offsets;
  vtkm::cont::ArrayHandle<vtkm::Id> counts;

  vtkm::cont::Algorithm::CopyIf(inOffsets, filter, offsets);
  vtkm::cont::Algorithm::CopyIf(inCounts, filter, counts);

  vtkm::Id totalStreams = vtkm::cont::Algorithm::Reduce(filter, static_cast<vtkm::Id>(0));
  vtkm::Id totalPoints  = vtkm::cont::Algorithm::Reduce(counts, static_cast<vtkm::Id>(0));

  vtkm::cont::ArrayHandle<vtkm::Vec3f> outCoords;
  vtkm::cont::ArrayHandle<vtkm::Id> outConnectivity;
  outConnectivity.Allocate(totalPoints);
  auto offsetsPortal = offsets.ReadPortal();
  auto countsPortal  = counts.ReadPortal();
  vtkm::Id runningCount = 0;
  for(vtkm::Id index = 0; index < totalStreams; index++)
  {
    auto copyOffset = offsetsPortal.Get(index);
    auto copyCount = countsPortal.Get(index);
    vtkm::cont::Algorithm::CopySubRange(inConnectivity, copyOffset, copyCount, outConnectivity, runningCount);
    runningCount +=copyCount;
  }

  vtkm::cont::ArrayHandle<vtkm::UInt8> outCellTypes;
  auto polyLineShape =
    vtkm::cont::make_ArrayHandleConstant<vtkm::UInt8>(vtkm::CELL_SHAPE_POLY_LINE, totalStreams);
  vtkm::cont::ArrayCopy(polyLineShape, outCellTypes);

  // Compress Coordinate System
  {
    vtkm::cont::ArrayHandle<vtkm::Id> _outConnectivity;
    vtkm::cont::Algorithm::Copy(outConnectivity, _outConnectivity);
    vtkm::cont::Algorithm::Unique(_outConnectivity);

    invoker(detail::Gather{}, _outConnectivity, input.GetCoordinateSystem().GetData(), outCoords);

    vtkm::cont::Algorithm::SortByKey(_outConnectivity, outCoords);

    vtkm::cont::ArrayHandle<vtkm::Id> newConnectivity;
    vtkm::cont::Algorithm::LowerBounds(_outConnectivity, outConnectivity, newConnectivity);
    vtkm::cont::Algorithm::Copy(newConnectivity, outConnectivity);
  }

  vtkm::cont::CellSetExplicit<> outStreams;
  {
    vtkm::cont::Algorithm::ScanExclusive(counts, offsets);
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
  output.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coords", outCoords));
  output.SetCellSet(outStreams);

  return output;
}
