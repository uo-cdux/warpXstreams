#ifndef seeding_generator_hxx
#define seeding_generator_hxx

#include <vtkm/Particle.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/Invoker.h>
#include <vtkm/cont/RuntimeDeviceTracker.h>
#include <vtkm/worklet/WorkletMapField.h>

#include <vtkm/Particle.h>

#include "SeedingConfig.h"

namespace seeding
{

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

class SingleSeed : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  SingleSeed(vtkm::Vec3f point)
  : Point(point)
  {}

  using ControlSignature = void(FieldIn, FieldOut);

  VTKM_EXEC
  void operator()(const vtkm::Id index,
                  vtkm::Particle& particle) const
  {
    particle.ID = index;
    particle.Pos = this->Point;
  }

private:
  vtkm::Vec3f Point;
};

void MakeUniformSeeds(vtkm::cont::CoordinateSystem& coords,
                      vtkm::cont::ArrayHandle<vtkm::Particle>& seeds)
{
  vtkm::cont::Invoker invoker;
  invoker(SeedsFromCoordinates{}, coords.GetData(), seeds);
}

// TODO : do better(lol)
/*void MakeSparseUniformSeeds(vtkm::Bounds& bounds,
                            vtkm::Id3 sparsity,
                            vtkm::Id3 )*/

void MakeSingleSeed(vtkm::Id seedCount,
                    vtkm::Vec3f& point,
                    vtkm::cont::ArrayHandle<vtkm::Particle>& seeds)
{
  vtkm::cont::Invoker invoker;
  vtkm::cont::ArrayHandleIndex indices(seedCount);
  SingleSeed singleSeedWorklet(point);
  invoker(singleSeedWorklet, indices, seeds);
}

void MakeRandomSeeds(vtkm::Id seedCount,
                     vtkm::cont::ArrayHandle<vtkm::Particle>& seeds)
{
  (void)seedCount;
  (void)seeds;
  // Not implemented yet.
}

class GetElectrons : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn position,
                                FieldIn mass,
                                FieldIn charge,
                                FieldIn ux,
                                FieldIn uy,
                                FieldIn uz,
                                FieldIn weighting,
                                FieldOut);

  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8);

  void operator()(const vtkm::Id index,
                  const vtkm::Vec3f& position,
                  const vtkm::FloatDefault& mass,
                  const vtkm::FloatDefault& charge,
                  const vtkm::FloatDefault& ux,
                  const vtkm::FloatDefault& uy,
                  const vtkm::FloatDefault& uz,
                  const vtkm::FloatDefault& w,
                  vtkm::Electron& electron) const
  {
    electron = vtkm::Electron(position, index, mass, charge, w, vtkm::Vec3f(ux,uy,uz));
  }
};

void GenerateElectrons(vtkm::cont::DataSet& dataset,
                       vtkm::cont::ArrayHandle<vtkm::Electron>& seeds)
{
  vtkm::cont::Invoker invoker;
  //vtkm::cont::ArrayHandle<vtkm::Vec3f> positions;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> mass, charge, mom_x, mom_y, mom_z, weighting;
  auto positions = dataset.GetCoordinateSystem().GetData();
  dataset.GetField("mass").GetData().CopyTo(mass);
  dataset.GetField("charge").GetData().CopyTo(charge);
  dataset.GetField("ux").GetData().CopyTo(mom_x);
  dataset.GetField("uy").GetData().CopyTo(mom_y);
  dataset.GetField("uz").GetData().CopyTo(mom_z);
  dataset.GetField("w").GetData().CopyTo(weighting);
  invoker(GetElectrons{}, positions, mass, charge, mom_x, mom_y, mom_z, weighting, seeds);
}


void GenerateSeeds(SeedingConfig& config,
                   vtkm::cont::DataSet& dataset,
                   vtkm::cont::ArrayHandle<vtkm::Particle>& seeds)
{
  Options option = config.GetOption();
  switch(option)
  {
    case Options::UNIFORM:
    {
      vtkm::cont::CoordinateSystem coords = dataset.GetCoordinateSystem();
      MakeUniformSeeds(coords, seeds);
    }
    break;

    case Options::UNIFORM_SPARSE:
    // Needs coordinate system
    {
      //TODO: Take care later
    }
    break;

    case Options::SINGLE:
    // Needs a point
    {
      vtkm::Vec3f point = config.GetPoint();
      MakeSingleSeed(static_cast<vtkm::Id>(1), point, seeds);
    }
    break;

    case Options::SINGLE_COPIES:
    // Needs a point
    {
      vtkm::Id seedCount = config.GetSeedCount();
      vtkm::Vec3f point = config.GetPoint();
      MakeSingleSeed(seedCount, point, seeds);
    }
    break;

    case Options::RANDOM:
    // Needs a random number generator
    break;
  }
}

}

#endif
