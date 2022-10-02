#ifndef seeding_generator_hxx
#define seeding_generator_hxx

#include <random>

#include <vtkm/Particle.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/Invoker.h>
#include <vtkm/cont/RuntimeDeviceTracker.h>
#include <vtkm/worklet/WorkletMapField.h>

#include "Config.h"

namespace seeding
{

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

void MakeUniformSeeds(vtkm::Bounds bounds,
                      vtkm::Id3 dimensions,
                      vtkm::cont::ArrayHandle<vtkm::Particle>& seeds)
{
  std::cout << "Making " << dimensions << " uniform seeds" << std::endl;
  std::cout << "Bounds : " << bounds << std::endl;
  vtkm::Vec3f spacing;
  spacing[0] = bounds.X.Length() / (dimensions[0] - 1);
  spacing[1] = bounds.Y.Length() / (dimensions[1] - 1);
  spacing[2] = bounds.Z.Length() / (dimensions[2] - 1);

  std::vector<vtkm::FloatDefault> Xs;
  for(vtkm::Id i = 0; i < dimensions[0]; i++)
    Xs.push_back(bounds.X.Min + i * spacing[0]);
  std::vector<vtkm::FloatDefault> Ys;
  for(vtkm::Id i = 0; i < dimensions[1]; i++)
    Ys.push_back(bounds.Y.Min + i * spacing[1]);
  std::vector<vtkm::FloatDefault> Zs;
  for(vtkm::Id i = 0; i < dimensions[2]; i++)
    Zs.push_back(bounds.Z.Min + i * spacing[2]);

  seeds.Allocate(dimensions[0]*dimensions[1]*dimensions[2]);
  auto portal = seeds.WritePortal();
  vtkm::Id index = 0;
  for(vtkm::Id zi = 0; zi < dimensions[2]; zi++)
  {
    for(vtkm::Id yi = 0; yi < dimensions[1]; yi++)
    {
      for(vtkm::Id xi = 0; xi < dimensions[0]; xi++)
      {
         vtkm::Particle particle(vtkm::Vec3f(Xs[xi], Ys[yi], Zs[zi]), index);
         portal.Set(index, particle);
         ++index;
      }
    }
  }
}

void MakeSingleSeed(vtkm::Id seedCount,
                    vtkm::Vec3f& point,
                    vtkm::cont::ArrayHandle<vtkm::Particle>& seeds)
{
  std::cout << "Making Single Seed(s)" << std::endl;
  vtkm::cont::Invoker invoker;
  vtkm::cont::ArrayHandleIndex indices(seedCount);
  SingleSeed singleSeedWorklet(point);
  invoker(singleSeedWorklet, indices, seeds);
}

void MakeRandomSeeds(vtkm::Id seedCount,
                     vtkm::Bounds& bounds,
                     vtkm::cont::ArrayHandle<vtkm::Particle>& seeds)
{
  std::cout << "Making " << seedCount << " random seeds" << std::endl;
  std::cout << "Bounds : " << bounds << std::endl;
  std::random_device device;
  std::default_random_engine generator(static_cast<vtkm::UInt32>(255));
  vtkm::FloatDefault zero(0), one(1);
  std::uniform_real_distribution<vtkm::FloatDefault> distribution(zero, one);
  std::vector<vtkm::Particle> points;
  points.resize(0);
  for (vtkm::Id i = 0; i < seedCount; i++)
  {
    vtkm::FloatDefault rx = distribution(generator);
    vtkm::FloatDefault ry = distribution(generator);
    vtkm::FloatDefault rz = distribution(generator);
    vtkm::Vec3f p;
    p[0] = static_cast<vtkm::FloatDefault>(bounds.X.Min + rx * bounds.X.Length());
    p[1] = static_cast<vtkm::FloatDefault>(bounds.Y.Min + ry * bounds.Y.Length());
    p[2] = static_cast<vtkm::FloatDefault>(bounds.Z.Min + rz * bounds.Z.Length());
    points.push_back(vtkm::Particle(p, static_cast<vtkm::Id>(i)));
  }
  vtkm::cont::ArrayHandle<vtkm::Particle> tmp = vtkm::cont::make_ArrayHandle(points, vtkm::CopyFlag::Off);
  vtkm::cont::ArrayCopy(tmp, seeds);
}

class GetChargedParticles : public vtkm::worklet::WorkletMapField
{
public:
  GetChargedParticles(vtkm::Bounds& samplingBounds)
  : SamplingBounds(samplingBounds)
  {}

  using ControlSignature = void(FieldIn x, FieldIn y, FieldIn z,
                                FieldIn mass,
                                FieldIn charge,
                                FieldIn ux, FieldIn uy, FieldIn uz,
                                FieldIn weighting,
                                FieldOut electron,
                                FieldOut filter);

  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11);

  void operator()(const vtkm::Id index,
                  const vtkm::FloatDefault& x,
                  const vtkm::FloatDefault& y,
                  const vtkm::FloatDefault& z,
                  const vtkm::FloatDefault& mass,
                  const vtkm::FloatDefault& charge,
                  const vtkm::FloatDefault& ux,
                  const vtkm::FloatDefault& uy,
                  const vtkm::FloatDefault& uz,
                  const vtkm::FloatDefault& w,
                  vtkm::ChargedParticle& electron,
                  vtkm::Id& filter) const
  {
    constexpr static vtkm::FloatDefault SPEED_OF_LIGHT =
      static_cast<vtkm::FloatDefault>(2.99792458e8);
    auto position = vtkm::Vec3f(x, y, z);
    auto momentum = vtkm::Vec3f(ux, uy, uz);
    // Change momentum to SI units
    momentum = momentum * mass * SPEED_OF_LIGHT;
    electron = vtkm::ChargedParticle(position, index, mass, charge, w, momentum);
    if(this->SamplingBounds.Contains(position))
    {
      filter = 1;
    }
    else
    {
      filter = 0;
    }
  }

private :
  vtkm::Bounds SamplingBounds;
};

class GetChargedParticles2 : public vtkm::worklet::WorkletMapField
{
public:
  GetChargedParticles2() {}

  using ControlSignature = void(FieldIn pos,
                                FieldIn mom,
                                FieldIn mass,
                                FieldIn charge,
                                FieldIn weighting,
                                FieldOut electrons);

  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6);

  void operator()(const vtkm::Id index,
                  const vtkm::Vec3f& pos,
                  const vtkm::Vec3f& mom,
                  const vtkm::FloatDefault& mass,
                  const vtkm::FloatDefault& charge,
                  const vtkm::FloatDefault& w,
                  vtkm::ChargedParticle& electron) const
  {
    /*constexpr static vtkm::FloatDefault SPEED_OF_LIGHT =
      static_cast<vtkm::FloatDefault>(2.99792458e8);
    auto position = vtkm::Vec3f(x, y, z);
    auto momentum = vtkm::Vec3f(ux, uy, uz);
    // Change momentum to SI units
    momentum = momentum * mass * SPEED_OF_LIGHT;*/
    electron = vtkm::ChargedParticle(pos, index, mass, charge, w, mom);
  }
};


void GenerateChargedParticles(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& pos,
                              const vtkm::cont::ArrayHandle<vtkm::Vec3f>& mom,
                              const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& mass,
                              const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& charge,
                              const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& weight,
                              vtkm::cont::ArrayHandle<vtkm::ChargedParticle>& seeds)
{
  vtkm::cont::Invoker invoker;
  GetChargedParticles2 worklet;
  invoker(worklet, pos, mom, mass, charge, weight, seeds);
}

void GenerateChargedParticles(const config::Config& config,
                       const vtkm::cont::DataSet& dataset,
                       vtkm::cont::ArrayHandle<vtkm::ChargedParticle>& seeds,
                       vtkm::cont::ArrayHandle<vtkm::Id>& filter)
{
  vtkm::cont::Invoker invoker;
  vtkm::Bounds samplingBounds = config.GetBounds();
  vtkm::Id3 useSamplingBounds = config.GetUserExtents();

  vtkm::Bounds dataBounds = dataset.GetCoordinateSystem().GetBounds();
   if(useSamplingBounds[0] == 0)
     samplingBounds.X = dataBounds.X;
   if(useSamplingBounds[1] == 0)
     samplingBounds.Y = dataBounds.Y;
   if(useSamplingBounds[2] == 0)
     samplingBounds.Z = dataBounds.Z;

  GetChargedParticles worklet(samplingBounds);
  std::cout << "Sampling Bounds : " << samplingBounds << std::endl;
  //vtkm::cont::ArrayHandle<vtkm::Vec3f> positions;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> mass, charge, weighting;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> x, y, z, mom_x, mom_y, mom_z;
  dataset.GetField("x").GetData().AsArrayHandle(x);
  dataset.GetField("y").GetData().AsArrayHandle(y);
  dataset.GetField("z").GetData().AsArrayHandle(z);
  dataset.GetField("mass").GetData().AsArrayHandle(mass);
  dataset.GetField("charge").GetData().AsArrayHandle(charge);
  dataset.GetField("ux").GetData().AsArrayHandle(mom_x);
  dataset.GetField("uy").GetData().AsArrayHandle(mom_y);
  dataset.GetField("uz").GetData().AsArrayHandle(mom_z);
  dataset.GetField("w").GetData().AsArrayHandle(weighting);
  invoker(worklet, x, y, z, mass, charge, mom_x, mom_y, mom_z, weighting, seeds, filter);
}

void GenerateSeeds(const config::Config& config,
                   const vtkm::cont::DataSet& dataset,
                   vtkm::cont::ArrayHandle<vtkm::Particle>& seeds)
{
  config::SeedingOption option = config.GetSeedingOption();

  switch(option)
  {
    case config::SeedingOption::UNIFORM:
    {
      vtkm::Bounds userBounds = config.GetBounds();
      vtkm::Id3 userExtents = config.GetUserExtents();
      vtkm::Bounds dataBounds = dataset.GetCoordinateSystem().GetBounds();
      if(userExtents[0] == 0)
        userBounds.X = dataBounds.X;
      if(userExtents[1] == 0)
        userBounds.Y = dataBounds.Y;
      if(userExtents[2] == 0)
        userBounds.Z = dataBounds.Z;

      vtkm::Id3 userDimensions = config.GetDimensions();
      using CellSetType = vtkm::cont::CellSetStructured<3>;
      CellSetType cellSet;
      dataset.GetCellSet().CopyTo(cellSet);
      vtkm::Id3 dataDimensions
        = cellSet.GetSchedulingRange(vtkm::TopologyElementTagPoint());
      if(userDimensions[0] == -1)
        userDimensions[0] = dataDimensions[0];
      if(userDimensions[1] == -1)
        userDimensions[1] = dataDimensions[1];
      if(userDimensions[2] == -1)
        userDimensions[2] = dataDimensions[2];
      MakeUniformSeeds(userBounds, userDimensions, seeds);
    }
    break;

    case config::SeedingOption::RANDOM:
    // Needs a random number generator
    {
      vtkm::Id seedCount = config.GetNumSeeds();
      vtkm::Bounds bounds = dataset.GetCoordinateSystem().GetBounds();
      MakeRandomSeeds(seedCount, bounds, seeds);
    }
    break;

    case config::SeedingOption::SINGLE:
    // Needs a point
    {
      vtkm::Vec3f point = config.GetPoint();
      MakeSingleSeed(static_cast<vtkm::Id>(1), point, seeds);
    }
    break;
  }
}

}
#endif
