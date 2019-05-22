#ifndef SPHEREGENERATEDIR_H
#define SPHEREGENERATEDIR_H
#include "GenerateDir.h"
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/rendering/raytracing/Ray.h>
class SphereGenerateDir : public GenerateDir
{
public:
  SphereGenerateDir(vtkm::cont::ArrayHandle<vtkm::Id> _light_pointids,
                  vtkm::cont::ArrayHandle<vtkm::Float32> &radii)
    : light_pointids(_light_pointids)
  , SphereRadii(radii){}

  vtkm::cont::ArrayHandle<vtkm::Id>  light_pointids;
  vtkm::cont::ArrayHandle<vtkm::Float32> SphereRadii;

  void apply(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays);

};

#endif
