#ifndef QUADGENERATEDIR_H
#define QUADGENERATEDIR_H
#include "GenerateDir.h"
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/rendering/raytracing/Ray.h>
class QuadGenerateDir : public GenerateDir
{
public:
  QuadGenerateDir(vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>> _light_pointids)
    : light_pointids(_light_pointids){}

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>>  light_pointids;

  void apply(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays);

};

#endif
