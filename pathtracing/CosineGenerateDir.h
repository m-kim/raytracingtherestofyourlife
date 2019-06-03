#ifndef COSINEGENERATEDIR_H
#define COSINEGENERATEDIR_H
#include "GenerateDir.h"
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/rendering/raytracing/Ray.h>
class CosineGenerateDir : public GenerateDir
{
public:
  CosineGenerateDir(){}

  void apply(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays);

};

#endif
