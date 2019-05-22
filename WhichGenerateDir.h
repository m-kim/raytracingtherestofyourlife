#ifndef WHICHGENERATEDIR_H
#define WHICHGENERATEDIR_H
#include "GenerateDir.h"
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/rendering/raytracing/Ray.h>
class WhichGenerateDir : public GenerateDir
{
public:
  WhichGenerateDir(){}

  void apply(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays);

};

#endif
