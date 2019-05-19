#ifndef CORNELLBOX_H
#define CORNELLBOX_H
#include <vtkm/cont/ArrayHandle.h>
#include "vec3.h"

class CornellBox
{
public:
  vtkm::cont::ArrayHandle<vec3> tex;
  vtkm::cont::ArrayHandle<vtkm::Id> matIdx[2];
  vtkm::cont::ArrayHandle<vtkm::Id> texIdx[2];
  vtkm::cont::ArrayHandle<int> matType, texType;
  vtkm::cont::ArrayHandle<vec3> pts1;

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>> QuadIds;
  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,3>> SphereIds;

  void invert(vtkm::Vec<vec3,4> &pts);
  void build();
};

#endif
