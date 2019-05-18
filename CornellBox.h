#ifndef CORNELLBOX_H
#define CORNELLBOX_H
#include <vtkm/cont/ArrayHandle.h>
#include "vec3.h"

class CornellBox
{
public:
  std::vector<int> cellTypeArray;
  vtkm::cont::ArrayHandle<vec3> tex;
  vtkm::cont::ArrayHandle<vtkm::Id> matIdx[5];
  vtkm::cont::ArrayHandle<vtkm::Id> texIdx[5];
  vtkm::cont::ArrayHandle<int> matType, texType;
  vtkm::cont::ArrayHandle<vec3> pts1[5];

  void invert(vtkm::Vec<vec3,4> &pts);
  void build();
};

#endif
