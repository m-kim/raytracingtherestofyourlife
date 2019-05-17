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
  vtkm::cont::ArrayHandle<int> matType, texType, flipped[5];
  vtkm::cont::ArrayHandle<vec3> pts1[5], pts2[5], translateOffset[5];
  vtkm::cont::ArrayHandle<float> angleArray[5];
  vtkm::cont::ArrayHandle<vec3> norms;
  void build();
};

#endif
