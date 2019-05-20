#ifndef CORNELLBOX_H
#define CORNELLBOX_H
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSet.h>
#include "vec3.h"
#include <vtkm/cont/CoordinateSystem.h>

class CornellBox
{
public:
  vtkm::cont::ArrayHandle<vec3> tex;
  vtkm::cont::ArrayHandle<vtkm::Id> matIdx[2];
  vtkm::cont::ArrayHandle<vtkm::Id> texIdx[2];
  vtkm::cont::ArrayHandle<int> matType, texType;
  vtkm::cont::CoordinateSystem coord;

  vtkm::cont::ArrayHandle<vtkm::UInt8> shapes[2];
  vtkm::cont::ArrayHandle<vtkm::Id> ptsIdx[2];
  vtkm::cont::ArrayHandle<vtkm::IdComponent> numindices[2];

  vtkm::cont::ArrayHandle<vtkm::Id> conn[2];

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>> QuadIds;
  vtkm::cont::ArrayHandle<vtkm::Id> SphereIds;

  vtkm::cont::ArrayHandle<vtkm::Float32> SphereRadii;

  void invert(vtkm::Vec<vec3,4> &pts);
  void build();

  vtkm::cont::DataSet buildDataSet();
};

#endif
