#ifndef SpherePdf_H
#define SpherePdf_H
#include "Pdf.h"
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/rendering/raytracing/Ray.h>
class SpherePdf : public Pdf
{
public:
  SpherePdf(vtkm::cont::ArrayHandle<vtkm::Id>  _SphereIds,
          vtkm::cont::ArrayHandle<vtkm::Id> _light_pointids,
            vtkm::cont::ArrayHandle<vtkm::Float32> radii)
    :SphereIds(_SphereIds)
  , light_pointids(_light_pointids)
  , SphereRadii(radii){}

  vtkm::cont::ArrayHandle<vtkm::Id>  SphereIds, light_pointids;
  vtkm::cont::ArrayHandle<vtkm::Float32>  SphereRadii;

  void apply(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays);

};

#endif
