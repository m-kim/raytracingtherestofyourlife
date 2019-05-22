#ifndef QUADPDF_H
#define QUADPDF_H
#include "Pdf.h"
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/rendering/raytracing/Ray.h>
class QuadPdf : public Pdf
{
public:
  QuadPdf(vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>>  _QuadIds,
          vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>> _light_pointids)
    :QuadIds(_QuadIds)
  , light_pointids(_light_pointids){}

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>>  QuadIds, light_pointids;

  void apply(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays);

};

#endif
