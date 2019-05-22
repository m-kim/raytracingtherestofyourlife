#include "QuadPdf.h"
#include "PdfWorklet.h"
#include <vtkm/worklet/Invoker.h>


void QuadPdf::apply(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays)
{
  vtkm::worklet::Invoker Invoke;
  QuadPDFWorklet quadPDFWorklet(this->lightables);
  QuadExecWrapper quadSurf(QuadIds, MatIdx, TexIdx);
  using vec3CompositeType = vtkm::cont::ArrayHandleCompositeVector<
    vtkm::cont::ArrayHandle<vtkm::Float32>,
    vtkm::cont::ArrayHandle<vtkm::Float32>,
    vtkm::cont::ArrayHandle<vtkm::Float32>>;
  auto generated_dir = vec3CompositeType(
        rays.GetBuffer("generated_dirX").Buffer,rays.GetBuffer("generated_dirY").Buffer,rays.GetBuffer("generated_dirZ").Buffer);

  auto sum_values = rays.GetBuffer("sum_values").Buffer;
  auto hrecs = HitRecord(rays.U, rays.V, rays.Distance, rays.NormalX, rays.NormalY, rays.NormalZ, rays.IntersectionX, rays.IntersectionY, rays.IntersectionZ);

  Invoke(quadPDFWorklet, rays.Origin, rays.Dir,hrecs,
         rays.Status, sum_values, generated_dir, seeds, quadSurf,
          light_pointids, this->light_indices, this->coordsHandle);

}


