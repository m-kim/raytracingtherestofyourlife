#include "QuadGenerateDir.h"
#include "PdfWorklet.h"
#include <vtkm/worklet/Invoker.h>


void QuadGenerateDir::apply(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays)
{
  vtkm::worklet::Invoker Invoke;
  QuadWorkletGenerateDir quadGenDir(2);

  using vec3CompositeType = vtkm::cont::ArrayHandleCompositeVector<
    vtkm::cont::ArrayHandle<vtkm::Float32>,
    vtkm::cont::ArrayHandle<vtkm::Float32>,
    vtkm::cont::ArrayHandle<vtkm::Float32>>;
  auto generated_dir = vec3CompositeType(
        rays.GetBuffer("generated_dirX").Buffer,rays.GetBuffer("generated_dirY").Buffer,rays.GetBuffer("generated_dirZ").Buffer);

  auto hrecs = HitRecord(rays.U, rays.V, rays.Distance, rays.NormalX, rays.NormalY, rays.NormalZ, rays.IntersectionX, rays.IntersectionY, rays.IntersectionZ);

  Invoke(quadGenDir, this->whichPdf, hrecs, generated_dir, seeds, light_pointids, light_indices, coordsHandle);

}


