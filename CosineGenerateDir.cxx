#include "CosineGenerateDir.h"
#include "PdfWorklet.h"
#include <vtkm/worklet/Invoker.h>


void CosineGenerateDir::apply(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays)
{
  vtkm::worklet::Invoker Invoke;
  CosineWorketletGenerateDir cosGenDir(1);

  using vec3CompositeType = vtkm::cont::ArrayHandleCompositeVector<
    vtkm::cont::ArrayHandle<vtkm::Float32>,
    vtkm::cont::ArrayHandle<vtkm::Float32>,
    vtkm::cont::ArrayHandle<vtkm::Float32>>;
  auto generated_dir = vec3CompositeType(
        rays.GetBuffer("generated_dirX").Buffer,rays.GetBuffer("generated_dirY").Buffer,rays.GetBuffer("generated_dirZ").Buffer);

  auto hrecs = HitRecord(rays.U, rays.V, rays.Distance, rays.NormalX, rays.NormalY, rays.NormalZ, rays.IntersectionX, rays.IntersectionY, rays.IntersectionZ);

  Invoke(cosGenDir, whichPdf, hrecs, generated_dir, seeds);

}

