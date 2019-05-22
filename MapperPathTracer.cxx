//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2016 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
//  Copyright 2016 UT-Battelle, LLC.
//  Copyright 2016 Los Alamos National Security.
//
//  Under the terms of Contract DE-NA0003525 with NTESS,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================

#include "MapperPathTracer.h"

#include <vtkm/cont/Timer.h>
#include <vtkm/cont/TryExecute.h>

#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/Cylinderizer.h>
#include "Camera.h"
#include <vtkm/rendering/raytracing/Logger.h>
#include <vtkm/worklet/Invoker.h>
#include <vtkm/rendering/raytracing/RayOperations.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/Storage.h>
#include "PathTracer.h"
#include "PathAlgorithms.h"
#include "SurfaceWorklets.h"
#include "QuadIntersector.h"
#include "SphereIntersector.h"
#include "CosineGenerateDir.h"
#include "WhichGenerateDir.h"
#include "SphereGenerateDir.h"
#include "QuadGenerateDir.h"
#include "EmitWorklet.h"
#include "ScatterWorklet.h"
#include "QuadPdf.h"
#include "SpherePdf.h"
#include <vtkm/rendering/raytracing/SphereExtractor.h>
#include <vtkm/rendering/raytracing/QuadExtractor.h>

#include <vtkm/cont/ArrayCopy.h>

namespace vtkm
{
namespace rendering
{


struct MapperPathTracer::InternalsType
{
  vtkm::rendering::CanvasRayTracer* Canvas;
  vtkm::rendering::pathtracing::PathTracer Tracer;
  vtkm::rendering::pathtracing::Camera RayCamera;
  vtkm::rendering::raytracing::Ray<vtkm::Float32> Rays;
  bool CompositeBackground;
  VTKM_CONT
  InternalsType()
    : Canvas(nullptr)
    , CompositeBackground(true)
  {
  }
};

MapperPathTracer::MapperPathTracer(int sc, int dc,
                                   vtkm::cont::ArrayHandle<vtkm::Id> *matIdx,
                                   vtkm::cont::ArrayHandle<vtkm::Id> *texIdx,
                                   vtkm::cont::ArrayHandle<int> &matType,
                                   vtkm::cont::ArrayHandle<int> &texType)
  : Internals(new InternalsType)
  , samplecount(sc)
  , depthcount(dc)
  , MatIdx(matIdx)
  , TexIdx(texIdx)
  , MatType(matType)
  , TexType(texType)
{
  auto rays = &this->Internals->Rays;

  rays->EnableIntersectionData();
  rays->AddBuffer(1, "specular_Ox");
  rays->AddBuffer(1, "specular_Oy");
  rays->AddBuffer(1, "specular_Oz");
  rays->AddBuffer(1, "specular_Dx");
  rays->AddBuffer(1, "specular_Dy");
  rays->AddBuffer(1, "specular_Dz");
  rays->AddBuffer(1, "specular_Ax");
  rays->AddBuffer(1, "specular_Ay");
  rays->AddBuffer(1, "specular_Az");



  rays->AddBuffer(depthcount, "attenuationX");
  rays->AddBuffer(depthcount, "attenuationY");
  rays->AddBuffer(depthcount, "attenuationZ");
  rays->AddBuffer(depthcount, "emittedX");
  rays->AddBuffer(depthcount, "emittedY");
  rays->AddBuffer(depthcount, "emittedZ");
  rays->AddBuffer(1, "generated_dirX");
  rays->AddBuffer(1, "generated_dirY");
  rays->AddBuffer(1, "generated_dirZ");
  rays->AddBuffer(1, "sumtotlx");
  rays->AddBuffer(1, "sumtotly");
  rays->AddBuffer(1, "sumtotlz");
  rays->AddBuffer(1, "sum_values");

  rays->AddBuffer(depthcount, "alphaChannelAE");
  rays->AddBuffer(1, "alphaChannel");
}

MapperPathTracer::~MapperPathTracer()
{
}

void MapperPathTracer::SetCanvas(vtkm::rendering::Canvas* canvas)
{
  if (canvas != nullptr)
  {
    this->Internals->Canvas = dynamic_cast<CanvasRayTracer*>(canvas);
    if (this->Internals->Canvas == nullptr)
    {
      throw vtkm::cont::ErrorBadValue("Ray Tracer: bad canvas type. Must be CanvasRayTracer");
    }
    auto canvasSize = canvas->GetWidth() * canvas->GetHeight();
    MatIdArray.Allocate(canvasSize);
    TexIdArray.Allocate(canvasSize);
    whichPDF.Allocate(canvasSize);

  }
  else
  {
    this->Internals->Canvas = nullptr;
  }
}

vtkm::rendering::Canvas* MapperPathTracer::GetCanvas() const
{
  return this->Internals->Canvas;
}
auto MapperPathTracer::extract(const vtkm::cont::DynamicCellSet &cellset)
{

  vtkm::rendering::raytracing::SphereExtractor sphereExtractor;
  sphereExtractor.ExtractCells(cellset, 90);
  auto SphereIds = sphereExtractor.GetPointIds();
  auto SphereRadii = sphereExtractor.GetRadii();
  auto ShapeOffset = cellset.Cast<vtkm::cont::CellSetExplicit<>>().GetIndexOffsetArray(vtkm::TopologyElementTagPoint(), vtkm::TopologyElementTagCell());
  for (int i=0; i<SphereIds.GetNumberOfValues(); i++){
    std::cout << SphereIds.GetPortalConstControl().Get(i) << std::endl;
    std::cout << SphereRadii.GetPortalConstControl().Get(i) << std::endl;

  }
  vtkm::rendering::raytracing::QuadExtractor quadExtractor;
  quadExtractor.ExtractCells(cellset);
  auto QuadIds = quadExtractor.GetQuadIds();

  return std::make_tuple(SphereIds, SphereRadii, ShapeOffset, QuadIds);
}

void MapperPathTracer::RenderCells(const vtkm::cont::DynamicCellSet& cellset,
                             const vtkm::cont::CoordinateSystem& coords,
                             const vtkm::cont::Field& scalarField,
                             const vtkm::cont::ColorTable& vtkmNotUsed(colorTable),
                             const vtkm::rendering::Camera& camera,
                             const vtkm::Range& vtkmNotUsed(scalarRange))
{

  raytracing::Logger* logger = raytracing::Logger::GetInstance();
  logger->OpenLogEntry("mapper_ray_tracer");
  vtkm::cont::Timer<> tot_timer;
  vtkm::cont::Timer<> timer;

  auto canvasSize = this->Internals->Canvas->GetHeight() * this->Internals->Canvas->GetWidth();

  //
  // Add supported shapes
  //
  vtkm::Bounds shapeBounds(vtkm::Range(0,555), vtkm::Range(0,555), vtkm::Range(0,555));

  auto tup = extract(cellset);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>> QuadIds = std::get<3>(tup);
  vtkm::cont::ArrayHandle<vtkm::Id> SphereIds = std::get<0>(tup);
  vtkm::cont::ArrayHandle<vtkm::Float32> SphereRadii = std::get<1>(tup);
  auto ShapeOffset = std::get<2>(tup);

//  if (quadExtractor.GetNumberOfQuads() > 0)
//  {
//    raytracing::QuadIntersector* quadIntersector = new raytracing::QuadIntersector();
//    quadIntersector->SetData(coords, quadExtractor.GetQuadIds());
//    this->Internals->Tracer.AddShapeIntersector(quadIntersector);
//    shapeBounds.Include(quadIntersector->GetShapeBounds());
//  }

  //
  // Create rays
  //
  for (int s =0; s<samplecount; s++){
//    UVGen uvgen(nx, ny, s);
//    Invoke(uvgen, seeds, uvs);

    vtkm::rendering::pathtracing::Camera& cam = this->Internals->Tracer.GetCamera();
    cam.SetParameters(camera, *this->Internals->Canvas);
    this->Internals->RayCamera.SetParameters(camera, *this->Internals->Canvas);

    this->Internals->RayCamera.CreateRays(this->Internals->Rays, shapeBounds);
    auto rays = &this->Internals->Rays;

    //static functions in other libraries are not exported
    //so this is not callable
  //  vtkm::rendering::raytracing::RayOperations::MapCanvasToRays(
  //    this->Internals->Rays, camera, *this->Internals->Canvas);


    auto sum_values = rays->GetBuffer("sum_values").Buffer;
    auto hrecs = vtkm::rendering::pathtracing::QuadIntersector::HitRecord(rays->U, rays->V, rays->Distance, rays->NormalX, rays->NormalY, rays->NormalZ, rays->IntersectionX, rays->IntersectionY, rays->IntersectionZ);
    auto srecs = vtkm::rendering::pathtracing::QuadIntersector::ScatterRecord(
          rays->GetBuffer("specular_Ox").Buffer,
           rays->GetBuffer("specular_Oy").Buffer,
           rays->GetBuffer("specular_Oz").Buffer,
           rays->GetBuffer("specular_Dx").Buffer,
           rays->GetBuffer("specular_Dy").Buffer,
           rays->GetBuffer("specular_Dz").Buffer,
          rays->GetBuffer("specular_Ax").Buffer,
          rays->GetBuffer("specular_Ay").Buffer,
          rays->GetBuffer("specular_Az").Buffer
          );
    auto hids = vtkm::rendering::pathtracing::QuadIntersector::HitId(MatIdArray, TexIdArray);
    auto tmin = rays->MinDistance;
    using vec3CompositeType = vtkm::cont::ArrayHandleCompositeVector<
      vtkm::cont::ArrayHandle<vtkm::Float32>,
      vtkm::cont::ArrayHandle<vtkm::Float32>,
      vtkm::cont::ArrayHandle<vtkm::Float32>>;

    auto attenuation = vec3CompositeType(
        rays->GetBuffer("attenuationX").Buffer, rays->GetBuffer("attenuationY").Buffer, rays->GetBuffer("attenuationZ").Buffer);
    auto emitted = vec3CompositeType(
          rays->GetBuffer("emittedX").Buffer,rays->GetBuffer("emittedY").Buffer,rays->GetBuffer("emittedZ").Buffer);

    auto generated_dir = vec3CompositeType(
          rays->GetBuffer("generated_dirX").Buffer,rays->GetBuffer("generated_dirY").Buffer,rays->GetBuffer("generated_dirZ").Buffer);
    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,3>> tex = scalarField.GetData().Cast<vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,3>>>();
    constexpr int lightables = 2;

    for (int depth=0; depth<depthcount; depth++){
      MyAlgos::Copy<float, float, vtkm::cont::StorageTagBasic>(0, sum_values);


      intersect(coords,
                *rays,
                QuadIds,
                SphereIds,
                SphereRadii,
                MatIdArray,
                TexIdArray,
                MatIdx,
                TexIdx,
                tmin,
                emitted,
                attenuation,
                depth);

      applyMaterials(*rays,
                     hrecs,
                     hids,
                     srecs,
                     tex,
                     MatType,
                     TexType,
                     emitted,
                     cam.seeds,
                     canvasSize,
                     depth);

      vtkm::Vec<vtkm::Id, 5> tmp[] = {vtkm::Vec<vtkm::Id,5>(0,8,9,10,11)};
      auto light_box_pointids = vtkm::cont::make_ArrayHandle(tmp,1);
      vtkm::Id tmp2[] = {0};
      auto light_box_indices = vtkm::cont::make_ArrayHandle(tmp2,1);
      vtkm::Id sphere_tmp[] = {vtkm::Id(4*12)};
      auto light_sphere_pointids = vtkm::cont::make_ArrayHandle(sphere_tmp,1);
      vtkm::Id sphere_tmp2[] = {0};
      auto light_sphere_indices = vtkm::cont::make_ArrayHandle(sphere_tmp2,1);

      generateRays(coords,
                   SphereRadii,
                   whichPDF,
                   *rays,
                   this->Internals->RayCamera.seeds,
                   light_box_pointids,
                   light_box_indices,
                   light_sphere_pointids,
                   light_sphere_indices);

      applyPDFs(coords,
                QuadIds,
                SphereIds,
                SphereRadii,
                MatIdx,
                TexIdx,
                *rays,
                hrecs,
                srecs,
                sum_values,
                generated_dir,
                attenuation,
                this->Internals->RayCamera.seeds,
                light_box_pointids,
                light_box_indices,
                light_sphere_pointids,
                light_sphere_indices,
                lightables, canvasSize, depth);

      vtkm::cont::ArrayHandleCast<vtkm::Int32, vtkm::cont::ArrayHandle<vtkm::UInt8>> castedStatus(rays->Status);

    }

    using CountType = vtkm::cont::ArrayHandleCounting<vtkm::Id>;


    using vec4CompositeType = vtkm::cont::ArrayHandleCompositeVector<
      vtkm::cont::ArrayHandle<vtkm::Float32>,
      vtkm::cont::ArrayHandle<vtkm::Float32>,
      vtkm::cont::ArrayHandle<vtkm::Float32>,
      vtkm::cont::ArrayHandle<vtkm::Float32>>;

    auto attenuation4 = vec4CompositeType(
        rays->GetBuffer("attenuationX").Buffer, rays->GetBuffer("attenuationY").Buffer, rays->GetBuffer("attenuationZ").Buffer
          ,rays->GetBuffer("alphaChannelAE").Buffer);
    auto emitted4 = vec4CompositeType(
          rays->GetBuffer("emittedX").Buffer,rays->GetBuffer("emittedY").Buffer,rays->GetBuffer("emittedZ").Buffer
          ,rays->GetBuffer("alphaChannelAE").Buffer);

    auto  sumtotl = vec4CompositeType(
        rays->GetBuffer("sumtotlx").Buffer,rays->GetBuffer("sumtotly").Buffer,rays->GetBuffer("sumtotlz").Buffer
        ,rays->GetBuffer("alphaChannel").Buffer);

    vtkm::cont::ArrayHandleConstant<vtkm::Vec<vtkm::Float32,4>> zero(vtkm::Vec<vtkm::Float32,4>(0.0f), canvasSize);

    MyAlgos::SliceTransform<
        decltype(emitted4),
        decltype(zero),
        decltype(sumtotl),
        decltype(vtkm::Sum())>
        (std::make_tuple((depthcount-1)*canvasSize, (depthcount-1)*canvasSize + canvasSize), emitted4,
         std::make_tuple(0, canvasSize), zero,
         std::make_tuple(0, canvasSize), sumtotl, vtkm::Sum());

    for (int depth = depthcount-2; depth >=0; depth--){
      MyAlgos::SliceTransform<
          decltype(attenuation4),
          decltype(sumtotl),
          decltype(sumtotl),
          decltype(vtkm::Multiply())>
          (std::make_tuple(depth*canvasSize, depth*canvasSize + canvasSize), attenuation4,
           std::make_tuple(0, canvasSize), sumtotl,
           std::make_tuple(0, canvasSize), sumtotl, vtkm::Multiply());

      MyAlgos::SliceTransform<
          decltype(emitted4),
          decltype(sumtotl),
          decltype(sumtotl),
          decltype(vtkm::Sum())>
          (std::make_tuple(depth*canvasSize, depth*canvasSize + canvasSize), emitted4,
           std::make_tuple(0, canvasSize), sumtotl,
           std::make_tuple(0, canvasSize), sumtotl, vtkm::Sum());



    }
    std::cout << vtkm::cont::Algorithm::Reduce( sumtotl, vtkm::Vec<vtkm::Float32,4>(0.0)) << std::endl;

    auto cols = this->Internals->Canvas->GetColorBuffer();
    MyAlgos::Copy<vtkm::Vec<vtkm::Float32,4>, vtkm::Vec<vtkm::Float32,4>,vtkm::cont::StorageTagBasic>
        (vtkm::Vec<vtkm::Float32,4>(0.0), cols);

    vtkm::cont::Algorithm::Transform(cols, sumtotl, cols, vtkm::Sum());

    //this->Internals->Tracer.SetField(scalarField, scalarRange);

    //this->Internals->Tracer.SetColorMap(this->ColorMap);
    //this->Internals->Tracer.Render(this->Internals->Rays);

    timer.Reset();
    this->Internals->Canvas->WriteToCanvas(
      this->Internals->Rays, this->Internals->Rays.GetBuffer("default").Buffer, camera);

    if (this->Internals->CompositeBackground)
    {
      this->Internals->Canvas->BlendBackground();
    }


  }
  vtkm::Float64 time = timer.GetElapsedTime();
  logger->AddLogData("write_to_canvas", time);
  time = tot_timer.GetElapsedTime();
  logger->CloseLogEntry(time);

}

void MapperPathTracer::SetCompositeBackground(bool on)
{
  this->Internals->CompositeBackground = on;
}

void MapperPathTracer::StartScene()
{
  // Nothing needs to be done.
  auto rays = &this->Internals->Rays;
  MyAlgos::Copy<vtkm::UInt8, vtkm::UInt8, vtkm::cont::StorageTagBasic>((1UL << 3), rays->Status);

}

void MapperPathTracer::EndScene()
{
  // Nothing needs to be done.
}

vtkm::rendering::Mapper* MapperPathTracer::NewCopy() const
{
  return new vtkm::rendering::MapperPathTracer(*this);
}

template<typename emittedType,
         typename attenType>
void MapperPathTracer::intersect(const vtkm::cont::CoordinateSystem &coord,
                               vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays,
                               vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>> &QuadIds,
                               vtkm::cont::ArrayHandle<vtkm::Id> &SphereIds,
                               vtkm::cont::ArrayHandle<vtkm::Float32> &SphereRadii,
                               vtkm::cont::ArrayHandle<vtkm::Int32> &matIdArray,
                               vtkm::cont::ArrayHandle<vtkm::Int32> &texIdArray,
                               vtkm::cont::ArrayHandle<vtkm::Id> *matIdx,
                               vtkm::cont::ArrayHandle<vtkm::Id> *texIdx,
                               vtkm::cont::ArrayHandle<float> &tmin,
                               emittedType &emitted,
                               attenType &attenuation,
                               const vtkm::Id depth) const
{
  using StorageTag = vtkm::cont::StorageTagBasic;
  using MyAlgos = ::details::PathAlgorithms<vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>, VTKM_DEFAULT_DEVICE_ADAPTER_TAG>;

  MyAlgos::Copy<float, float, StorageTag>(std::numeric_limits<float>::max(), rays.Distance);
  MyAlgos::Copy<float, float, StorageTag>(0.001, tmin);
  vtkm::worklet::Invoker Invoke;

  vtkm::Id canvasSize = rays.DirX.GetNumberOfValues();



  vtkm::cont::ArrayHandle<vtkm::Int32> nodes;
  MyAlgos::Copy<vtkm::Int32, vtkm::Int32>(0, nodes);

  vtkm::rendering::pathtracing::QuadIntersector quadIntersector;

  quadIntersector.SetData(coord, QuadIds, matIdx[0], texIdx[0], matIdArray, texIdArray);
  quadIntersector.IntersectRays(rays);

  vtkm::rendering::pathtracing::SphereIntersector sphereIntersector;

  sphereIntersector.SetData(coord, SphereIds, SphereRadii, matIdx[1], texIdx[1], matIdArray, texIdArray);
  sphereIntersector.IntersectRays(rays);

  CollectIntersecttWorklet collectIntersect(canvasSize, depth);
  Invoke(collectIntersect, rays.Status, emitted, attenuation);


}

template<typename HitRecord, typename HitId, typename ScatterRecord,
         typename emittedType>
void MapperPathTracer::applyMaterials(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays,
                                    HitRecord &hrecs,
                                    HitId &hids,
                                    ScatterRecord &srecs,
                                    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,3>> tex,
                                    vtkm::cont::ArrayHandle<int> matType,
                                    vtkm::cont::ArrayHandle<int> texType,
                                    emittedType &emitted,
                                    vtkm::cont::ArrayHandle<unsigned int> &seeds,
                                    vtkm::Id canvasSize,
                                    vtkm::Id depth) const
{
  LambertianWorklet lmbWorklet( canvasSize, depth);
  DiffuseLightWorklet dlWorklet(canvasSize ,depth);
  DielectricWorklet deWorklet( canvasSize ,depth, 1.5, canvasSize);

  vtkm::worklet::Invoker Invoke;
  Invoke(lmbWorklet, rays.Origin, rays.Dir, hrecs, hids, srecs, rays.Status,
              tex, matType, texType, emitted);

  Invoke(dlWorklet, rays.Origin, rays.Dir, hrecs, hids, srecs, rays.Status,
              tex, matType, texType, emitted);

  Invoke(deWorklet, seeds, rays.Origin, rays.Dir, hrecs, hids, srecs, rays.Status,
                tex, matType, texType, emitted);

}

void MapperPathTracer::generateRays(const vtkm::cont::CoordinateSystem &coord,
                                    vtkm::cont::ArrayHandle<vtkm::Float32> &SphereRadii,
                                    vtkm::cont::ArrayHandle<int> &whichPDF,
                                    vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays,
                                    vtkm::cont::ArrayHandle<vtkm::UInt32> &seeds,
                                    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>>  light_box_pointids,
                                    vtkm::cont::ArrayHandle<vtkm::Id> light_box_indices,
                                    vtkm::cont::ArrayHandle<vtkm::Id> &light_sphere_pointids,
                                    vtkm::cont::ArrayHandle<vtkm::Id> & light_sphere_indices
                                    ) const
{
  WhichGenerateDir whichGen;
  whichGen.SetData(coord, seeds, light_sphere_indices, whichPDF);
  whichGen.apply(rays);

  CosineGenerateDir cosGen;
  cosGen.SetData(coord, seeds, light_box_indices, whichPDF);
  cosGen.apply(rays);

  QuadGenerateDir quadGen(light_box_pointids);
  quadGen.SetData(coord, seeds, light_box_indices, whichPDF);
  quadGen.apply(rays);

  SphereGenerateDir sphereGen(light_sphere_pointids, SphereRadii);
  sphereGen.SetData(coord, seeds, light_sphere_indices, whichPDF);
  sphereGen.apply(rays);
  }

template<typename HitRecord, typename ScatterRecord,
         typename attenType, typename GenDirType>
void MapperPathTracer::applyPDFs(const vtkm::cont::CoordinateSystem &coord,
                                 vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>> QuadIds,
                                 vtkm::cont::ArrayHandle<vtkm::Id> &SphereIds,
                                 vtkm::cont::ArrayHandle<vtkm::Float32> &SphereRadii,
                                 vtkm::cont::ArrayHandle<vtkm::Id> *matIdx,
                                 vtkm::cont::ArrayHandle<vtkm::Id> *texIdx,
                                 vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays,
                                 HitRecord &hrecs,
                                 ScatterRecord srecs,
                                 vtkm::cont::ArrayHandle<vtkm::Float32> &sum_values,
                                 GenDirType generated_dir,
                                 attenType &attenuation,
                                 vtkm::cont::ArrayHandle<unsigned int> &seeds,
                                 vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>>  light_box_pointids,
                                 vtkm::cont::ArrayHandle<vtkm::Id> light_box_indices,
                                 vtkm::cont::ArrayHandle<vtkm::Id> &light_sphere_pointids,
                                 vtkm::cont::ArrayHandle<vtkm::Id> & light_sphere_indices,
                                 int lightables,
                                 vtkm::Id canvasSize,
                                 vtkm::Id depth
                                 ) const
{
  vtkm::worklet::Invoker Invoke;
  QuadPdf quadPdf(QuadIds, light_box_pointids);
  quadPdf.SetData(coord, matIdx[0], texIdx[0],
      seeds, light_box_indices, lightables);
  quadPdf.apply(rays);

  SpherePdf spherePdf(SphereIds, light_sphere_pointids, SphereRadii);
  spherePdf.SetData(coord, matIdx[1], texIdx[1],
      seeds, light_sphere_indices, lightables);
  spherePdf.apply(rays);

  PDFCosineWorklet pdfWorklet(canvasSize, depth, canvasSize, lightables);
  Invoke(pdfWorklet, rays.Origin, rays.Dir, hrecs, srecs, rays.Status, sum_values, generated_dir,  rays.Origin, rays.Dir, attenuation);
}


}
} // vtkm::rendering