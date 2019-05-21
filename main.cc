//==================================================================================================
// Written in 2019 by Mark Kim
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is distributed
// without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication along
// with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==================================================================================================

#include <iostream>
#include <limits>
#include <vector>
#include <tuple>
#include <vtkm/cont/ArrayHandleConstant.h>
#include "Camera.h"
#include <vtkm/cont/internal/DeviceAdapterAlgorithmGeneral.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleCompositeVector.h>
#include <vtkm/rendering/raytracing/Ray.h>
#include <vtkm/rendering/raytracing/RayOperations.h>
#include <omp.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/worklet/Invoker.h>
#include <vtkm/rendering/raytracing/BoundingVolumeHierarchy.h>
#include "BVHTraverser.h"

#include <fstream>
#include "Worklets.h"
#include "SurfaceWorklets.h"
#include "EmitWorklet.h"
#include "ScatterWorklet.h"
#include "PdfWorklet.h"

#include "CornellBox.h"
#include "PathAlgorithms.h"
#include "AABBSurface.h"


using ArrayType = vtkm::cont::ArrayHandle<vec3>;

inline vec3 de_nan(const vec3& c) {
    vec3 temp = c;
    if (!(temp[0] == temp[0])) temp[0] = 0;
    if (!(temp[1] == temp[1])) temp[1] = 0;
    if (!(temp[2] == temp[2])) temp[2] = 0;
    return temp;
}


struct Append
{
  template <typename T>
  VTKM_EXEC_CONT T operator()(const T& x, const T& y) const
  {
    return x + y;
  }
};


vtkm::rendering::raytracing::LinearBVH bvhSphere, bvhQuad;


using RayType = vtkm::rendering::raytracing::Ray<float>;

void buildBVH(
              CornellBox &cb)
{

  vtkm::cont::CoordinateSystem coords = cb.coord;
  vtkm::rendering::raytracing::AABBs aabbQuads;
  vtkm::worklet::DispatcherMapField<::detail::FindQuadAABBs>(::detail::FindQuadAABBs())
    .Invoke(cb.QuadIds,
            aabbQuads.xmins,
            aabbQuads.ymins,
            aabbQuads.zmins,
            aabbQuads.xmaxs,
            aabbQuads.ymaxs,
            aabbQuads.zmaxs,
            cb.coord);

  bvhQuad.SetData(aabbQuads);
  bvhQuad.Construct();


  vtkm::rendering::raytracing::AABBs aabbSpheres;
  vtkm::worklet::DispatcherMapField<::detail::FindSphereAABBs>(::detail::FindSphereAABBs())
    .Invoke(cb.SphereIds,
            cb.SphereRadii,
            aabbSpheres.xmins,
            aabbSpheres.ymins,
            aabbSpheres.zmins,
            aabbSpheres.xmaxs,
            aabbSpheres.ymaxs,
            aabbSpheres.zmaxs,
            cb.coord);

  bvhSphere.SetData(aabbSpheres);
  bvhSphere.Construct();


  //ShapeBounds = bvh.TotalBounds;


}
template<typename HitRecord,
         typename HitId,
         typename emittedType,
         typename attenType>
void intersect(CornellBox &cb,
               RayType &rays,
               HitRecord &hrecs,
               HitId &hids,
               vtkm::cont::ArrayHandle<float> &tmin,
               emittedType &emitted,
               attenType &attenuation,
               const vtkm::Id depth)
{
  using MyAlgos = details::PathAlgorithms<vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>, VTKM_DEFAULT_DEVICE_ADAPTER_TAG>;
  using StorageTag = vtkm::cont::StorageTagBasic;

  MyAlgos::Copy<float, float, StorageTag>(std::numeric_limits<float>::max(), rays.Distance);
  MyAlgos::Copy<float, float, StorageTag>(0.001, tmin);
  vtkm::worklet::Invoker Invoke;

  vtkm::Id canvasSize = rays.DirX.GetNumberOfValues();



  vtkm::cont::ArrayHandle<vtkm::Int32> nodes;
  nodes.Allocate(rays.Dir.GetNumberOfValues());
  for (int j=0; j<nodes.GetNumberOfValues(); j++){
    nodes.GetPortalControl().Set(j, 0);
  }
  QuadIntersect quad;
  SphereIntersecttWorklet sphereWorklet(canvasSize, depth);


  QuadExecWrapper quadIntersect(cb.QuadIds, cb.matIdx[0], cb.texIdx[0]);
  SphereExecWrapper sphereIntersect(cb.SphereIds, cb.SphereRadii, cb.matIdx[1], cb.texIdx[1]);

#if 0
   vtkm::cont::ArrayHandle<vtkm::Id> leafs;
  leafs.Allocate(13);
  leafs.GetPortalControl().Set(0, 12);
  for (int i=0; i<12; i++){

    leafs.GetPortalControl().Set(i+1, i);
  }

  Invoke(quad,
         nodes,
         rays.Origin,
         rays.Dir,
         hrecs,
         hids,
         tmin,
         rays.Distance,
         rays.Status,
         quadIntersect,
         cb.coord,
         leafs);
  vtkm::cont::ArrayHandle<vtkm::Id> sphereleafs;
  sphereleafs.Allocate(2);
  sphereleafs.GetPortalControl().Set(0, 1);
  sphereleafs.GetPortalControl().Set(1, 0);

  Invoke(sphereWorklet,
         nodes,
         rays.Origin,
         rays.Dir,
         hrecs,
         hids,
         tmin,
         rays.Distance,
         rays.Status,
         surf,
         cb.coord,
         sphereleafs);
#else
  vtkm::rendering::pathtracing::BVHTraverser traverser;
  traverser.IntersectRays(rays, bvhQuad, hrecs, hids, tmin, quadIntersect, cb.coord);

  vtkm::rendering::pathtracing::BVHTraverser traverser2;
  traverser2.IntersectRays(rays, bvhSphere, hrecs, hids, tmin, sphereIntersect, cb.coord);

#endif



  CollectIntersecttWorklet collectIntersect(canvasSize, depth);
  Invoke(collectIntersect, rays.Status, emitted, attenuation);

}

template<typename HitRecord, typename HitId, typename ScatterRecord,
         typename emittedType>
void applyMaterials(RayType &rays,
                    HitRecord &hrecs,
                    HitId &hids,
                    ScatterRecord &srecs,
                    vtkm::cont::ArrayHandle<vec3> tex,
                    vtkm::cont::ArrayHandle<int> matType,
                    vtkm::cont::ArrayHandle<int> texType,
                    emittedType &emitted,
                    vtkm::cont::ArrayHandle<unsigned int> &seeds,
                    vtkm::Id canvasSize,
                    vtkm::Id depth)
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

template<typename HitRecord,
         typename GenDirType>
void generateRays(
    vtkm::cont::ArrayHandle<int> &whichPDF,
    HitRecord &hrecs,
    GenDirType &generated_dir,
    vtkm::cont::ArrayHandle<unsigned int> &seeds,
    std::vector<ArrayType> &light_box_pts,
    std::vector<ArrayType> &light_sphere_pts,
    std::vector<vtkm::cont::ArrayHandle<vtkm::Float32>> &light_sphere_radii
    )
{
  vtkm::worklet::Invoker Invoke;
  GenerateDir genDir(3);
  CosineGenerateDir cosGenDir(1);
  XZRectGenerateDir xzRectGenDir(2);
  SphereGenerateDir sphereGenDir(3);

  Invoke(genDir, seeds, whichPDF);
  Invoke(cosGenDir, whichPDF, hrecs, generated_dir, seeds);
  Invoke(xzRectGenDir, whichPDF, hrecs, generated_dir, seeds, light_box_pts[0], light_box_pts[1]);
  Invoke(sphereGenDir, whichPDF, hrecs, generated_dir, seeds, light_sphere_pts[0], light_sphere_radii[0]);
  }

template<typename HitRecord, typename ScatterRecord,
         typename attenType, typename GenDirType>
void applyPDFs(CornellBox &cb,
    RayType &rays,
    HitRecord &hrecs,
    ScatterRecord srecs,
    vtkm::cont::ArrayHandle<vtkm::Float32> &sum_values,
    GenDirType generated_dir,
    attenType &attenuation,
    vtkm::cont::ArrayHandle<unsigned int> &seeds,
    std::vector<ArrayType> light_box_pts,
    std::vector<ArrayType> light_sphere_pts,
    std::vector<vtkm::cont::ArrayHandle<vtkm::Float32>> light_sphere_radii,
    int lightables,
    vtkm::Id canvasSize,
    vtkm::Id depth
    )
{
  vtkm::worklet::Invoker Invoke;
  XZRectPDFWorklet xzPDFWorklet(lightables);
  SpherePDFWorklet spherePDFWorklet(lightables);
  PDFCosineWorklet pdfWorklet(canvasSize, depth, canvasSize, lightables);
  XZRectExecWrapper xzsurf;
  Invoke(xzPDFWorklet, rays.Origin, rays.Dir,hrecs, rays.Status, sum_values, generated_dir, seeds,xzsurf, light_box_pts[0], light_box_pts[1]);
  SphereExecWrapper surf(cb.SphereIds, cb.SphereRadii, cb.matIdx[1], cb.texIdx[1]);

  Invoke(spherePDFWorklet, rays.Origin, rays.Dir,hrecs, rays.Status, sum_values, generated_dir, seeds, surf, light_sphere_pts[0], light_sphere_radii[0]);

  Invoke(pdfWorklet, rays.Origin, rays.Dir, hrecs, srecs, rays.Status, sum_values, generated_dir,  rays.Origin, rays.Dir, attenuation);

}
template <typename T, typename U, class CIn, class COut, class BinaryFunctor>
VTKM_CONT static T MyScanInclusive(const vtkm::cont::ArrayHandle<T, CIn>& input,
                                 vtkm::cont::ArrayHandle<U, COut>& output)
{
  vtkm::cont::detail::ScanInclusiveResultFunctor<U> functor;
  vtkm::cont::TryExecute(functor, input, output);
  return functor.result;
}

template<typename DeviceAdapterTag>
void resizeRays(vtkm::rendering::raytracing::Ray<float> &rays, int canvasSize)
{
  rays.Resize(canvasSize);
}

const auto
parse(int argc, char **argv){
  int x = 128;
  int y = 128;
  int s = 10;
  int depth = 5;

  for (int i=1; i<argc; i++){
    if (!strcmp(argv[i], "-x")){
      if (i+1 < argc){
        x = atoi(argv[i+1]);
        i += 1;
      }

    }
    else if (!strcmp(argv[i], "-y")){
      if (i+1 < argc){
        y = atoi(argv[i+1]);
        i += 1;
      }
    }
    else if (!strcmp(argv[i], "-samplecount")){
      if (i+1 < argc){
        s = atoi(argv[i+1]);
        i += 1;
      }
    }
    else if (!strcmp(argv[i], "-raydepth")){
      if (i+1 < argc){
        depth = atoi(argv[i+1]);
        i += 1;
      }
    }
  }

  return std::make_tuple(x,y, s, depth);
}

int main(int argc, char *argv[]) {
  using MyAlgos = details::PathAlgorithms<vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>, VTKM_DEFAULT_DEVICE_ADAPTER_TAG>;
  using StorageTag = vtkm::cont::StorageTagBasic;
  using Device = VTKM_DEFAULT_DEVICE_ADAPTER_TAG;

  const auto tup = parse(argc, argv);
  const int nx = std::get<0>(tup);
  const int ny = std::get<1>(tup);
  const int ns = std::get<2>(tup);
  const int depthcount = std::get<3>(tup);

  CornellBox cb;

  auto canvasSize = nx*ny;

  constexpr int lightables = 2;
  std::vector<ArrayType> light_box_pts(2), light_sphere_pts(1);
  std::vector<vtkm::cont::ArrayHandle<vtkm::Float32>> light_sphere_radii(1);
  light_box_pts[0].Allocate(1);
  light_box_pts[0].GetPortalControl().Set(0, vec3(213, 554, 227));
  light_box_pts[1].Allocate(1);
  light_box_pts[1].GetPortalControl().Set(0, vec3(343, 554, 332));

  light_sphere_pts[0].Allocate(1);
  light_sphere_pts[0].GetPortalControl().Set(0, vec3(190, 90, 190));
  light_sphere_radii[0].Allocate(1);
  light_sphere_radii[0].GetPortalControl().Set(0, 90);

  auto ds = cb.buildDataSet();

  cb.extract();
  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 2>> uvs;
  uvs.Allocate(nx*ny);

  vtkm::cont::ArrayHandle<vtkm::Id> rays_index;

  vtkm::rendering::raytracing::Ray<float> rays;
  rays.EnableIntersectionData();
  vtkm::rendering::raytracing::RayOperations::Resize(rays, canvasSize, Device());


  MyAlgos::Copy<vtkm::UInt8, vtkm::UInt8, StorageTag>((1UL << 3), rays.Status);

  vtkm::cont::ArrayHandleConstant<vec3> zero(vec3(0.0f), nx*ny);
  vtkm::cont::ArrayHandle<vec3> cols;
  cols.Allocate(nx*ny);
  MyAlgos::Copy<vec3, vec3, StorageTag>(vec3(0.0f), cols);

  vtkm::cont::ArrayHandle<vtkm::Float32> sum_values;

  vtkm::cont::ArrayHandle<vtkm::Int32> matIdArray, texIdArray;
  matIdArray.Allocate(canvasSize);
  texIdArray.Allocate(canvasSize);

  rays.AddBuffer(1, "specular_Ox");
  rays.AddBuffer(1, "specular_Oy");
  rays.AddBuffer(1, "specular_Oz");
  rays.AddBuffer(1, "specular_Dx");
  rays.AddBuffer(1, "specular_Dy");
  rays.AddBuffer(1, "specular_Dz");
  rays.AddBuffer(1, "specular_Ax");
  rays.AddBuffer(1, "specular_Ay");
  rays.AddBuffer(1, "specular_Az");
  using ScatterRecord = vtkm::cont::ArrayHandleCompositeVector<vtkm::cont::ArrayHandle<vtkm::Float32>,
  vtkm::cont::ArrayHandle<vtkm::Float32>,
  vtkm::cont::ArrayHandle<vtkm::Float32>,
  vtkm::cont::ArrayHandle<vtkm::Float32>,
  vtkm::cont::ArrayHandle<vtkm::Float32>,
  vtkm::cont::ArrayHandle<vtkm::Float32>,
  vtkm::cont::ArrayHandle<vtkm::Float32>,
  vtkm::cont::ArrayHandle<vtkm::Float32>,
  vtkm::cont::ArrayHandle<vtkm::Float32>>;
  auto srecs = ScatterRecord(
        rays.GetBuffer("specular_Ox").Buffer,
         rays.GetBuffer("specular_Oy").Buffer,
         rays.GetBuffer("specular_Oz").Buffer,
         rays.GetBuffer("specular_Dx").Buffer,
         rays.GetBuffer("specular_Dy").Buffer,
         rays.GetBuffer("specular_Dz").Buffer,
        rays.GetBuffer("specular_Ax").Buffer,
        rays.GetBuffer("specular_Ay").Buffer,
        rays.GetBuffer("specular_Az").Buffer
        );
  using HitRecord = vtkm::cont::ArrayHandleCompositeVector<decltype(rays.U),
  decltype(rays.V),
  decltype(rays.Distance),
  decltype(rays.NormalX),
  decltype(rays.NormalY),
  decltype(rays.NormalZ),
  decltype(rays.IntersectionX),
  decltype(rays.IntersectionY),
  decltype(rays.IntersectionZ)>;

  auto hrecs = HitRecord(rays.U, rays.V, rays.Distance, rays.NormalX, rays.NormalY, rays.NormalZ, rays.IntersectionX, rays.IntersectionY, rays.IntersectionZ);

  using HitId = vtkm::cont::ArrayHandleCompositeVector<decltype(matIdArray), decltype(texIdArray)>;
  auto hids = HitId(matIdArray, texIdArray);

  rays.AddBuffer(depthcount, "attenuationX");
  rays.AddBuffer(depthcount, "attenuationY");
  rays.AddBuffer(depthcount, "attenuationZ");
  rays.AddBuffer(depthcount, "emittedX");
  rays.AddBuffer(depthcount, "emittedY");
  rays.AddBuffer(depthcount, "emittedZ");
  rays.AddBuffer(1, "generated_dirX");
  rays.AddBuffer(1, "generated_dirY");
  rays.AddBuffer(1, "generated_dirZ");
  rays.AddBuffer(1, "sumtotlx");
  rays.AddBuffer(1, "sumtotly");
  rays.AddBuffer(1, "sumtotlz");
  rays.AddBuffer(1, "sum_values");
  using vec3CompositeType = vtkm::cont::ArrayHandleCompositeVector<
    vtkm::cont::ArrayHandle<vtkm::Float32>,
    vtkm::cont::ArrayHandle<vtkm::Float32>,
    vtkm::cont::ArrayHandle<vtkm::Float32>>;
  auto attenuation = vec3CompositeType(
      rays.GetBuffer("attenuationX").Buffer, rays.GetBuffer("attenuationY").Buffer, rays.GetBuffer("attenuationZ").Buffer);
  auto emitted = vec3CompositeType(
        rays.GetBuffer("emittedX").Buffer,rays.GetBuffer("emittedY").Buffer,rays.GetBuffer("emittedZ").Buffer);

  auto generated_dir = vec3CompositeType(
        rays.GetBuffer("generated_dirX").Buffer,rays.GetBuffer("generated_dirY").Buffer,rays.GetBuffer("generated_dirZ").Buffer);
  vtkm::cont::ArrayHandle<int> whichPDF;

  auto  sumtotl = vec3CompositeType(
      rays.GetBuffer("sumtotlx").Buffer,rays.GetBuffer("sumtotly").Buffer,rays.GetBuffer("sumtotlz").Buffer);

  auto tmin = rays.MinDistance;

  sum_values = rays.GetBuffer("sum_values").Buffer;
  whichPDF.Allocate(nx*ny);
  vtkm::rendering::pathtracing::Camera cam;
  cam.SetPosition(vec3(278,278,-800));
  cam.SetWidth(nx);
  cam.SetHeight(ny);
  cam.SetFieldOfView(40);
  cam.SetUp(vec3(0,1,0));
  cam.SetLookAt(vec3(278,278,-799));
  vtkm::cont::ArrayHandle<unsigned int> seeds = cam.seeds;

  buildBVH(cb);
  for (unsigned int i=0; i<canvasSize; i++){

    unsigned int idx = i;
    auto val = xorshiftWang::getWang32(idx);
    idx++;
    val = xorshiftWang::getWang32(val);
    idx++;
    val = xorshiftWang::getWang32(val);
    idx++;
    val = xorshiftWang::getWang32(val);
    seeds.GetPortalControl().Set(i, val);
  }

  vtkm::worklet::Invoker Invoke;
  for (int s =0; s<ns; s++){
    UVGen uvgen(nx, ny, s);
    Invoke(uvgen, seeds, uvs);

    vtkm::Bounds bounds(vtkm::Range(0,555), vtkm::Range(0,555), vtkm::Range(-800,555));
    cam.CreateRays(rays, bounds);

    MyAlgos::Copy<vtkm::UInt8, vtkm::UInt8, StorageTag>((1UL << 3), rays.Status);


    for (int depth=0; depth<depthcount; depth++){
      MyAlgos::Copy<float, float, StorageTag>(0, sum_values);

      intersect(cb, rays, hrecs,hids, tmin, emitted, attenuation, depth);

      applyMaterials(rays, hrecs, hids, srecs, cb.tex, cb.matType, cb.texType, emitted, seeds, canvasSize, depth);
      generateRays(whichPDF, hrecs, generated_dir, seeds, light_box_pts, light_sphere_pts, light_sphere_radii);
      applyPDFs(cb, rays, hrecs, srecs, sum_values, generated_dir, attenuation, seeds, light_box_pts, light_sphere_pts, light_sphere_radii, lightables, canvasSize, depth);

      vtkm::cont::ArrayHandleCast<vtkm::Int32, vtkm::cont::ArrayHandle<vtkm::UInt8>> castedStatus(rays.Status);

    }

    using CountType = vtkm::cont::ArrayHandleCounting<vtkm::Id>;


    MyAlgos::SliceTransform<
        decltype(emitted),
        decltype(zero),
        decltype(sumtotl),
        decltype(vtkm::Sum())>
        (std::make_tuple((depthcount-1)*canvasSize, (depthcount-1)*canvasSize + canvasSize), emitted,
         std::make_tuple(0, canvasSize), zero,
         std::make_tuple(0, canvasSize), sumtotl, vtkm::Sum());

    for (int depth = depthcount-2; depth >=0; depth--){
      MyAlgos::SliceTransform<
          decltype(attenuation),
          decltype(sumtotl),
          decltype(sumtotl),
          decltype(vtkm::Multiply())>
          (std::make_tuple(depth*canvasSize, depth*canvasSize + canvasSize), attenuation,
           std::make_tuple(0, canvasSize), sumtotl,
           std::make_tuple(0, canvasSize), sumtotl, vtkm::Multiply());

      MyAlgos::SliceTransform<
          decltype(emitted),
          decltype(sumtotl),
          decltype(sumtotl),
          decltype(vtkm::Sum())>
          (std::make_tuple(depth*canvasSize, depth*canvasSize + canvasSize), emitted,
           std::make_tuple(0, canvasSize), sumtotl,
           std::make_tuple(0, canvasSize), sumtotl, vtkm::Sum());


    }

    vtkm::cont::Algorithm::Transform(cols, sumtotl, cols, vtkm::Sum());

    std::cout << "ns: " << s <<" " << vtkm::cont::Algorithm::Reduce(sumtotl, vec3(0.0)) << std::endl;

  }

  std::fstream fs;
  fs.open("output.pnm", std::fstream::out);
  if (fs.is_open()){
    fs << "P3\n" << nx << " "  << ny << " 255" << std::endl;
    for (int i=0; i<cols.GetNumberOfValues(); i++){
      auto col = cols.GetPortalConstControl().Get(i);
      col = de_nan(col);
      col = col / float(ns);
      col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
      int ir = int(255.99*col[0]);
      int ig = int(255.99*col[1]);
      int ib = int(255.99*col[2]);
      fs << ir << " " << ig << " " << ib << std::endl;
    }
    fs.close();
  }
  else
    std::cout << "Couldn't save pnm." << std::endl;
//  std::vector<std::uint8_t> PngBuffer;
}

