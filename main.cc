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
#include <sstream>
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
#include <vtkm/cont/ArrayHandleExtractComponent.h>

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
#include "QuadIntersector.h"
#include "SphereIntersector.h"

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
               vtkm::cont::ArrayHandle<vtkm::Int32> &matIdArray,
               vtkm::cont::ArrayHandle<vtkm::Int32> &texIdArray,
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


  vtkm::rendering::pathtracing::QuadIntersector quadIntersector;

  quadIntersector.SetData(cb.coord, cb.QuadIds, cb.matIdx[0], cb.texIdx[0], matIdArray, texIdArray);
  quadIntersector.IntersectRays(rays);

  vtkm::rendering::pathtracing::SphereIntersector sphereIntersector;

  sphereIntersector.SetData(cb.coord, cb.SphereIds, cb.SphereRadii, cb.matIdx[1], cb.texIdx[1], matIdArray, texIdArray);
  sphereIntersector.IntersectRays(rays);


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
void generateRays(CornellBox &cb,
                  vtkm::cont::ArrayHandle<int> &whichPDF,
                  HitRecord &hrecs,
                  GenDirType &generated_dir,
                  vtkm::cont::ArrayHandle<unsigned int> &seeds,
                  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>>  light_box_pointids,
                  vtkm::cont::ArrayHandle<vtkm::Id> light_box_indices,
                  vtkm::cont::ArrayHandle<vtkm::Id> &light_sphere_pointids,
                  vtkm::cont::ArrayHandle<vtkm::Id> & light_sphere_indices
                  )
{
  vtkm::worklet::Invoker Invoke;
  GenerateDir genDir(3);
  CosineGenerateDir cosGenDir(1);
  QuadGenerateDir quadGenDir(2);
  SphereGenerateDir sphereGenDir(3);

  Invoke(genDir, seeds, whichPDF);
  Invoke(cosGenDir, whichPDF, hrecs, generated_dir, seeds);
  Invoke(quadGenDir, whichPDF, hrecs, generated_dir, seeds, light_box_pointids, light_box_indices, cb.coord);
  Invoke(sphereGenDir, whichPDF, hrecs, generated_dir, seeds, light_sphere_pointids, light_sphere_indices, cb.coord, cb.SphereRadii);
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
               vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>>  light_box_pointids,
               vtkm::cont::ArrayHandle<vtkm::Id> light_box_indices,
              vtkm::cont::ArrayHandle<vtkm::Id> &light_sphere_pointids,
              vtkm::cont::ArrayHandle<vtkm::Id> & light_sphere_indices,
              int lightables,
              vtkm::Id canvasSize,
              vtkm::Id depth
    )
{
  vtkm::worklet::Invoker Invoke;
  QuadPDFWorklet quadPDFWorklet(lightables);
  SpherePDFWorklet spherePDFWorklet(lightables);
  PDFCosineWorklet pdfWorklet(canvasSize, depth, canvasSize, lightables);
  QuadExecWrapper quadSurf(cb.QuadIds, cb.matIdx[0], cb.texIdx[0]);
  Invoke(quadPDFWorklet, rays.Origin, rays.Dir,hrecs,
         rays.Status, sum_values, generated_dir, seeds,quadSurf,
          light_box_pointids, light_box_indices, cb.coord);
  SphereExecWrapper surf(cb.SphereIds, cb.SphereRadii, cb.matIdx[1], cb.texIdx[1]);


  Invoke(spherePDFWorklet,
         rays.Origin,
         rays.Dir,hrecs,
         rays.Status,
         sum_values,
         generated_dir,
         seeds,
         surf,
         light_sphere_pointids,
         light_sphere_indices,
         cb.coord,
         cb.SphereRadii);

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

  bool hemi = false;
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
    else if (!strcmp(argv[i], "-hemisphere"))
    {
      hemi = true;
    }
  }

  return std::make_tuple(x,y, s, depth, hemi);
}

ArrayType run(int nx, int ny, int samplecount, int depthcount,
              vtkm::rendering::Canvas &canvas, vtkm::rendering::pathtracing::Camera &rayCam)
{
  using MyAlgos = details::PathAlgorithms<vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>, VTKM_DEFAULT_DEVICE_ADAPTER_TAG>;
  using StorageTag = vtkm::cont::StorageTagBasic;
  using Device = VTKM_DEFAULT_DEVICE_ADAPTER_TAG;

  CornellBox cb;

  auto canvasSize = nx*ny;

  constexpr int lightables = 2;
  std::vector<ArrayType> light_sphere_pts(1);
  std::vector<vtkm::cont::ArrayHandle<vtkm::Float32>> light_sphere_radii(1);
  vtkm::Vec<vtkm::Id, 5> tmp[] = {vtkm::Vec<vtkm::Id,5>(0,8,9,10,11)};
  auto light_box_pointids = vtkm::cont::make_ArrayHandle(tmp,1);
  vtkm::Id tmp2[] = {0};
  auto light_box_indices = vtkm::cont::make_ArrayHandle(tmp2,1);
  vtkm::Id sphere_tmp[] = {vtkm::Id(4*12)};
  auto light_sphere_pointids = vtkm::cont::make_ArrayHandle(sphere_tmp,1);
  vtkm::Id sphere_tmp2[] = {0};
  auto light_sphere_indices = vtkm::cont::make_ArrayHandle(sphere_tmp2,1);

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


  rays.AddBuffer(1, "specular_Ox");
  rays.AddBuffer(1, "specular_Oy");
  rays.AddBuffer(1, "specular_Oz");
  rays.AddBuffer(1, "specular_Dx");
  rays.AddBuffer(1, "specular_Dy");
  rays.AddBuffer(1, "specular_Dz");
  rays.AddBuffer(1, "specular_Ax");
  rays.AddBuffer(1, "specular_Ay");
  rays.AddBuffer(1, "specular_Az");

  auto srecs = vtkm::rendering::pathtracing::QuadIntersector::ScatterRecord(
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
  vtkm::cont::ArrayHandle<vtkm::Int32> matIdArray, texIdArray;
  matIdArray.Allocate(canvasSize);
  texIdArray.Allocate(canvasSize);

  sum_values = rays.GetBuffer("sum_values").Buffer;
  vtkm::Bounds bounds(vtkm::Range(0,555), vtkm::Range(0,555), vtkm::Range(0,555));
  whichPDF.Allocate(nx*ny);
  vtkm::cont::ArrayHandle<unsigned int> seeds = rayCam.seeds;

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
  for (int s =0; s<samplecount; s++){
    UVGen uvgen(nx, ny, s);
    Invoke(uvgen, seeds, uvs);

    rayCam.CreateRays(rays, bounds);

    MyAlgos::Copy<vtkm::UInt8, vtkm::UInt8, StorageTag>((1UL << 3), rays.Status);


    for (int depth=0; depth<depthcount; depth++){
      MyAlgos::Copy<float, float, StorageTag>(0, sum_values);
      auto hrecs = vtkm::rendering::pathtracing::QuadIntersector::HitRecord(rays.U, rays.V, rays.Distance, rays.NormalX, rays.NormalY, rays.NormalZ, rays.IntersectionX, rays.IntersectionY, rays.IntersectionZ);
      auto hids = vtkm::rendering::pathtracing::QuadIntersector::HitId(matIdArray, texIdArray);

      intersect(cb, rays, hrecs,hids, matIdArray, texIdArray, tmin, emitted, attenuation, depth);

      applyMaterials(rays, hrecs, hids, srecs, cb.tex, cb.matType, cb.texType, emitted, seeds, canvasSize, depth);
      generateRays(cb, whichPDF, hrecs, generated_dir, seeds, light_box_pointids,light_box_indices, light_sphere_pointids, light_sphere_indices);

      applyPDFs(cb, rays, hrecs, srecs, sum_values, generated_dir, attenuation, seeds,
                light_box_pointids, light_box_indices, light_sphere_pointids, light_sphere_indices, lightables, canvasSize, depth);

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

  return cols;
}

void save(std::string fn,
          int nx, int ny, int samplecount,
          vtkm::cont::ArrayHandle<vec3> &cols)
{
  std::fstream fs;
  fs.open(fn.c_str(), std::fstream::out);
  if (fs.is_open()){
    fs << "P3\n" << nx << " "  << ny << " 255" << std::endl;
    for (int i=0; i<cols.GetNumberOfValues(); i++){
      auto col = cols.GetPortalConstControl().Get(i);
      col = de_nan(col);
      col = col / float(samplecount);
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

void generateHemisphere(int nx, int ny, int samplecount, int depthcount)
{
  vtkm::rendering::CanvasRayTracer canvas(nx,ny);
  vtkm::rendering::pathtracing::Camera rayCam;
  rayCam.SetPosition(vec3(278,278,-800));
  rayCam.SetWidth(nx);
  rayCam.SetHeight(ny);
  rayCam.SetFieldOfView(40);
  rayCam.SetUp(vec3(0,1,0));
  rayCam.SetLookAt(vec3(278,278,278));

  int numPhi = 30;
  int numTheta = 30;

  float rTheta = (2.0*M_PI)/float(numTheta);
  float rPhi = (M_PI/2.0)/float(numPhi);

  float r = -1078;
  for (int i=0; i<numTheta; i++){
    for (int j=0; j<numPhi; j++){
      auto x = r * cos(i * rTheta) * sin(j * rPhi);
      auto y = r * sin(i * rTheta) * sin(j * rPhi);
      auto z = r * cos(i*rPhi);

      vec3 pos(x+278, y+278, z+278 );
      rayCam.SetPosition(pos);
      auto cols = run(nx,ny, samplecount, depthcount, canvas, rayCam);
      std::stringstream sstr;
      sstr << "output-" << i*rTheta << "-" << j*rPhi << ".pnm";
      save(sstr.str(), nx, ny, samplecount, cols);
    }
  }
}
int main(int argc, char *argv[]) {

  const auto tup = parse(argc, argv);
  const int nx = std::get<0>(tup);
  const int ny = std::get<1>(tup);
  const int samplecount = std::get<2>(tup);
  const int depthcount = std::get<3>(tup);
  const bool hemi = std::get<4>(tup);
  if (hemi)
    generateHemisphere(nx,ny, samplecount, depthcount);
  else
  {
    vtkm::rendering::CanvasRayTracer canvas(nx,ny);
    vtkm::rendering::pathtracing::Camera rayCam;
    rayCam.SetPosition(vec3(278,278,-800));
    rayCam.SetWidth(nx);
    rayCam.SetHeight(ny);
    rayCam.SetFieldOfView(40);
    rayCam.SetUp(vec3(0,1,0));
    rayCam.SetLookAt(vec3(278,278,278));
    auto cols = run(nx,ny, samplecount, depthcount, canvas, rayCam);
    std::stringstream sstr;
    sstr << "output.pnm";
    save(sstr.str(), nx, ny, samplecount, cols);

  }
}

