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
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/raytracing/Camera.h>
#include <vtkm/cont/internal/DeviceAdapterAlgorithmGeneral.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleCompositeVector.h>
#include <vtkm/rendering/raytracing/Ray.h>
#include <vtkm/rendering/raytracing/RayOperations.h>
#include <omp.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/worklet/Invoker.h>
#include <fstream>
#include "Worklets.h"
#include "SurfaceWorklets.h"
#include "EmitWorklet.h"
#include "ScatterWorklet.h"
#include "PdfWorklet.h"

#include "CornellBox.h"
using ArrayType = vtkm::cont::ArrayHandle<vec3>;
CornellBox cb;


template <class IndexPortalType, class InputPortalType, class OutputPortalType>
struct SliceCopyKernel
{
  IndexPortalType IndexPortal;
  InputPortalType InputPortal;
  OutputPortalType OutputPortal;
  vtkm::Id InputOffset;
  vtkm::Id OutputOffset;

  VTKM_CONT
  SliceCopyKernel(IndexPortalType indexPortal,
             InputPortalType inputPortal,
             OutputPortalType outputPortal,
             vtkm::Id inputOffset = 0,
             vtkm::Id outputOffset = 0)
    : IndexPortal(indexPortal)
    , InputPortal(inputPortal)
    , OutputPortal(outputPortal)
    , InputOffset(inputOffset)
    , OutputOffset(outputOffset)
  {
  }

  VTKM_SUPPRESS_EXEC_WARNINGS
  VTKM_EXEC_CONT
  void operator()(vtkm::Id index) const
  {
    using ValueType = typename OutputPortalType::ValueType;
    auto newIndex = this->IndexPortal.Get(index);
    this->OutputPortal.Set(
      newIndex,
      static_cast<ValueType>(this->InputPortal.Get(newIndex)));
  }

  VTKM_CONT
  void SetErrorMessageBuffer(const vtkm::exec::internal::ErrorMessageBuffer&) {}
};

template <class T, class OutputPortalType>
struct CopyKernel
{
  const T Value;
  OutputPortalType OutputPortal;

  VTKM_CONT
  CopyKernel(T val,
             OutputPortalType outputPortal)
    : Value(val)
    , OutputPortal(outputPortal)
  {
  }

  VTKM_SUPPRESS_EXEC_WARNINGS
  VTKM_EXEC_CONT
  void operator()(vtkm::Id index) const
  {
    this->OutputPortal.Set(
      index,
      this->Value);
  }

  VTKM_CONT
  void SetErrorMessageBuffer(const vtkm::exec::internal::ErrorMessageBuffer&) {}
};
template <typename IndexType1,
          typename InPortalType1,
          typename IndexType2,
          typename InPortalType2,
          typename IndexOutType,
          typename OutPortalType,
          typename BinaryFunctor>
struct BinarySliceTransformKernel : vtkm::exec::FunctorBase
{
  IndexType1     Index1;
  IndexType2     Index2;
  IndexOutType     OutIndex;
  InPortalType1 InPortal1;
  InPortalType2 InPortal2;
  OutPortalType OutPortal;
  BinaryFunctor BinaryOperator;

  VTKM_CONT
  BinarySliceTransformKernel(const IndexType1 &index1,
                        const InPortalType1& inPortal1,
                             const IndexType2 &index2,
                        const InPortalType2& inPortal2,
                             const IndexOutType &outIndex,
                        const OutPortalType& outPortal,
                        BinaryFunctor binaryOperator)
    : Index1(index1)
    , Index2(index2)
    , InPortal1(inPortal1)
    , InPortal2(inPortal2)
    , OutIndex(outIndex)
    , OutPortal(outPortal)
    , BinaryOperator(binaryOperator)
  {
  }

  VTKM_SUPPRESS_EXEC_WARNINGS
  VTKM_EXEC
  void operator()(vtkm::Id index) const
  {
    auto newIndex1 = this->Index1.Get(index);
    auto newIndex2 = this->Index2.Get(index);
    auto outIndex = this->OutIndex.Get(index);
    this->OutPortal.Set(
      outIndex, this->BinaryOperator(this->InPortal1.Get(newIndex1), this->InPortal2.Get(newIndex2)));
  }
};
template <
          typename InPortalType1,
          typename InPortalType2,
          typename OutPortalType,
          typename BinaryFunctor>
struct BinarySliceTransformIdxKernel : vtkm::exec::FunctorBase
{
  InPortalType1 InPortal1;
  InPortalType2 InPortal2;
  OutPortalType OutPortal;
  BinaryFunctor BinaryOperator;
  vtkm::Id idxBegin1, idxBegin2, idxEnd1, idxEnd2, idxOutBegin, idxOutEnd;
  VTKM_CONT
  BinarySliceTransformIdxKernel(const std::tuple<vtkm::Id, vtkm::Id> idx1,
                        const InPortalType1& inPortal1,
                             const std::tuple<vtkm::Id, vtkm::Id> idx2,
                        const InPortalType2& inPortal2,
                             const std::tuple<vtkm::Id, vtkm::Id> outIdx,
                        const OutPortalType& outPortal,
                        BinaryFunctor binaryOperator)
    : idxBegin1(std::get<0>(idx1))
    , idxBegin2(std::get<0>(idx2))
    , idxOutBegin(std::get<0>(outIdx))
    , idxEnd1(std::get<1>(idx1))
    , idxEnd2(std::get<1>(idx2))
    , idxOutEnd(std::get<1>(outIdx))
    , InPortal1(inPortal1)
    , InPortal2(inPortal2)
    , OutPortal(outPortal)
    , BinaryOperator(binaryOperator)
  {
  }

  VTKM_SUPPRESS_EXEC_WARNINGS
  VTKM_EXEC
  void operator()(vtkm::Id index) const
  {
    auto newIndex1 = idxBegin1 + index;
    auto newIndex2 = idxBegin2 + index;
    auto outIndex = idxOutBegin + index;
    this->OutPortal.Set(
      outIndex, this->BinaryOperator(this->InPortal1.Get(newIndex1), this->InPortal2.Get(newIndex2)));
  }
};

template <class DerivedAlgorithm, class DeviceAdapterTag>
struct MyAlgorithms
{
private:
  template <typename T, class CIn>
  VTKM_CONT static T GetExecutionValue(const vtkm::cont::ArrayHandle<T, CIn>& input, vtkm::Id index)
  {
    using OutputArrayType = vtkm::cont::ArrayHandle<T, vtkm::cont::StorageTagBasic>;

    OutputArrayType output;
    auto inputPortal = input.PrepareForInput(DeviceAdapterTag());
    auto outputPortal = output.PrepareForOutput(1, DeviceAdapterTag());

    vtkm::cont::internal::CopyKernel<decltype(inputPortal), decltype(outputPortal)> kernel(
      inputPortal, outputPortal, index);

    DerivedAlgorithm::Schedule(kernel, 1);

    return output.GetPortalConstControl().Get(0);
  }
public:
  template <typename IndexType1,
            typename InputType1,
            typename IndexType2,
            typename InputType2,
            typename OutIndexType,
            typename OutputType,
            typename BinaryFunctorType>
  VTKM_CONT static void SliceTransform(
                            const IndexType1 &index1,
                            const InputType1& input1,
                            const IndexType2 &index2,
                            const InputType2& input2,
                            const OutIndexType &outIndex,
                             OutputType& output,
                            BinaryFunctorType binaryOperator)
  {
    const vtkm::Id inSize = index1.GetNumberOfValues();
    auto index1Portal  = index1.PrepareForInput(DeviceAdapterTag());
    auto index2Portal  = index2.PrepareForInput(DeviceAdapterTag());
    auto input1Portal = input1.PrepareForInput(DeviceAdapterTag());
    auto input2Portal = input2.PrepareForInput(DeviceAdapterTag());
    auto outIndexPortal = outIndex.PrepareForInput(DeviceAdapterTag());
    auto outputPortal = output.PrepareForOutput(inSize, DeviceAdapterTag());

    BinarySliceTransformKernel<decltype(index1Portal),
                               decltype(input1Portal),
                                decltype(index2Portal),
                               decltype(input2Portal),
                               decltype(outIndexPortal),
                               decltype(outputPortal),
                               BinaryFunctorType> kernel(index1Portal, input1Portal, index2Portal, input2Portal, outIndexPortal, outputPortal, binaryOperator);
    DerivedAlgorithm::Schedule(kernel, inSize);
  }
  template <
            typename InputType1,
            typename InputType2,
            typename OutputType,
            typename BinaryFunctorType>
  VTKM_CONT static void SliceTransform(
                            std::tuple<vtkm::Id, vtkm::Id> idx1,
                            const InputType1& input1,
                            std::tuple<vtkm::Id, vtkm::Id> idx2,
                            const InputType2& input2,
                            std::tuple<vtkm::Id, vtkm::Id> outIdx,
                             OutputType& output,
                            BinaryFunctorType binaryOperator)
  {
    const vtkm::Id inSize = output.GetNumberOfValues();
    auto input1Portal  = input1.PrepareForInput(DeviceAdapterTag());
    auto input2Portal = input2.PrepareForInput(DeviceAdapterTag());
    auto outputPortal = output.PrepareForOutput(inSize, DeviceAdapterTag());

    BinarySliceTransformIdxKernel<
                               decltype(input1Portal),
                               decltype(input2Portal),
                               decltype(outputPortal),
                               BinaryFunctorType> kernel(idx1, input1Portal, idx2, input2Portal, outIdx, outputPortal, binaryOperator);
    DerivedAlgorithm::Schedule(kernel, inSize);
  }

  template <typename T, typename U, class CIn, class COut>
  VTKM_CONT static void SliceCopy(
                            const vtkm::cont::ArrayHandle<vtkm::Id> &index,
                            const vtkm::cont::ArrayHandle<T, CIn>& input,
                             vtkm::cont::ArrayHandle<U, COut>& output)
  {
    const vtkm::Id inSize = index.GetNumberOfValues();
    auto indexPortal  = index.PrepareForInput(DeviceAdapterTag());
    auto inputPortal = input.PrepareForInput(DeviceAdapterTag());
    auto outputPortal = output.PrepareForOutput(inSize, DeviceAdapterTag());

    SliceCopyKernel<decltype(indexPortal),
                             decltype(inputPortal),
                             decltype(outputPortal)> kernel(indexPortal, inputPortal, outputPortal);
    DerivedAlgorithm::Schedule(kernel, inSize);
  }

  template <typename T, typename U, class COut>
  VTKM_CONT static void Copy(
                            T input,
                             vtkm::cont::ArrayHandle<U, COut>& output)
  {
    const vtkm::Id inSize = output.GetNumberOfValues();
    auto outputPortal = output.PrepareForOutput(inSize, DeviceAdapterTag());

    CopyKernel<T, decltype(outputPortal)> kernel(input, outputPortal);
    DerivedAlgorithm::Schedule(kernel, inSize);
  }
};
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


std::vector<vtkm::cont::ArrayHandle<vec3>> ptsArray;


using RayType = vtkm::rendering::raytracing::Ray<float>;

template<typename HitRecord,
         typename HitId>
void intersect(RayType &rays,
               HitRecord &hrecs,
               HitId &hids,
               vtkm::cont::ArrayHandle<float> &tmin,
               ArrayType &emitted,
               ArrayType &attenuation,
               const vtkm::Id depth)
{
  using MyAlgos = MyAlgorithms<vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>, VTKM_DEFAULT_DEVICE_ADAPTER_TAG>;
  using StorageTag = vtkm::cont::StorageTagBasic;

  MyAlgos::Copy<float, float, StorageTag>(std::numeric_limits<float>::max(), rays.Distance);
  MyAlgos::Copy<float, float, StorageTag>(0.001, tmin);
  vtkm::worklet::Invoker Invoke;

  vtkm::Id canvasSize = rays.DirX.GetNumberOfValues();
//  for (int i=0; i<cb.cellTypeArray.size(); i++){
//    if (cb.cellTypeArray[i] == 0){

//    }
//    else if (cb.cellTypeArray[i] == 1){
//      //xz
//      QuadIntersect quad;

//      vtkm::cont::ArrayHandle<vtkm::Int32> nodes;
//      nodes.Allocate(rays.Dir.GetNumberOfValues());
//      for (int j=0; j<nodes.GetNumberOfValues(); j++){
//        nodes.GetPortalControl().Set(j, 0);
//      }

//       vtkm::cont::ArrayHandle<vtkm::Id> leafs;
//      leafs.Allocate(6);
//      leafs.GetPortalControl().Set(0, 5);
//      leafs.GetPortalControl().Set(1, 0);
//      leafs.GetPortalControl().Set(2, 1);
//      leafs.GetPortalControl().Set(3, 2);
//      leafs.GetPortalControl().Set(4, 3);
//      leafs.GetPortalControl().Set(5, 4);

//      vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>> QuadIds;
//      QuadIds.Allocate(5);
//      QuadIds.GetPortalControl().Set(0, vtkm::Vec<vtkm::Id, 5>(0,0,1,2,3));
//      QuadIds.GetPortalControl().Set(1, vtkm::Vec<vtkm::Id, 5>(0,4,5,6,7));
//      QuadIds.GetPortalControl().Set(2, vtkm::Vec<vtkm::Id, 5>(0,8,9,10,11));
//      QuadIds.GetPortalControl().Set(3, vtkm::Vec<vtkm::Id, 5>(0,12,13,14,15));
//      QuadIds.GetPortalControl().Set(4, vtkm::Vec<vtkm::Id, 5>(0,16,17,18,19));
//      Invoke(quad, nodes, rays.Origin, rays.Dir, hrecs, hids, tmin, rays.Distance, rays.Status, cb.pts1[i], leafs,
//             cb.matIdx[i], cb.texIdx[i], QuadIds);

//    }
//    else if (cb.cellTypeArray[i] == 2){
//      //xy
//      QuadIntersect quad;

//      vtkm::cont::ArrayHandle<vtkm::Int32> nodes;
//      nodes.Allocate(rays.Dir.GetNumberOfValues());
//      for (int j=0; j<nodes.GetNumberOfValues(); j++){
//        nodes.GetPortalControl().Set(j, 0);
//      }

//       vtkm::cont::ArrayHandle<vtkm::Id> leafs;
//      leafs.Allocate(4);
//      leafs.GetPortalControl().Set(0, 3);
//      leafs.GetPortalControl().Set(1, 0);
//      leafs.GetPortalControl().Set(2, 1);
//      leafs.GetPortalControl().Set(3, 2);

//      vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,5>> QuadIds;
//      QuadIds.Allocate(3);
//      QuadIds.GetPortalControl().Set(0, vtkm::Vec<vtkm::Id, 5>(0,0,1,2,3));
//      QuadIds.GetPortalControl().Set(1, vtkm::Vec<vtkm::Id, 5>(0,4,5,6,7));
//      QuadIds.GetPortalControl().Set(2, vtkm::Vec<vtkm::Id, 5>(0,8,9,10,11));

//      Invoke(quad, nodes, rays.Origin, rays.Dir, hrecs, hids, tmin, rays.Distance, rays.Status, cb.pts1[i], leafs,
//             cb.matIdx[i], cb.texIdx[i], QuadIds);

//    }
//    else if (cb.cellTypeArray[i] == 3){

//      vtkm::cont::ArrayHandle<vtkm::Int32> nodes;
//      nodes.Allocate(rays.Dir.GetNumberOfValues());
//      for (int j=0; j<nodes.GetNumberOfValues(); j++){
//        nodes.GetPortalControl().Set(j, 0);
//      }

//      vtkm::cont::ArrayHandle<vtkm::Id> leafs;
//      leafs.Allocate(2);
//      leafs.GetPortalControl().Set(0, 1);
//      leafs.GetPortalControl().Set(1, 0);

//      vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id,3>> SphereIds;
//      SphereIds.Allocate(1);
//      SphereIds.GetPortalControl().Set(0, vtkm::Vec<vtkm::Id, 3>(0,0,1));

//      SphereExecWrapper surf;
//      SphereIntersecttWorklet sphereIntersect(canvasSize, depth);
//      Invoke(sphereIntersect, nodes, rays.Origin, rays.Dir, hrecs, hids, tmin, rays.Distance, rays.Status, surf, cb.pts1[i], leafs,
//             cb.matIdx[i], cb.texIdx[i], SphereIds);

//    }
//  }

    QuadIntersect quad;

    vtkm::cont::ArrayHandle<vtkm::Int32> nodes;
    nodes.Allocate(rays.Dir.GetNumberOfValues());
    for (int j=0; j<nodes.GetNumberOfValues(); j++){
      nodes.GetPortalControl().Set(j, 0);
    }

     vtkm::cont::ArrayHandle<vtkm::Id> leafs;
    leafs.Allocate(13);
    leafs.GetPortalControl().Set(0, 12);
    for (int i=0; i<12; i++){

      leafs.GetPortalControl().Set(i+1, i);
    }

    Invoke(quad, nodes, rays.Origin, rays.Dir, hrecs, hids, tmin, rays.Distance, rays.Status, cb.pts1, leafs,
           cb.matIdx, cb.texIdx, cb.QuadIds);

  CollectIntersecttWorklet collectIntersect(canvasSize, depth);
  Invoke(collectIntersect, rays.Status, emitted, attenuation);

}

template<typename HitRecord, typename HitId>
void applyMaterials(RayType &rays,
                    HitRecord &hrecs,
                    HitId &hids,
                    vtkm::cont::ArrayHandle<ScatterRecord> &srecs,
                    vtkm::cont::ArrayHandle<vec3> tex,
                    vtkm::cont::ArrayHandle<int> matType,
                    vtkm::cont::ArrayHandle<int> texType,
                    ArrayType &emitted,
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

template<typename HitRecord>
void generateRays(
    vtkm::cont::ArrayHandle<int> &whichPDF,
    HitRecord &hrecs,
    ArrayType &generated_dir,
    vtkm::cont::ArrayHandle<unsigned int> &seeds,
    std::vector<ArrayType> &light_box_pts,
    std::vector<ArrayType> &light_sphere_pts
    )
{
  vtkm::worklet::Invoker Invoke;
  GenerateDir genDir(3);
  CosineGenerateDir cosGenDir(1);
  XZRectGenerateDir xzRectGenDir(2);
  SphereGenerateDir sphereGenDir(3);

  Invoke(genDir, seeds, whichPDF);
  Invoke(cosGenDir, whichPDF, hrecs, generated_dir, seeds,  light_sphere_pts[0], light_sphere_pts[1]);
  Invoke(xzRectGenDir, whichPDF, hrecs, generated_dir, seeds, light_box_pts[0], light_box_pts[1]);
  Invoke(sphereGenDir, whichPDF, hrecs, generated_dir, seeds, light_sphere_pts[0], light_sphere_pts[1]);
  }

template<typename HitRecord>
void applyPDFs(
    RayType &rays,
    HitRecord &hrecs,
    vtkm::cont::ArrayHandle<ScatterRecord> srecs,
    vtkm::cont::ArrayHandle<vtkm::Float32> &sum_values,
    ArrayType generated_dir,
    ArrayType &attenuation,
    vtkm::cont::ArrayHandle<unsigned int> &seeds,
    std::vector<ArrayType> light_box_pts,
    std::vector<ArrayType> light_sphere_pts,
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
  SphereExecWrapper surf;

  Invoke(spherePDFWorklet, rays.Origin, rays.Dir,hrecs, rays.Status, sum_values, generated_dir, seeds, surf, light_sphere_pts[0], light_sphere_pts[1]);

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
  using MyAlgos = MyAlgorithms<vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>, VTKM_DEFAULT_DEVICE_ADAPTER_TAG>;
  using StorageTag = vtkm::cont::StorageTagBasic;
  using Device = VTKM_DEFAULT_DEVICE_ADAPTER_TAG;

  const auto tup = parse(argc, argv);
  const int nx = std::get<0>(tup);
  const int ny = std::get<1>(tup);
  const int ns = std::get<2>(tup);
  const int depthcount = std::get<3>(tup);

  auto canvasSize = nx*ny;

  constexpr int lightables = 2;
  std::vector<ArrayType> light_box_pts(2), light_sphere_pts(2);
  light_box_pts[0].Allocate(1);
  light_box_pts[0].GetPortalControl().Set(0, vec3(213, 554, 227));
  light_box_pts[1].Allocate(1);
  light_box_pts[1].GetPortalControl().Set(0, vec3(343, 554, 332));

  light_sphere_pts[0].Allocate(1);
  light_sphere_pts[0].GetPortalControl().Set(0, vec3(190, 90, 190));
  light_sphere_pts[1].Allocate(1);
  light_sphere_pts[1].GetPortalControl().Set(0, vec3(90,0,0));

  cb.build();

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 2>> uvs;
  uvs.Allocate(nx*ny);

  vtkm::cont::ArrayHandle<vtkm::Id> rays_index;

  vtkm::cont::ArrayHandle<unsigned int> seeds;
  seeds.Allocate(canvasSize);
  vtkm::rendering::raytracing::Ray<float> rays;
  rays.EnableIntersectionData();
  vtkm::rendering::raytracing::RayOperations::Resize(rays, canvasSize, Device());


  MyAlgos::Copy<vtkm::UInt8, vtkm::UInt8, StorageTag>((1UL << 3), rays.Status);

  vtkm::cont::ArrayHandleConstant<vec3> zero(vec3(0.0f), nx*ny);
  vtkm::cont::ArrayHandle<vec3> cols;
  cols.Allocate(nx*ny);
  MyAlgos::Copy<vec3, vec3, StorageTag>(vec3(0.0f), cols);

  vtkm::cont::ArrayHandle<vtkm::Float32> sum_values;

  vtkm::cont::ArrayHandle<ScatterRecord> srecs;

  vtkm::cont::ArrayHandle<vtkm::Int32> matIdArray, texIdArray;
  matIdArray.Allocate(canvasSize);
  texIdArray.Allocate(canvasSize);

  vtkm::cont::ArrayHandle<vtkm::Float32> pxArray,pyArray,pzArray;
  pxArray.Allocate(rays.U.GetNumberOfValues());
  pyArray.Allocate(rays.U.GetNumberOfValues());
  pzArray.Allocate(rays.U.GetNumberOfValues());

  using HitRecord = vtkm::cont::ArrayHandleCompositeVector<decltype(rays.U),
  decltype(rays.V),
  decltype(rays.Distance),
  decltype(rays.NormalX),
  decltype(rays.NormalY),
  decltype(rays.NormalZ),
  decltype(pxArray),
  decltype(pyArray),
  decltype(pzArray)>;

  auto hrecs = HitRecord(rays.U, rays.V, rays.Distance, rays.NormalX, rays.NormalY, rays.NormalZ, pxArray, pyArray, pzArray);

  using HitId = vtkm::cont::ArrayHandleCompositeVector<decltype(matIdArray), decltype(texIdArray)>;
  auto hids = HitId(matIdArray, texIdArray);

  vtkm::cont::ArrayHandle<float> tmin;
  ArrayType attenuation;
  ArrayType emitted;
  ArrayType generated_dir;
  vtkm::cont::ArrayHandle<int> whichPDF;

  attenuation.Allocate(canvasSize* depthcount);
  emitted.Allocate(canvasSize * depthcount);
  srecs.Allocate(canvasSize);
  ArrayType sumtotl;
  sumtotl.Allocate(canvasSize);
  tmin.Allocate(canvasSize);
  sum_values.Allocate(nx*ny);
  generated_dir.Allocate(nx*ny);
  whichPDF.Allocate(nx*ny);

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
    RayGen raygen(nx,ny, 40,40,
                  vtkm::Vec<vtkm::Float32,3>(0,0,1),
                  vtkm::Vec<vtkm::Float32,3>(0,1,0),
                  0, nx, 0, 0, ns);
    Invoke(raygen, rays.DirX, rays.DirY, rays.DirZ, seeds, rays.PixelIdx);

    vec3 lookfrom(278, 278, -800);

    MyAlgos::Copy<float, float, StorageTag>(lookfrom[0], rays.OriginX);
    MyAlgos::Copy<float, float, StorageTag>(lookfrom[1], rays.OriginY);
    MyAlgos::Copy<float, float, StorageTag>(lookfrom[2], rays.OriginZ);


    MyAlgos::Copy<vtkm::UInt8, vtkm::UInt8, StorageTag>((1UL << 3), rays.Status);


    for (int depth=0; depth<depthcount; depth++){
      MyAlgos::Copy<float, float, StorageTag>(0, sum_values);


      intersect(rays, hrecs,hids, tmin, emitted, attenuation, depth);

      applyMaterials(rays, hrecs, hids, srecs, cb.tex, cb.matType, cb.texType, emitted, seeds, canvasSize, depth);
      generateRays(whichPDF, hrecs, generated_dir, seeds, light_box_pts, light_sphere_pts);
      applyPDFs(rays, hrecs, srecs, sum_values, generated_dir, attenuation, seeds, light_box_pts, light_sphere_pts, lightables, canvasSize, depth);

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

