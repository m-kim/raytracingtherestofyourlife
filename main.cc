//==================================================================================================
// Written in 2016 by Peter Shirley <ptrshrl@gmail.com>
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
#include <JSONPNGConvert.h>
#include <lodepng.h>
#include <vtkm/cont/ArrayHandleConstant.h>
#include <vtkm/rendering/raytracing/Camera.h>
#include <vtkm/cont/internal/DeviceAdapterAlgorithmGeneral.h>
#include <omp.h>
#include <vtkm/cont/Algorithm.h>
#include <fstream>
#include "Worklets.h"
#include "SurfaceWorklets.h"
#include "sphere.h"
#include "moving_sphere.h"
#include "hitable_list.h"
#include "float.h"
#include "camera.h"
#include "material.h"
#include "bvh.h"
#include "box.h"
#include "surface_texture.h"
#include "aarect.h"
#include "texture.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "pdf.h"
#include "EmitWorklet.h"
#include "ScatterWorklet.h"


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

vtkm::cont::ArrayHandle<vec3> tex;
vtkm::cont::ArrayHandle<vtkm::Id> matIdx[5];
vtkm::cont::ArrayHandle<vtkm::Id> texIdx[5];
vtkm::cont::ArrayHandle<int> matType, texType, flipped[5];
vtkm::cont::ArrayHandle<vec3> pts1[5], pts2[5];
vtkm::cont::ArrayHandle<vec3> norms;

std::vector<vtkm::cont::ArrayHandle<vec3>> ptsArray;
std::vector<int> cellTypeArray;

void cornell_box(hitable **scene, camera **cam, float aspect) {
  tex.Allocate(4);
  tex.GetPortalControl().Set(0, vec3(0.65, 0.05, 0.05));
  tex.GetPortalControl().Set(1, vec3(0.73, 0.73, 0.73));
  tex.GetPortalControl().Set(2, vec3(0.12, 0.45, 0.15));
  tex.GetPortalControl().Set(3, vec3(15, 15, 15));

  matType.Allocate(5);
  matType.GetPortalControl().Set(0, 0); //lambertian
  matType.GetPortalControl().Set(1, 0); //lambertian
  matType.GetPortalControl().Set(2, 0); //lambertian
  matType.GetPortalControl().Set(3, 1); //light
  matType.GetPortalControl().Set(4, 2); //dielectric

  texType.Allocate(5);
  texType.GetPortalControl().Set(0, 0); //red
  texType.GetPortalControl().Set(1, 1); //white
  texType.GetPortalControl().Set(2, 2); //green
  texType.GetPortalControl().Set(3, 3); //super bright
  texType.GetPortalControl().Set(4, 0); //dielectric


  cellTypeArray.resize(5);
  norms.Allocate(16);
    int i = 0;
    hitable **list = new hitable*[8];

    flipped[0].Allocate(2);
    matIdx[0].Allocate(2);
    texIdx[0].Allocate(2);
    pts1[0].Allocate(2);
    pts2[0].Allocate(2);
    cellTypeArray[0] = 0;
    matIdx[0].GetPortalControl().Set(0, 2);
    texIdx[0].GetPortalControl().Set(0, 2);
    pts1[0].GetPortalControl().Set(0, vec3(555,0,0));
    pts2[0].GetPortalControl().Set(0, vec3(555,555,555));
    flipped[0].GetPortalControl().Set(0, 1);
    list[i++] = new flip_normals(new yz_rect(0, 555, 0, 555, 555, 2,2));


    matIdx[0].GetPortalControl().Set(1, 0);
    texIdx[0].GetPortalControl().Set(1, 0);
    pts1[0].GetPortalControl().Set(1, vec3(0,0,0));
    pts2[0].GetPortalControl().Set(1, vec3(0,555,555));
    flipped[0].GetPortalControl().Set(1, 0);
    list[i++] = new yz_rect(0, 555, 0, 555, 0, 0,0);

    flipped[1].Allocate(3);
    matIdx[1].Allocate(3);
    texIdx[1].Allocate(3);
    pts1[1].Allocate(3);
    pts2[1].Allocate(3);
    cellTypeArray[1] =  1;
    matIdx[1].GetPortalControl().Set(0, 3);
    texIdx[1].GetPortalControl().Set(0, 3);
    pts1[1].GetPortalControl().Set(0, vec3(213,554,227));
    pts2[1].GetPortalControl().Set(0, vec3(343,554,332));
    flipped[1].GetPortalControl().Set(0, 1);
    list[i++] = new flip_normals(new xz_rect(213, 343, 227, 332, 554, 3,3));

    matIdx[1].GetPortalControl().Set(1, 1);
    texIdx[1].GetPortalControl().Set(1, 1);
    pts1[1].GetPortalControl().Set(1, vec3(0,555,0));
    pts2[1].GetPortalControl().Set(1, vec3(555,555,555));
    flipped[1].GetPortalControl().Set(1, 1);
    list[i++] = new flip_normals(new xz_rect(0, 555, 0, 555, 555, 1,1));

    matIdx[1].GetPortalControl().Set(2, 1);
    texIdx[1].GetPortalControl().Set(2, 1);
    pts1[1].GetPortalControl().Set(2, vec3(0,0,0));
    pts2[1].GetPortalControl().Set(2, vec3(555,0,555));
    flipped[1].GetPortalControl().Set(2, 0);
    list[i++] = new xz_rect(0, 555, 0, 555, 0, 1,1);

    flipped[2].Allocate(1);
    matIdx[2].Allocate(1);
    texIdx[2].Allocate(1);
    pts1[2].Allocate(1);
    pts2[2].Allocate(1);
    cellTypeArray[2] =  2;
    matIdx[2].GetPortalControl().Set(0, 1);
    texIdx[2].GetPortalControl().Set(0, 1);
    pts1[2].GetPortalControl().Set(0, vec3(0,0,555));
    pts2[2].GetPortalControl().Set(0, vec3(555,555,555));
    flipped[2].GetPortalControl().Set(0, 1);
    list[i++] = new flip_normals(new xy_rect(0, 555, 0, 555, 555, 1,1));


    flipped[3].Allocate(1);
    matIdx[3].Allocate(1);
    texIdx[3].Allocate(1);
    pts1[3].Allocate(1);
    pts2[3].Allocate(1);
    cellTypeArray[3] =  3;
    matIdx[3].GetPortalControl().Set(0, 4);
    texIdx[3].GetPortalControl().Set(0, 0);
    pts1[3].GetPortalControl().Set(0, vec3(190,90,190));
    pts2[3].GetPortalControl().Set(0, vec3(90,0,0));
    flipped[3].GetPortalControl().Set(0, 0);
    list[i++] = new sphere(vec3(190, 90, 190),90 , 4, 0);

    flipped[4].Allocate(1);
    matIdx[4].Allocate(1);
    texIdx[4].Allocate(1);
    pts1[4].Allocate(1);
    pts2[4].Allocate(1);
    cellTypeArray[4] =  4;
    matIdx[4].GetPortalControl().Set(0, 1);
    texIdx[4].GetPortalControl().Set(0, 1);
    flipped[4].GetPortalControl().Set(0, 0);
    list[i++] = new translate(new rotate_y(
                    new box(vec3(0, 0, 0), vec3(165, 330, 165), 1,1),  15), vec3(265,0,295));
    *scene = new hitable_list(list,i);

    vec3 lookfrom(278, 278, -800);
    vec3 lookat(278,278,0);
    float dist_to_focus = 10.0;
    float aperture = 0.0;
    float vfov = 40.0;
    *cam = new camera(lookfrom, lookat, vec3(0,1,0),
                      vfov, aspect, aperture, dist_to_focus, 0.0, 1.0);
}


int main() {
  using MyAlgos = MyAlgorithms<vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>, VTKM_DEFAULT_DEVICE_ADAPTER_TAG>;
  using StorageTag = vtkm::cont::StorageTagBasic;

  constexpr int nx = 128;
  constexpr int ny = 128;
  constexpr int ns = 1;

  constexpr int depthcount = 2;
  auto canvasSize = nx*ny;

  //std::cout << "P3\n" << nx << " " << ny << "\n255\n";
  hitable *world;
  camera *cam;
  constexpr float aspect = float(ny) / float(nx);
  cornell_box(&world, &cam, aspect);
  hitable *light_shape = new xz_rect(213, 343, 227, 332, 554, 0,0);
  hitable *glass_sphere = new sphere(vec3(190, 90, 190), 90, 0,0);
  hitable *a[2];
  a[0] = light_shape;
  a[1] = glass_sphere;
  hitable_list hlist(a,2);

  vtkm::cont::ArrayHandle<ray> rays;
  rays.Allocate(nx*ny);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 2>> uvs;
  uvs.Allocate(nx*ny);

  vtkm::cont::ArrayHandleConstant<vec3> zero(vec3(0.0f), nx*ny);
  vtkm::cont::ArrayHandle<vec3> cols;
  cols.Allocate(nx*ny);
  MyAlgos::Copy<vec3, vec3, StorageTag>(vec3(0.0f), cols);

  vtkm::cont::ArrayHandle<vtkm::Float32> DirX, DirY, DirZ;
  DirX.Allocate(nx*ny); DirY.Allocate(nx*ny); DirZ.Allocate(nx*ny);
  vtkm::cont::ArrayHandle<vtkm::Id> PixelIdx;
  PixelIdx.Allocate(nx*ny);

  vtkm::cont::ArrayHandle<scatter_record> srecs;
  vtkm::cont::ArrayHandle<vtkm::Int8> scattered, hitArray, finished;

  vtkm::cont::ArrayHandle<hit_record> hrecs;

  vtkm::cont::ArrayHandle<float> closest, tmin;
  using ArrayType = vtkm::cont::ArrayHandle<vec3>;
  ArrayType attenuation;
  ArrayType emitted;
  attenuation.Allocate(rays.GetNumberOfValues() * depthcount);
  emitted.Allocate(rays.GetNumberOfValues() * depthcount);
  scattered.Allocate(rays.GetNumberOfValues());
  hitArray.Allocate(rays.GetNumberOfValues());
  finished.Allocate(rays.GetNumberOfValues());
  hrecs.Allocate(rays.GetNumberOfValues());
  srecs.Allocate(rays.GetNumberOfValues());
  ArrayType sumtotl;
  sumtotl.Allocate(rays.GetNumberOfValues());
  closest.Allocate(rays.GetNumberOfValues());
  tmin.Allocate(rays.GetNumberOfValues());
  for (int s =0; s<ns; s++){
    UVGen uvgen(nx, ny, s);
    vtkm::worklet::AutoDispatcherMapField<UVGen>(
          uvgen).Invoke(uvs);

    RayGen raygen(nx,ny, 40,40,
                  vtkm::Vec<vtkm::Float32,3>(0,0,1),
                  vtkm::Vec<vtkm::Float32,3>(0,1,0),
                  0, nx, 0, 0, ns);
    vtkm::worklet::AutoDispatcherMapField<RayGen>(raygen)
          .Invoke(DirX, DirY, DirZ, PixelIdx);

    vec3 lookfrom(278, 278, -800);
    RayLook rl(lookfrom);
    vtkm::worklet::AutoDispatcherMapField<RayLook>(rl)
          .Invoke(DirX, DirY, DirZ, uvs, rays);



    MyAlgos::Copy<vtkm::Int8, vtkm::Int8, StorageTag>(1, scattered);
    MyAlgos::Copy<vtkm::Int8, vtkm::Int8, StorageTag>(0, finished);


    for (int depth=0; depth<depthcount; depth++){
      MyAlgos::Copy<float, float, StorageTag>(std::numeric_limits<float>::max(), closest);
      MyAlgos::Copy<float, float, StorageTag>(0.001, tmin);
//      float tmin = 0.001;
//      float tmax = std::numeric_limits<float>::max();
      RayShade rs(world, canvasSize, depth);
//      vtkm::worklet::AutoDispatcherMapField<RayShade>(rs)
//            .Invoke(rays, hrecs, scattered, tex,  attenuation, emitted);

      MyAlgos::Copy<vtkm::Int8, vtkm::Int8, StorageTag>(0, hitArray);
      LambertianWorklet lmbWorklet( canvasSize, depth);
      DiffuseLightWorklet dlWorklet(canvasSize ,depth);
      DielectricWorklet deWorklet( canvasSize ,depth, 1.5, rays.GetNumberOfValues());
      PDFCosineWorklet pdfWorklet(canvasSize, depth, &hlist, rays.GetNumberOfValues());

      for (int i=0; i<cellTypeArray.size(); i++){
        if (cellTypeArray[i] == 0){
          YZRectWorklet yz(canvasSize, depth);
            vtkm::worklet::AutoDispatcherMapField<YZRectWorklet>(yz)
                .Invoke(rays, hrecs, tmin, closest, scattered, hitArray, pts1[i], pts2[i],
                        matIdx[i], texIdx[i],flipped[i]);
        }
        else if (cellTypeArray[i] == 1){
          XZRectWorklet xz(canvasSize, depth);
          vtkm::worklet::AutoDispatcherMapField<XZRectWorklet>(xz)
              .Invoke(rays, hrecs, tmin, closest, scattered, hitArray, pts1[i], pts2[i],
                      matIdx[i], texIdx[i],flipped[i]);


        }
        else if (cellTypeArray[i] == 2){
          //xy
          XYRectWorklet xy(canvasSize, depth);
          vtkm::worklet::AutoDispatcherMapField<XYRectWorklet>(xy)
              .Invoke(rays, hrecs, tmin, closest, scattered, hitArray, pts1[i], pts2[i],
                      matIdx[i], texIdx[i],flipped[i]);

        }
        else if (cellTypeArray[i] == 3){
          SphereIntersecttWorklet sphereIntersect(canvasSize, depth);
          vtkm::worklet::AutoDispatcherMapField<SphereIntersecttWorklet>(sphereIntersect)
              .Invoke(rays, hrecs, tmin, closest, scattered, hitArray, pts1[i], pts2[i],
                      matIdx[i], texIdx[i]);

        }
      }
#if 0
      for (int i=0; i<rays.GetNumberOfValues(); i++){

        auto r_in = rays.GetPortalConstControl().Get(i);
        auto hrec = hrecs.GetPortalConstControl().Get(i);
        auto sctr = scattered.GetPortalConstControl().Get(i);
        auto hit = hitArray.GetPortalControl().Get(i);

//        rs.operator()(i, r_in, hrec, sctr, tex.GetPortalControl(), attenuation.GetPortalControl(), emitted.GetPortalControl());
//        hrecs.GetPortalControl().Set(i, hrec);
//        scattered.GetPortalControl().Set(i,sctr);

        auto _tmin = tmin.GetPortalControl().Get(i);
        auto _tmax = closest.GetPortalControl().Get(i);

        auto pt1 = pts.GetPortalConstControl().Get(12);
        auto pt2 = pts.GetPortalConstControl().Get(13);
        auto mId = matIdx.GetPortalConstControl().Get(6);
        auto tId = texIdx.GetPortalConstControl().Get(6);
        sphereIntersect.operator ()(i, r_in, hrec,  _tmin, _tmax,sctr,hit, pt1, pt2, mId, tId,
                                    matType.GetPortalControl(),
                                    texType.GetPortalControl());

        pt1 = pts.GetPortalConstControl().Get(10);
        pt2 = pts.GetPortalConstControl().Get(11);
        mId = matIdx.GetPortalConstControl().Get(5);
        tId = texIdx.GetPortalConstControl().Get(5);
        XYRectWorklet xy(canvasSize, depth, true);
        xy.operator ()(i, r_in, hrec, _tmin, _tmax, sctr,hit, pt1, pt2, mId, tId,
            matType.GetPortalConstControl(),
            texType.GetPortalConstControl());

        pt1 = pts.GetPortalConstControl().Get(0);
        pt2 = pts.GetPortalConstControl().Get(1);
        mId = matIdx.GetPortalConstControl().Get(0);
        tId = texIdx.GetPortalConstControl().Get(0);
        YZRectWorklet yz(canvasSize, depth, true);
        yz.operator()(i, r_in, hrec, _tmin, _tmax, sctr, hit, pt1, pt2, mId, tId,
            matType.GetPortalConstControl(),
            texType.GetPortalConstControl());
        pt1 = pts.GetPortalConstControl().Get(2);
        pt2 = pts.GetPortalConstControl().Get(3);
        mId = matIdx.GetPortalConstControl().Get(1);
        tId = texIdx.GetPortalConstControl().Get(1);
        YZRectWorklet yz2(canvasSize, depth, false);

        yz2.operator()(i, r_in, hrec, _tmin, _tmax, sctr, hit, pt1, pt2, mId, tId,
            matType.GetPortalConstControl(),
            texType.GetPortalConstControl());

        pt1 = pts.GetPortalConstControl().Get(4);
        pt2 = pts.GetPortalConstControl().Get(5);
        mId = matIdx.GetPortalConstControl().Get(2);
        tId = texIdx.GetPortalConstControl().Get(2);
        XZRectWorklet xz(canvasSize, depth, true);
        xz.operator()(i, r_in, hrec, _tmin, _tmax, sctr, hit, pt1,pt2, mId, tId,
            matType.GetPortalConstControl(),
            texType.GetPortalConstControl());
        pt1 = pts.GetPortalConstControl().Get(6);
        pt2 = pts.GetPortalConstControl().Get(7);
        mId = matIdx.GetPortalConstControl().Get(3);
        tId = texIdx.GetPortalConstControl().Get(3);
        xz.operator()(i, r_in, hrec, _tmin, _tmax,sctr, hit, pt1,pt2, mId, tId,
            matType.GetPortalConstControl(),
            texType.GetPortalConstControl());
        pt1 = pts.GetPortalConstControl().Get(8);
        pt2 = pts.GetPortalConstControl().Get(9);
        mId = matIdx.GetPortalConstControl().Get(4);
        tId = texIdx.GetPortalConstControl().Get(4);
        XZRectWorklet xz2(canvasSize, depth, false);
        xz2.operator()(i, r_in, hrec, _tmin, _tmax,sctr, hit, pt1,pt2, mId, tId,
            matType.GetPortalConstControl(),
            texType.GetPortalConstControl());
#endif
      for (int i=0; i<rays.GetNumberOfValues(); i++){
        auto sctr = scattered.GetPortalControl().Get(i);
        auto hit = hitArray.GetPortalControl().Get(i);
        if (!(sctr && hit)){
          sctr = false;
          attenuation.GetPortalControl().Set(i + canvasSize * depth, vec3(1.0));
          emitted.GetPortalControl().Set(i + canvasSize * depth, vec3(0.0f));
        }

        scattered.GetPortalControl().Set(i,sctr);
      }
      vtkm::worklet::AutoDispatcherMapField<LambertianWorklet>(lmbWorklet)
          .Invoke(rays, hrecs, srecs, finished, scattered,
                  tex, matType, texType, emitted);

      vtkm::worklet::AutoDispatcherMapField<DiffuseLightWorklet>(dlWorklet)
          .Invoke(rays, hrecs, srecs, finished, scattered,
                  tex, matType, texType, emitted);

      vtkm::worklet::AutoDispatcherMapField<DielectricWorklet>(deWorklet)
            .Invoke(rays, hrecs, srecs, finished, scattered,
                    tex, matType, texType, emitted);

      vtkm::worklet::AutoDispatcherMapField<PDFCosineWorklet>(pdfWorklet)
            .Invoke(rays, hrecs, srecs, finished, scattered, rays, attenuation);
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

     std::cout << vtkm::cont::Algorithm::Reduce(emitted, vec3(0.0)) << std::endl;

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
  }

  std::vector<std::uint8_t> ImageBuffer;
  ImageBuffer.reserve(nx*ny*4);

  for (int i=0; i<cols.GetNumberOfValues(); i++){
    auto col = cols.GetPortalConstControl().Get(i);
    col = col / float(ns);
    col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
    int ir = int(255.99*col[0]);
    int ig = int(255.99*col[1]);
    int ib = int(255.99*col[2]);
    ImageBuffer.push_back(ir);
    ImageBuffer.push_back(ig);
    ImageBuffer.push_back(ib);
    ImageBuffer.push_back(255);
  }

//  std::vector<std::uint8_t> PngBuffer;
//  lodepng::encode(PngBuffer, ImageBuffer, nx, ny);
  lodepng::encode("output.png", ImageBuffer, nx,ny);
}

