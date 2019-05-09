#ifndef EMITWORKLET_H
#define EMITWORKLET_H
#include <vtkm/worklet/AutoDispatcherMapField.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/rendering/xorShift.h>
#include "ray.h"
#include "hitable.h"
#include "hitable_list.h"
#include "material.h"

class LambertianWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  LambertianWorklet(int id, int dc, int d)
      :pdfIdx(id)
      , depthcount(dc)
      , depth(d)
  {
  }

  VTKM_EXEC
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, vec3 albedo) const {
      srec.is_specular = false;
      srec.attenuation = albedo;
      srec.pdfIdx = pdfIdx;//std::make_shared<cosine_pdf>(hrec.normal);
      return true;
  }
  VTKM_EXEC
  vec3 emit(const ray& r_in, const hit_record& rec, vec3 emit) const { return vec3(0,0,0); }

  using ControlSignature = void(FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, WholeArrayInOut<>, WholeArrayInOut<>, WholeArrayInOut<>, WholeArrayInOut<>);
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8);
  VTKM_EXEC
  template<typename VecArrayType,
          typename ColorArrayType,
          typename MatTypeArray>
  void operator()(vtkm::Id idx,
                  ray &r_in,
                  hit_record &hrec,
                  scatter_record &srec,
                  vtkm::Int8 &scattered,
                  ColorArrayType col,
                  MatTypeArray matType,
                  VecArrayType emitted,
                  VecArrayType attenuation) const
  {

    if (scattered){
      auto mt = matType.Get(hrec.matId);
      if (mt == 0){
        vec3 em = emit(r_in, hrec, col.Get(hrec.texId));
        scattered = scatter(r_in, hrec, srec, col.Get(hrec.texId));
        emitted.Set(idx * depthcount + depth, em);

      }
    }
    else{
      emitted.Set(idx * depthcount + depth, vec3(0.0f));
      attenuation.Set(idx * depthcount + depth, vec3(1.0f));
    }

  }

  const int depth, depthcount;
  int pdfIdx;
};

class DiffuseLightWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  DiffuseLightWorklet(int id, int dc, int d)
    :pdfIdx(id)
    , depthcount(dc)
    , depth(d)
  {
  }

  VTKM_EXEC
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, vec3) const {
        return false;}
  VTKM_EXEC
  vec3 emit(const ray& r_in, const hit_record& rec, vec3 emit) const {
      if (dot(rec.normal, r_in.direction()) < 0.0)
          return emit;
      else
          return vec3(0,0,0);
  }

  using ControlSignature = void(FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, WholeArrayInOut<>, WholeArrayInOut<>, WholeArrayInOut<>, WholeArrayInOut<>);
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8);
  VTKM_EXEC
  template<typename VecArrayType,
          typename ColorArrayType,
          typename MatTypeArray>
  void operator()(vtkm::Id idx,
                  ray &r_in,
                  hit_record &hrec,
                  scatter_record &srec,
                  vtkm::Int8 &scattered,
                  ColorArrayType col,
                  MatTypeArray matType,
                  VecArrayType emitted,
                  VecArrayType attenuation) const
  {
    if (scattered){
      auto mt = matType.Get(hrec.matId);
      if (mt == 1){
        vec3 em = emit(r_in, hrec, col.Get(hrec.texId));
        scattered = scatter(r_in, hrec, srec, col.Get(hrec.texId));
        emitted.Set(idx * depthcount + depth, em);
      }
    }
    else{
      emitted.Set(idx * depthcount + depth, vec3(0.0f));
      attenuation.Set(idx * depthcount + depth, vec3(1.0f));
    }
  }

  const int depth, depthcount;
  int pdfIdx;

};

#endif
