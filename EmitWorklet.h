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

  using ControlSignature = void(FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, WholeArrayInOut<>, WholeArrayInOut<>, WholeArrayInOut<>);
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8);
  VTKM_EXEC
  template<typename VecArrayType,
          typename ColorArrayType,
          typename MatTypeArray>
  void operator()(vtkm::Id idx,
                  ray &r_in,
                  hit_record &hrec,
                  scatter_record &srec,
                  vtkm::Int8 &fin,
                  vtkm::Int8 &scattered,
                  ColorArrayType col,
                  MatTypeArray matType,
                  VecArrayType emitted) const
  {
    if (!fin){
      if (scattered){
        auto mt = matType.Get(hrec.matId);
        if (mt == 0){
          vec3 em = emit(r_in, hrec, col.Get(hrec.texId));
          scattered = scatter(r_in, hrec, srec, col.Get(hrec.texId));
          emitted.Set(idx * depthcount + depth, em);

        }
      }
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

  using ControlSignature = void(FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, WholeArrayInOut<>, WholeArrayInOut<>, WholeArrayInOut<>);
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8);
  VTKM_EXEC
  template<typename VecArrayType,
          typename ColorArrayType,
          typename MatTypeArray>
  void operator()(vtkm::Id idx,
                  ray &r_in,
                  hit_record &hrec,
                  scatter_record &srec,
                  vtkm::Int8 &fin,
                  vtkm::Int8 &scattered,
                  ColorArrayType col,
                  MatTypeArray matType,
                  VecArrayType emitted) const
  {
    if(!fin){
      if (scattered){
        auto mt = matType.Get(hrec.matId);
        if (mt == 1){
          vec3 em = emit(r_in, hrec, col.Get(hrec.texId));
          scattered = scatter(r_in, hrec, srec, col.Get(hrec.texId));
          emitted.Set(idx * depthcount + depth, em);
        }
      }
    }
  }

  const int depth, depthcount;
  int pdfIdx;

};
class DielectricWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  DielectricWorklet(int id, int dc, int d, float rid)
      :pdfIdx(id)
      , depthcount(dc)
      , depth(d)
      , ref_idx(rid)
  {
  }

  VTKM_EXEC
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, vec3 albedo) const {
      srec.is_specular = true;
      srec.pdfIdx = pdfIdx;
      srec.attenuation = vec3(1.0, 1.0, 1.0);
      vec3 outward_normal;
       vec3 reflected = reflect(r_in.direction(), hrec.normal);
       vec3 refracted;
       float ni_over_nt;
       float reflect_prob;
       float cosine;
       if (dot(r_in.direction(), hrec.normal) > 0) {
            outward_normal = -hrec.normal;
            ni_over_nt = ref_idx;
            cosine = ref_idx * dot(r_in.direction(), hrec.normal) * vtkm::RMagnitude(r_in.direction());
       }
       else {
            outward_normal = hrec.normal;
            ni_over_nt = 1.0 / ref_idx;
            cosine = -dot(r_in.direction(), hrec.normal) * vtkm::RMagnitude(r_in.direction());
       }
       if (refract(r_in.direction(), outward_normal, ni_over_nt, refracted)) {
          reflect_prob = schlick(cosine, ref_idx);
       }
       else {
          reflect_prob = 1.0;
       }
       if (drand48() < reflect_prob) {
          srec.specular_ray = ray(hrec.p, reflected);
       }
       else {
          srec.specular_ray = ray(hrec.p, refracted);
       }
       return true;
  }
  VTKM_EXEC
  vec3 emit(const ray& r_in, const hit_record& rec, vec3 emit) const { return vec3(0,0,0); }

  using ControlSignature = void(FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, WholeArrayInOut<>, WholeArrayInOut<>, WholeArrayInOut<>);
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8);
  VTKM_EXEC
  template<typename VecArrayType,
          typename ColorArrayType,
          typename MatTypeArray>
  void operator()(vtkm::Id idx,
                  ray &r_in,
                  hit_record &hrec,
                  scatter_record &srec,
                  vtkm::Int8 &fin,
                  vtkm::Int8 &scattered,
                  ColorArrayType col,
                  MatTypeArray matType,
                  VecArrayType emitted) const
  {
    if (!fin){
      if (scattered){
        auto mt = matType.Get(hrec.matId);
        if (mt == 2){
          vec3 em = emit(r_in, hrec, col.Get(hrec.texId));
          scattered = scatter(r_in, hrec, srec, col.Get(hrec.texId));
          emitted.Set(idx * depthcount + depth, em);

        }
      }
    }

  }

  const int depth, depthcount;
  int pdfIdx;
  float ref_idx;
};

#endif
