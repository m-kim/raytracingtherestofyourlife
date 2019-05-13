#ifndef EMITWORKLET_H
#define EMITWORKLET_H
#include <vtkm/worklet/AutoDispatcherMapField.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/rendering/xorShift.h>

class LambertianWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  LambertianWorklet(int cs, int d)
      : canvasSize(cs)
      , depth(d)
  {
  }

  VTKM_EXEC
  bool scatter(const ray& r_in, const HitRecord& hrec, ScatterRecord& srec, vec3 albedo) const {
      srec.is_specular = false;
      srec.attenuation = albedo;
      return true;
  }
  VTKM_EXEC
  vec3 emit(const ray& r_in, const HitRecord& rec, vec3 emit) const { return vec3(0,0,0); }

  using ControlSignature = void(FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, WholeArrayInOut<>,WholeArrayInOut<>, WholeArrayInOut<>, WholeArrayInOut<>);
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8,_9);
  VTKM_EXEC
  template<typename VecArrayType,
          typename ColorArrayType,
          typename MatTypeArray,
          typename TexTypeArray>
  void operator()(vtkm::Id idx,
                  ray &r_in,
                  HitRecord &hrec,
                  ScatterRecord &srec,
                  vtkm::Int8 &fin,
                  vtkm::Int8 &scattered,
                  ColorArrayType col,
                  MatTypeArray matType,
                  TexTypeArray texType,
                  VecArrayType emitted) const
  {
    if (!fin){
      if (scattered){
        auto mt = matType.Get(hrec.matId);
        if (mt == 0){
          vec3 em = emit(r_in, hrec, col.Get(texType.Get(hrec.texId)));
          scattered = scatter(r_in, hrec, srec, col.Get(texType.Get(hrec.texId)));
          emitted.Set(canvasSize * depth + idx, em);

        }
      }
    }

  }

  const int depth, canvasSize;
};

class DiffuseLightWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  DiffuseLightWorklet(int cs, int d)
    : canvasSize(cs)
    , depth(d)
  {
  }

  VTKM_EXEC
  bool scatter(const ray& r_in, const HitRecord& hrec, ScatterRecord& srec, vec3) const {
        return false;}
  VTKM_EXEC
  vec3 emit(const ray& r_in, const HitRecord& rec, vec3 emit) const {
      if (dot(rec.normal, r_in.direction()) < 0.0)
          return emit;
      else
          return vec3(0,0,0);
  }

  using ControlSignature = void(FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, WholeArrayInOut<>, WholeArrayInOut<>, WholeArrayInOut<>, WholeArrayInOut<>);
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9);
  VTKM_EXEC
  template<typename VecArrayType,
          typename ColorArrayType,
          typename MatTypeArray,
          typename TexTypeArray>
  void operator()(vtkm::Id idx,
                  ray &r_in,
                  HitRecord &hrec,
                  ScatterRecord &srec,
                  vtkm::Int8 &fin,
                  vtkm::Int8 &scattered,
                  ColorArrayType col,
                  MatTypeArray matType,
                  TexTypeArray texType,
                  VecArrayType emitted) const
  {
    if(!fin){
      if (scattered){
        auto mt = matType.Get(hrec.matId);
        if (mt == 1){
          vec3 em = emit(r_in, hrec, col.Get(texType.Get(hrec.texId)));
          scattered = scatter(r_in, hrec, srec, col.Get(texType.Get(hrec.texId)));
          emitted.Set(canvasSize * depth + idx, em);
        }
      }
    }
  }

  const int depth, canvasSize;

};
class DielectricWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  DielectricWorklet(int cs, int d, float rid, vtkm::UInt32 rc)
      : canvasSize(cs)
      , depth(d)
      , ref_idx(rid)
      , RayCount(rc)
  {
  }

  VTKM_EXEC
  float schlick(float cosine, float ref_idx) const {
      float r0 = (1-ref_idx) / (1+ref_idx);
      r0 = r0*r0;
      return r0 + (1-r0)*pow((1 - cosine),5);
  }

  VTKM_EXEC
  bool refract(const vec3& v, const vec3& n, float ni_over_nt, vec3& refracted) const {
      vec3 uv = unit_vector(v);
      float dt = dot(uv, n);
      float discriminant = 1.0 - ni_over_nt*ni_over_nt*(1-dt*dt);
      if (discriminant > 0) {
          refracted = ni_over_nt*(uv - n*dt) - n*sqrt(discriminant);
          return true;
      }
      else
          return false;
  }

  VTKM_EXEC
  vec3 reflect(const vec3& v, const vec3& n) const {
       return v - 2*dot(v,n)*n;
  }

  VTKM_EXEC
  bool scatter(const ray& r_in, const HitRecord& hrec, ScatterRecord& srec, vec3 albedo, double _rand) const {
      srec.is_specular = true;
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
       if (_rand < reflect_prob) {
          srec.specular_ray = ray(hrec.p, reflected);
       }
       else {
          srec.specular_ray = ray(hrec.p, refracted);
       }
       return true;
  }
  VTKM_EXEC
  vec3 emit(const ray& r_in, const HitRecord& rec, vec3 emit) const { return vec3(0,0,0); }

  using ControlSignature = void(FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, WholeArrayInOut<>, WholeArrayInOut<>, WholeArrayInOut<>, WholeArrayInOut<>);
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9);
  VTKM_EXEC
  template<typename VecArrayType,
          typename ColorArrayType,
          typename MatTypeArray,
          typename TexTypeArray>
  void operator()(vtkm::Id idx,
                  ray &r_in,
                  HitRecord &hrec,
                  ScatterRecord &srec,
                  vtkm::Int8 &fin,
                  vtkm::Int8 &scattered,
                  ColorArrayType col,
                  MatTypeArray matType,
                  TexTypeArray texType,
                  VecArrayType emitted) const
  {
    if (!fin){
      if (scattered){
        auto mt = matType.Get(hrec.matId);
        if (mt == 2){
          vtkm::Vec<vtkm::UInt32, 4> randState;
          randState[0] = vtkm::random::xorshift::getRand32(idx*1) + 1;
          randState[1] = vtkm::random::xorshift::getRand32(idx*2) + 2;
          randState[2] = vtkm::random::xorshift::getRand32(idx*3) + 3;
          randState[3] = vtkm::random::xorshift::getRand32(idx*4) + 4; //arbitrary random state based off number of rays being shot through

          vec3 em = emit(r_in, hrec, col.Get(texType.Get(hrec.texId)));
          scattered = scatter(r_in, hrec, srec, col.Get(texType.Get(hrec.texId)), drand48());//vtkm::random::xorshift::getRandF(randState));
          emitted.Set(canvasSize * depth + idx, em);

        }
      }
    }

  }

  const int depth, canvasSize;
  float ref_idx;
  const vtkm::UInt32 RayCount;

};

#endif
