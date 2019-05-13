#ifndef SCATTERWORKLET_H
#define SCATTERWORKLET_H
#include <vtkm/worklet/AutoDispatcherMapField.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/rendering/xorShift.h>
#include "onb.h"

class cosine_pdf {
    public:
  inline vec3 random_cosine_direction(float r1, float r2) const {
      float z = sqrt(1-r2);
      float phi = 2*M_PI*r1;
      float x = cos(phi)*2*sqrt(r2);
      float y = sin(phi)*2*sqrt(r2);
      return vec3(x, y, z);
  }
  void SetW(const vec3& w) { uvw.build_from_w(w); }
  virtual float value(const vec3& direction) const {
      float cosine = dot(unit_vector(direction), uvw.w());
      if (cosine > 0)
          return cosine/M_PI;
      else
          return 0;
  }
  virtual vec3 generate(float r1, float r2, float r3) const  {
      return uvw.local(random_cosine_direction(r1,r2));
  }
  onb uvw;
};

class PDFCosineWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  PDFCosineWorklet( int cs, int d, vtkm::UInt32 rc, int ts)
      : canvasSize(cs)
      , depth(d)
      , RayCount(rc)
      , type_size(ts)
  {
  }



  VTKM_EXEC
  float scattering_pdf(const ray& r_in, const HitRecord& rec, const ray& scattered) const {
      float cosine = dot(rec.normal, unit_vector(scattered.direction()));
      if (cosine < 0)
          return 0;
      return cosine / M_PI;
  }

  using ControlSignature = void(FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, WholeArrayInOut<>);
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9);
  VTKM_EXEC
  template<typename VecArrayType>
  void operator()(vtkm::Id idx,
                  ray &r_in,
                  HitRecord &hrec,
                  ScatterRecord &srec,
                  vtkm::Int8 &fin,
                  vtkm::Int8 &is_scattered,
                  float &sum_value,
                  vec3 &direction,
                  ray &r_out,
                  VecArrayType attenuation) const
  {
    if (!fin){
      vec3 atten(1.0);
      r_out = r_in;

      if (is_scattered){
        if (srec.is_specular) {
          atten = srec.attenuation;
          r_out = srec.specular_ray;
        }
        else {

          cosine_pdf newPdf;
          newPdf.SetW(hrec.normal);

          r_out = ray(hrec.p, direction, r_in.time());

          //float pdf_val = p.value(r_out.direction());
          auto pdf_val = 0.5 * sum_value + 0.5 * newPdf.value(direction);

          auto sctr = scattering_pdf(r_in, hrec, r_out)/pdf_val;
          atten = srec.attenuation * sctr;
        }
      }
      attenuation.Set(canvasSize * depth + idx, atten);
    }
    fin = !is_scattered;
  }

  const int depth, canvasSize;
  vtkm::UInt32 RayCount;

  int type_size;
};
#endif
