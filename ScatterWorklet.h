#ifndef SCATTERWORKLET_H
#define SCATTERWORKLET_H
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include "Record.h"
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
  float value(const vec3& direction) const {
      float cosine = dot(unit_vector(direction), uvw.w());
      if (cosine > 0)
          return cosine/M_PI;
      else
          return 0;
  }
  vec3 generate(float r1, float r2, float r3) const  {
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
  template<typename HitRecord>
  float scattering_pdf(const vec3 &origin, const vec3 &direction, const HitRecord& rec, const vec3& scattered_origin, const vec3 &scattered_direction) const {
    vec3 n(rec[static_cast<vtkm::Id>(HR::Nx)], rec[static_cast<vtkm::Id>(HR::Ny)], rec[static_cast<vtkm::Id>(HR::Nz)]);
      float cosine = dot(n, unit_vector(scattered_direction));
      if (cosine < 0)
          return 0;
      return cosine / M_PI;
  }

  using ControlSignature = void(FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, WholeArrayInOut<>);
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10);
  VTKM_EXEC
  template<typename VecArrayType, typename HitRecord,
  int FinBitIdx = 1,
  int ScatterBitIdx= 3>

  void operator()(vtkm::Id idx,
                  vec3 &r_origin,
                  vec3 &r_direction,
                  HitRecord &hrec,
                  ScatterRecord &srec,
                  vtkm::UInt8 &fin,
                  float &sum_value,
                  vec3 &generated_dir,
                  vec3 &out_origin,
                  vec3 &out_direction,
                  VecArrayType attenuation) const
  {
    if (! (fin & (1UL << FinBitIdx))){
      vec3 atten(1.0);
      out_origin = r_origin;
      out_direction = r_direction;

      if (fin & (1UL << ScatterBitIdx)){
        if (srec.is_specular) {
          atten = srec.attenuation;
          out_origin = srec.o;
          out_direction = srec.dir;
        }
        else {

          cosine_pdf newPdf;
          vec3 n(hrec[static_cast<vtkm::Id>(HR::Nx)], hrec[static_cast<vtkm::Id>(HR::Ny)], hrec[static_cast<vtkm::Id>(HR::Nz)]);
          newPdf.SetW(n);


          out_origin = vec3(hrec[static_cast<vtkm::Id>(HR::Px)], hrec[static_cast<vtkm::Id>(HR::Py)], hrec[static_cast<vtkm::Id>(HR::Pz)]);
          out_direction = generated_dir;
          //float pdf_val = p.value(r_out.direction());
          auto pdf_val = 0.5 * sum_value + 0.5 * newPdf.value(generated_dir);

          auto sctr = scattering_pdf(r_origin, r_direction, hrec, out_origin, out_direction)/pdf_val;
          atten = srec.attenuation * sctr;
        }
      }
      attenuation.Set(canvasSize * depth + idx, atten);
    }
    fin |= (fin >> ScatterBitIdx);
  }

  const int depth, canvasSize;
  vtkm::UInt32 RayCount;

  int type_size;
};
#endif
