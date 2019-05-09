#ifndef SCATTERWORKLET_H
#define SCATTERWORKLET_H
#include <vtkm/worklet/AutoDispatcherMapField.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/rendering/xorShift.h>
#include "ray.h"
#include "hitable.h"
#include "hitable_list.h"
#include "material.h"

class PDFCosineWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  PDFCosineWorklet( int dc, int d, hitable *ls)
      : depthcount(dc)
      , depth(d)
      , light_shape(ls)
  {
  }

  VTKM_EXEC
  float scattering_pdf(const ray& r_in, const hit_record& rec, const ray& scattered) const {
      float cosine = dot(rec.normal, unit_vector(scattered.direction()));
      if (cosine < 0)
          return 0;
      return cosine / M_PI;
  }

  using ControlSignature = void(FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, FieldInOut<>, WholeArrayInOut<>);
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7);
  VTKM_EXEC
  template<typename VecArrayType>
  void operator()(vtkm::Id idx,
                  ray &r_in,
                  hit_record &hrec,
                  vtkm::Int8 &is_scattered,
                  scatter_record &srec,
                  ray &r_out,
                  VecArrayType attenuation) const
  {
    vec3 atten;
    if (is_scattered){
      if (srec.is_specular) {
        atten = srec.attenuation;
        r_out = srec.specular_ray;
        //return std::make_tuple(3, srec.attenuation, vec3(0,0,0), srec.specular_ray);
      }
      else {
        cosine_pdf newPdf;
        newPdf.SetW(hrec.normal);
        hitable_pdf plight(light_shape, hrec.p);
        mixture_pdf p(&plight, &newPdf);
        r_out = ray(hrec.p, p.generate(), r_in.time());
        float pdf_val = p.value(r_out.direction());
        atten = srec.attenuation * scattering_pdf(r_in, hrec, r_out)/pdf_val;
        //auto ret = scatter(r, srec, hrec, light_shape, mat);
        //return std::make_tuple(2, srec.attenuation*std::get<0>(ret), emitted, std::get<1>(ret));
      }
    }
    else{
      atten = vec3(1.0);
      r_out = r_in;
    }
    attenuation.Set(idx * depthcount + depth, atten);
  }

  const int depth, depthcount;
  hitable *light_shape;
};
#endif
