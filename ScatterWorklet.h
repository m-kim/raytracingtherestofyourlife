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
  PDFCosineWorklet( int cs, int d, hitable *ls, vtkm::UInt32 rc)
      : canvasSize(cs)
      , depth(d)
      , light_shape(ls)
      , RayCount(rc)
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
                  scatter_record &srec,
                  vtkm::Int8 &fin,
                  vtkm::Int8 &is_scattered,
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
          //return std::make_tuple(3, srec.attenuation, vec3(0,0,0), srec.specular_ray);
        }
        else {
//          vtkm::Vec<vtkm::UInt32, 4> randState;
//          randState[0] = vtkm::random::xorshift::getRand32(idx*1) + 1;
//          randState[1] = vtkm::random::xorshift::getRand32(idx*2) + 2;
//          randState[2] = vtkm::random::xorshift::getRand32(idx*3) + 3;
//          randState[3] = vtkm::random::xorshift::getRand32(idx*4) + 4; //arbitrary random state based off number of rays being shot through
          float r1 = drand48();
          float r2 = drand48();
          float r3 = drand48();
#if 0
          onb uvw;
          vec3 direction;

          //generate
          uvw.build_from_w(hrec.normal);
          if (r3 < 0.5)
            direction = uvw.local(random_cosine_direction(r1,r2));
          else
            direction = light_shape->random(hrec.p, r1,r2,r3);

          float pdf_val = 0.5 * value(direction, uvw);
          pdf_val += 0.5 * light_shape->pdf_value(hrec.p, direction);

          r_out = ray(hrec.p, direction, r_in.time());
          atten = srec.attenuation * scattering_pdf(r_in, hrec, r_out)/pdf_val;
#else
          cosine_pdf newPdf;

          newPdf.SetW(hrec.normal);
          hitable_pdf plight(light_shape, hrec.p);
          mixture_pdf p(&plight, &newPdf);
          r_out = ray(hrec.p, p.generate(r1,r2,r3), r_in.time());
          float pdf_val = p.value(r_out.direction());
          auto sctr = scattering_pdf(r_in, hrec, r_out)/pdf_val;
          atten = srec.attenuation * sctr;
          //auto ret = scatter(r, srec, hrec, light_shape, mat);
          //return std::make_tuple(2, srec.attenuation*std::get<0>(ret), emitted, std::get<1>(ret));
#endif
        }
      }
      attenuation.Set(canvasSize * depth + idx, atten);
    }
    fin = !is_scattered;
  }

  const int depth, canvasSize;
  hitable *light_shape;
  vtkm::UInt32 RayCount;
};
#endif
