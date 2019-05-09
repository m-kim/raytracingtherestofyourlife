#ifndef WORKLETS_H
#define WORKLETS_H

#include <vtkm/worklet/AutoDispatcherMapField.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/rendering/xorShift.h>
#include "ray.h"
#include "hitable.h"
#include "hitable_list.h"
#include "material.h"

using vec3 = vtkm::Vec<vtkm::Float32, 3>;

class UVGen : public vtkm::worklet::WorkletMapField
{
public:
  vtkm::UInt32 RayCount;
  VTKM_CONT

  UVGen(vtkm::Id _nx,
        vtkm::Id _ny,
        vtkm::UInt32 rayCount)
    : numx(_nx),
      numy(_ny),
      RayCount(rayCount)
  {
  }

  using ControlSignature = void(FieldOut<>);

  using ExecutionSignature = void(WorkIndex, _1);
  template <typename Precision>
  VTKM_EXEC void operator()(vtkm::Id &idx,
                            Precision& uv) const
  {
    vtkm::Vec<vtkm::UInt32, 4> randState;
    randState[0] = vtkm::random::xorshift::getRand32(RayCount*1) + 1;
    randState[1] = vtkm::random::xorshift::getRand32(RayCount*2) + 2;
    randState[2] = vtkm::random::xorshift::getRand32(RayCount*3) + 3;
    randState[3] = vtkm::random::xorshift::getRand32(RayCount*4) + 4; //arbitrary random state based off number of rays being shot through

    vtkm::random::xorshift::getRandF(randState);
    vtkm::random::xorshift::getRandF(randState);
    vtkm::random::xorshift::getRandF(randState);
    vtkm::random::xorshift::getRandF(randState);

    vtkm::Id i = idx % numx;
    vtkm::Id j = idx / numx;
    uv[0] = vtkm::Float32(i + vtkm::random::xorshift::getRandF(randState)) / numx;
    uv[1] = vtkm::Float32(j + vtkm::random::xorshift::getRandF(randState)) / numy;

  }

  const vtkm::Id numx, numy;
};


class RayGen : public vtkm::worklet::WorkletMapField
{

public:
  vtkm::Int32 numx;
  vtkm::Int32 numy;
  vtkm::Int32 Minx;
  vtkm::Int32 Miny;
  vtkm::Int32 SubsetWidth;
  vtkm::Vec<vtkm::Float32, 3> nlook; // normalized look
  vtkm::Vec<vtkm::Float32, 3> delta_x;
  vtkm::Vec<vtkm::Float32, 3> delta_y;
  vtkm::UInt32 RayCount;
  VTKM_CONT
  RayGen(vtkm::Int32 _nx,
                    vtkm::Int32 _ny,
                    vtkm::Float32 fovX,
                    vtkm::Float32 fovY,
                    vtkm::Vec<vtkm::Float32, 3> look,
                    vtkm::Vec<vtkm::Float32, 3> up,
                    vtkm::Float32 _zoom,
                    vtkm::Int32 subsetWidth,
                    vtkm::Int32 minx,
                    vtkm::Int32 miny,
                    vtkm::UInt32 rayCount)
    : numx(_nx)
    , numy(_ny)
    , Minx(minx)
    , Miny(miny)
    , SubsetWidth(subsetWidth)
    , RayCount(rayCount)
  {
    vtkm::Float32 thx = tanf((fovX * vtkm::Pi_180f()) * .5f);
    vtkm::Float32 thy = tanf((fovY * vtkm::Pi_180f()) * .5f);
    vtkm::Vec<vtkm::Float32, 3> u = vtkm::Cross(look, up);
    vtkm::Normalize(u);

    vtkm::Vec<vtkm::Float32, 3> v = vtkm::Cross(u, look);
    vtkm::Normalize(v);
    delta_x = u * (2 * thx / (float)numx);
    delta_y = v * (2 * thy / (float)numy);

    if (_zoom > 0)
    {
      delta_x[0] = delta_x[0] / _zoom;
      delta_x[1] = delta_x[1] / _zoom;
      delta_x[2] = delta_x[2] / _zoom;
      delta_y[0] = delta_y[0] / _zoom;
      delta_y[1] = delta_y[1] / _zoom;
      delta_y[2] = delta_y[2] / _zoom;
    }
    nlook = look;
    vtkm::Normalize(nlook);
  }

  using ControlSignature = void(FieldOut<>, FieldOut<>, FieldOut<>, FieldOut<>);

  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4);

  template <typename Precision>
  VTKM_EXEC void operator()(vtkm::Id idx,
                            Precision& rayDirX,
                            Precision& rayDirY,
                            Precision& rayDirZ,
                            vtkm::Id& pixelIndex) const
  {
    vtkm::Vec<vtkm::UInt32, 4> randState;
    randState[0] = vtkm::random::xorshift::getRand32(RayCount*1) + 1;
    randState[1] = vtkm::random::xorshift::getRand32(RayCount*2) + 2;
    randState[2] = vtkm::random::xorshift::getRand32(RayCount*3) + 3;
    randState[3] = vtkm::random::xorshift::getRand32(RayCount*4) + 4;

    vtkm::random::xorshift::getRandF(randState);
    vtkm::random::xorshift::getRandF(randState);
    vtkm::random::xorshift::getRandF(randState);
    vtkm::random::xorshift::getRandF(randState);

    vtkm::Vec<Precision, 3> ray_dir(rayDirX, rayDirY, rayDirZ);
    int i = vtkm::Int32(idx) % SubsetWidth;
    int j = vtkm::Int32(idx) / SubsetWidth;
    i += Minx;
    j += Miny;
    // Write out the global pixelId
    pixelIndex = static_cast<vtkm::Id>(j * numx + i);
    // ray_dir = nlook + delta_x * ((2.f * Precision(i) - Precision(w)) / 2.0f) +
      // delta_y * ((2.f * Precision(j) - Precision(h)) / 2.0f);

    Precision _randU = vtkm::random::xorshift::getRandF(randState);
    Precision _randV = vtkm::random::xorshift::getRandF(randState);

    if (RayCount < 2) {_randU = 0.5f; _randV = 0.5f;}

    ray_dir = nlook + delta_x * ((2.f * (Precision(i) + (1.f - _randU)) - Precision(numx)) / 2.0f) +
      delta_y * ((2.f * (Precision(j) + (_randV)) - Precision(numy)) / 2.0f);

      // if (idx == 15427)
      // {
        // DBGVAR(ray_dir);
      // }
    // avoid some numerical issues
    for (vtkm::Int32 d = 0; d < 3; ++d)
    {
      if (ray_dir[d] == 0.f)
        ray_dir[d] += 0.0000001f;
    }
    Precision dot = vtkm::Dot(ray_dir, ray_dir);
    Precision sq_mag = vtkm::Sqrt(dot);

    rayDirX = ray_dir[0] / sq_mag;
    rayDirY = ray_dir[1] / sq_mag;
    rayDirZ = ray_dir[2] / sq_mag;
  }

}; // class perspective ray gen

class RayLook : public vtkm::worklet::WorkletMapField
{
public:
  vtkm::UInt32 RayCount;
  VTKM_CONT

  RayLook(vec3 lf)
    : lookfrom(lf)
  {
  }

  using ControlSignature = void(FieldIn<>, FieldIn<>, FieldIn<>, FieldIn<>, FieldOut<>);

  using ExecutionSignature = void(_1, _2, _3, _4, _5);
  template <typename Precision, typename UVType>
  VTKM_EXEC void operator()(
                            Precision x,
                            Precision y,
                            Precision z,
                            UVType &uv,
                            ray &r) const
  {
    //auto uv = uvs.GetPortalConstControl().Get(i);
    vec3 dir(x,y,z);
    r = ray(lookfrom, dir);
  }
  vec3 lookfrom;

};

class RayShade : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  RayShade(hitable *w, hitable_list hl, int dc, int d)
    :world(w)
    ,hlist(hl)
    ,depthcount(dc)
    ,depth(d)
  {
  }


  VTKM_EXEC
  template<typename MatType, typename TextureType>
  auto color(const ray& r, const hit_record & hrec, const hitable *light_shape, MatType &mat, TextureType &tex) const {

    scatter_record srec;
    vec3 emitted = mat.emitted(r, hrec, tex.value(hrec.u, hrec.v, hrec.p));
    if (mat.scatter(r, hrec, srec, tex.value(hrec.u, hrec.v, hrec.p))) {
        if (srec.is_specular) {
          return std::make_tuple(3, srec.attenuation, vec3(0,0,0), srec.specular_ray);
        }
        else {
          pdf_ptrs[srec.pdfIdx]->SetW(hrec.normal);
            hitable_pdf plight(light_shape, hrec.p);
            mixture_pdf p(&plight,pdf_ptrs[srec.pdfIdx] );
            ray scattered = ray(hrec.p, p.generate(), r.time());
            float pdf_val = p.value(scattered.direction());
            return std::make_tuple(2, srec.attenuation*mat.scattering_pdf(r, hrec, scattered)/pdf_val,
                                   emitted, scattered);
        }
    }
    else
      return std::make_tuple(1, vec3(1.0f), emitted, r);

  }

  using ControlSignature = void(FieldInOut<>, FieldInOut<>, FieldInOut<>,
  WholeArrayInOut<>, WholeArrayInOut<>);

  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5);
  VTKM_EXEC
  template<typename VecArrayType>
  void operator()(vtkm::Id idx,
            ray &ray_io, hit_record &hrec, vtkm::UInt8 &fin,
            VecArrayType attenuation,
             VecArrayType emitted) const
  {
    vtkm::Int8 state;
    vec3 att, em;

    if (!fin && world->hit(ray_io, 0.001, std::numeric_limits<float>::max(), hrec)){
      std::tie(state, att, em, ray_io) = color(ray_io, hrec, &hlist, *mat_ptrs[hrec.matId], *tex_ptrs[hrec.texId]);
      attenuation.Set(idx * depthcount + depth, att);
      emitted.Set(idx * depthcount + depth, em);
      if (state < 2){
        fin = 1;
      }
    }
    else{
      attenuation.Set(idx * depthcount + depth, vec3(1.0));
      emitted.Set(idx * depthcount + depth, vec3(0.0f));
    }
  }

  std::vector<texture*> tex_ptrs;
  std::vector<material*> mat_ptrs;
  std::vector<pdf*> pdf_ptrs;

  hitable *world;
  hitable_list hlist;
  int depthcount;
  int depth;

};
#endif
