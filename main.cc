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
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/rendering/raytracing/Camera.h>

#include "Worklets.h"
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


inline vec3 de_nan(const vec3& c) {
    vec3 temp = c;
    if (!(temp[0] == temp[0])) temp[0] = 0;
    if (!(temp[1] == temp[1])) temp[1] = 0;
    if (!(temp[2] == temp[2])) temp[2] = 0;
    return temp;
}



auto color(const ray& r, const hit_record & hrec, hitable *light_shape) {

  scatter_record srec;
  vec3 emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v, hrec.p);
  if (hrec.mat_ptr->scatter(r, hrec, srec)) {
      if (srec.is_specular) {
        return std::make_tuple(3, srec.attenuation, vec3(0,0,0), srec.specular_ray);
      }
      else {
          hitable_pdf plight(light_shape, hrec.p);
          mixture_pdf p(&plight, srec.pdf_ptr.get());
          ray scattered = ray(hrec.p, p.generate(), r.time());
          float pdf_val = p.value(scattered.direction());
          return std::make_tuple(2, srec.attenuation*hrec.mat_ptr->scattering_pdf(r, hrec, scattered)/pdf_val,
                                 emitted, scattered);
      }
  }
  else
    return std::make_tuple(1, vec3(1.0f), emitted, r);

}

void cornell_box(hitable **scene, camera **cam, float aspect) {
    int i = 0;
    hitable **list = new hitable*[8];
    material *red = new lambertian( new constant_texture(vec3(0.65, 0.05, 0.05)) );
    material *white = new lambertian( new constant_texture(vec3(0.73, 0.73, 0.73)) );
    material *green = new lambertian( new constant_texture(vec3(0.12, 0.45, 0.15)) );
    material *light = new diffuse_light( new constant_texture(vec3(15, 15, 15)) );
    list[i++] = new flip_normals(new yz_rect(0, 555, 0, 555, 555, green));
    list[i++] = new yz_rect(0, 555, 0, 555, 0, red);
    list[i++] = new flip_normals(new xz_rect(213, 343, 227, 332, 554, light));
    list[i++] = new flip_normals(new xz_rect(0, 555, 0, 555, 555, white));
    list[i++] = new xz_rect(0, 555, 0, 555, 0, white);
    list[i++] = new flip_normals(new xy_rect(0, 555, 0, 555, 555, white));
    material *glass = new dielectric(1.5);
    list[i++] = new sphere(vec3(190, 90, 190),90 , glass);
    list[i++] = new translate(new rotate_y(
                    new box(vec3(0, 0, 0), vec3(165, 330, 165), white),  15), vec3(265,0,295));
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
  constexpr int nx = 128;
  constexpr int ny = 128;
  constexpr int ns = 10;

  constexpr int depthcount = 50;
  //std::cout << "P3\n" << nx << " " << ny << "\n255\n";
  hitable *world;
  camera *cam;
  constexpr float aspect = float(ny) / float(nx);
  cornell_box(&world, &cam, aspect);
  hitable *light_shape = new xz_rect(213, 343, 227, 332, 554, 0);
  hitable *glass_sphere = new sphere(vec3(190, 90, 190), 90, 0);
  hitable *a[2];
  a[0] = light_shape;
  a[1] = glass_sphere;
  hitable_list hlist(a,2);

  vtkm::cont::ArrayHandle<ray> rays;
  rays.Allocate(nx*ny);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 2>> uvs;
  uvs.Allocate(nx*ny);

  vtkm::cont::ArrayHandle<vec3> cols;
  cols.Allocate(nx*ny);
  vtkm::cont::ArrayHandleConstant<vec3> zero(vec3(0,0,0), nx*ny);
  vtkm::cont::ArrayCopy(zero, cols);

  vtkm::cont::ArrayHandle<vtkm::Float32> DirX, DirY, DirZ;
  DirX.Allocate(nx*ny); DirY.Allocate(nx*ny); DirZ.Allocate(nx*ny);
  vtkm::cont::ArrayHandle<vtkm::Id> PixelIdx;
  PixelIdx.Allocate(nx*ny);
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

    for (int i=0; i<rays.GetNumberOfValues(); i++){
      auto uv = uvs.GetPortalConstControl().Get(i);
      //rays.GetPortalControl().Set(i, cam->get_ray(uv));
      vec3 lookfrom(278, 278, -800);
      auto x= DirX.GetPortalControl().Get(i);
      auto y= DirY.GetPortalControl().Get(i);
      auto z= DirZ.GetPortalControl().Get(i);
      vec3 dir(x,y,z);
      rays.GetPortalControl().Set(i, ray(lookfrom, dir));
    }

    using ArrayType = vtkm::cont::ArrayHandle<vec3>;
    ArrayType attenuation;
    ArrayType emitted;
    attenuation.Allocate(rays.GetNumberOfValues() * depthcount);
    emitted.Allocate(rays.GetNumberOfValues() * depthcount);

    vtkm::cont::ArrayHandle<vtkm::UInt8> finished;
    finished.Allocate(rays.GetNumberOfValues());
    vtkm::cont::ArrayHandle<hit_record> hrecs;
    hrecs.Allocate(rays.GetNumberOfValues());
    for (int i=0; i<rays.GetNumberOfValues(); i++)
    {
      finished.GetPortalControl().Set(i, 0);
    }
    for (int depth=0; depth<depthcount; depth++){
      for (int i=0; i<rays.GetNumberOfValues(); i++)
      {
          //vec3 p = r.point_at_parameter(2.0);
        auto ray = rays.GetPortalConstControl().Get(i);
        vtkm::Int8 state;

        vec3 att, em;
        auto hrec = hrecs.GetPortalControl().Get(i);
        auto fin = finished.GetPortalConstControl().Get(i);
        if (!fin && world->hit(ray, 0.001, std::numeric_limits<float>::max(), hrec)){
          std::tie(state, att, em, ray) = color(ray, hrec, &hlist);
          attenuation.GetPortalControl().Set(i * depthcount + depth, att);
          emitted.GetPortalControl().Set(i * depthcount + depth, em);
          if (state < 2){
            finished.GetPortalControl().Set(i, 1 );
          }
        }
        else{
          attenuation.GetPortalControl().Set(i * depthcount + depth, vec3(1.0));
          emitted.GetPortalControl().Set(i * depthcount + depth, vec3(0.0f));
        }

        rays.GetPortalControl().Set(i,ray);
      }
    }

    ArrayType sumtotl;
    sumtotl.Allocate(rays.GetNumberOfValues());
    for (int i=0; i<sumtotl.GetNumberOfValues(); i++){
      sumtotl.GetPortalControl().Set(i, emitted.GetPortalConstControl().Get(i * depthcount + depthcount -1));
    }
    for (int depth = depthcount-2; depth >=0; depth--){
      for (int i=0; i<sumtotl.GetNumberOfValues(); i++){
        auto sum = sumtotl.GetPortalConstControl().Get(i);
        sumtotl.GetPortalControl().Set(i, emitted.GetPortalConstControl().Get(i*depthcount + depth) + attenuation.GetPortalConstControl().Get(i*depthcount + depth) * sum);
      }
    }
    for (int i=0; i<cols.GetNumberOfValues(); i++){
      auto col = cols.GetPortalConstControl().Get(i);
      cols.GetPortalControl().Set(i, col+sumtotl.GetPortalConstControl().Get(i));
    }
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

