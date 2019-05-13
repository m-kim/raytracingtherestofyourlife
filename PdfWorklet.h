#ifndef PDFWORKLET_H
#define PDFWORKLET_H
#include "Surface.h"
class GenerateDir : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_EXEC_CONT
  GenerateDir(int ts) : type_size(ts){}

  using ControlSignature = void(FieldInOut<>);
  using ExecutionSignature = void( _1);

  VTKM_EXEC
  void operator()(int &which) const {
    which = vtkm::Min(type_size,int(drand48() * type_size+1));
  }
  const int type_size;
};

class CosineGenerateDir : public vtkm::worklet::WorkletMapField
{
public:
  CosineGenerateDir(int cur = 0) : current(cur){}

  VTKM_EXEC
  vec3 random(const vec3& o, float r1, float r2,
              float x0, float x1, float z0, float z1, float k) const {
      vec3 random_point = vec3(x0 + r1*(x1-x0), k,  z0 + r2*(z1-z0));
      return random_point - o;
  }

  VTKM_EXEC
  inline vec3 random_cosine_direction(float r1, float r2) const {
      float z = sqrt(1-r2);
      float phi = 2*M_PI*r1;
      float x = cos(phi)*2*sqrt(r2);
      float y = sin(phi)*2*sqrt(r2);
      return vec3(x, y, z);
  }
  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(_1, _2, _3, _4, _5);

  VTKM_EXEC
  template<typename PtArrayType>
  void operator()(
      int which,
       HitRecord &hrec,
       vec3 &generated,
       PtArrayType pt1,
       PtArrayType pt2
    ) const
  {
    if (which <= current ){
      float r1 =drand48();
      float r2 =drand48();
      onb uvw;
      uvw.build_from_w(hrec.normal);
      generated = uvw.local(random_cosine_direction(r1,r2));
    }
  }
  int current;
};

class XZRectGenerateDir : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_EXEC_CONT
  XZRectGenerateDir(int cur = 1): current(cur) {}
  VTKM_EXEC
  vec3 random(const vec3& o, float r1, float r2,
              float x0, float x1, float z0, float z1, float k) const {
      vec3 random_point = vec3(x0 + r1*(x1-x0), k,  z0 + r2*(z1-z0));
      return random_point - o;
  }

  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(_1, _2, _3, _4, _5);

  VTKM_EXEC
  template<typename PtArrayType>
  void operator()(int which,
           HitRecord &hrec,
           vec3 &generated,
           PtArrayType pt1,
           PtArrayType pt2
           ) const
  {
    if (current == which){
      for (int i = 0; i < pt1.GetNumberOfValues(); i++){
        float x0 = pt1.Get(i)[0];
        float x1 = pt2.Get(i)[0];
        float z0 = pt1.Get(i)[2];
        float z1 = pt2.Get(i)[2];
        float k = pt1.Get(i)[1];
        generated = random(hrec.p, drand48(), drand48(),
                           x0,x1,z0,z1,k);
      }
    }
  }

  int current;
};
class SphereGenerateDir : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_EXEC_CONT
  SphereGenerateDir(int cur = 2): current(cur) {}

  VTKM_EXEC
  inline vec3 random_to_sphere(float radius, float distance_squared, float r1, float r2) const {
  //    float r1 = drand48();
  //    float r2 = drand48();
      float z = 1 + r2*(sqrt(1-radius*radius/distance_squared) - 1);
      float phi = 2*M_PI*r1;
      float x = cos(phi)*sqrt(1-z*z);
      float y = sin(phi)*sqrt(1-z*z);
      return vec3(x, y, z);
  }
  VTKM_EXEC
  vec3 random(const vec3& o, float r1, float r2,
              const vec3 &center, float radius) const {
       vec3 direction = center - o;
       float distance_squared = MagnitudeSquared(direction);
       onb uvw;
       uvw.build_from_w(direction);
       return uvw.local(random_to_sphere(radius, distance_squared, r1,r2));
  }
  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5);

  VTKM_EXEC
  template<typename PtArrayType>
  void operator()(vtkm::Id idx,
          int &which,
           HitRecord &hrec,
           vec3 &generated,
           PtArrayType pt1,
           PtArrayType pt2
           ) const
  {
    if (which == current){
      for (int i = 0; i < pt1.GetNumberOfValues(); i++){
        float radius = pt2.Get(i)[0];
        generated = random(hrec.p, drand48(), drand48(), pt1.Get(i), radius);
      }
    }
  }
  int current;
};

class XZRectPDFWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  XZRectPDFWorklet(int ls
      )
    :list_size(ls)
  {
  }
  VTKM_EXEC
  float  pdf_value(const vec3& o, const vec3& v,
                   float x0, float x1, float z0, float z1, float k) const {
    HitRecord rec;
    int matId,texId;
    if (surf.hit(ray(o,v), rec,0.001, std::numeric_limits<float>::max(),x0,x1,z0,z1,k, matId, texId)) {
        float area = (x1-x0)*(z1-z0);
        float distance_squared = rec.t * rec.t * vtkm::MagnitudeSquared(v);
        float cosine = fabs(dot(v, rec.normal) * vtkm::RMagnitude(v));
        return  distance_squared / (cosine * area);
    }
    else
        return 0;
  }


  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7);

  VTKM_EXEC
  template<typename PtArrayType,
          typename FlippedType>
  void operator()(vtkm::Id idx,
                  ray &ray_io,
                  HitRecord &hrec,
                  FlippedType &scattered,
                  float &sum_value,
                  vec3 &generated,
                  PtArrayType pt1,
                  PtArrayType pt2
                  ) const
  {
    if (scattered){
      float weight = 1.0/list_size;
      int index = int(drand48() * list_size);
      for (int i = 0; i < pt1.GetNumberOfValues(); i++){
        float x0 = pt1.Get(i)[0];
        float x1 = pt2.Get(i)[0];
        float z0 = pt1.Get(i)[2];
        float z1 = pt2.Get(i)[2];
        float k = pt1.Get(i)[1];
        sum_value += weight*pdf_value(hrec.p, generated, x0,x1,z0,z1,k);
        //if (!index)

      }
    }
  }

  float list_size;
  XZRect surf;
};

class SpherePDFWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  SpherePDFWorklet(int ls
      )
    : list_size(ls)
  {
  }
  VTKM_EXEC
  float  pdf_value(const vec3& o, const vec3& v,
                  vec3 center, float radius) const {
    HitRecord rec;
    int matId,texId;
    if (surf.hit(ray(o, v), rec, 0.001, std::numeric_limits<float>::max(), center, radius, matId, texId )) {
        float cos_theta_max = sqrt(1 - radius*radius/vtkm::MagnitudeSquared(center-o));
        float solid_angle = 2*M_PI*(1-cos_theta_max);
        return  1 / solid_angle;
    }
    else
        return 0;

  }


  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7);

  VTKM_EXEC
  template<typename PtArrayType,
          typename FlippedType>
  void operator()(vtkm::Id idx,
                  ray &ray_io,
                  HitRecord &hrec,
                  FlippedType &scattered,
                  float &sum_value,
                  vec3 &generated,
                  PtArrayType pt1,
                  PtArrayType pt2
                  ) const
  {
    if (scattered){
      float weight = 1.0/list_size;
      int index = int(drand48() *list_size);
      for (int i = 0; i < pt1.GetNumberOfValues(); i++){
        float radius = pt2.Get(i)[0];
        sum_value += weight*pdf_value(hrec.p, generated, pt1.Get(i), radius);
      }
    }
  }

  float list_size;
  Sphere surf;
};
#endif
