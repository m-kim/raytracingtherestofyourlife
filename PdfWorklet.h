#ifndef PDFWORKLET_H
#define PDFWORKLET_H
#include "Surface.h"
#include "Record.h"
#include "wangXor.h"
class GenerateDir : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_EXEC_CONT
  GenerateDir(int ts) : type_size(ts){}

  using ControlSignature = void(FieldInOut<>, FieldInOut<>);
  using ExecutionSignature = void( _1, _2);

  VTKM_EXEC
  void operator()(unsigned int &seed, int &which) const {
    which = vtkm::Min(type_size,int(xorshiftWang::getRandF(seed) * type_size+1));
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
  inline vec3 de_nan(const vec3& c) const {
      vec3 temp = c;
      if (!(temp[0] == temp[0])) temp[0] = 0;
      if (!(temp[1] == temp[1])) temp[1] = 0;
      if (!(temp[2] == temp[2])) temp[2] = 0;
      return temp;
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
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6);

  VTKM_EXEC
  template<typename PtArrayType, typename HitRecord>
  void operator()(
      int which,
       HitRecord &hrec,
       vec3 &generated,
      unsigned int &seed,
       PtArrayType pt1,
       PtArrayType pt2
    ) const
  {
    if (which <= current ){
      float r1 =xorshiftWang::getRandF(seed);
      float r2 =xorshiftWang::getRandF(seed);
      onb uvw;

      vec3 hrecn(hrec[static_cast<vtkm::Id>(HR::Nx)], hrec[static_cast<vtkm::Id>(HR::Ny)], hrec[static_cast<vtkm::Id>(HR::Nz)]);
      uvw.build_from_w(hrecn);
      generated = de_nan(uvw.local(random_cosine_direction(r1,r2)));
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
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6);

  VTKM_EXEC
  template<typename PtArrayType, typename HitRecord>
  void operator()(int which,
           HitRecord &hrec,
           vec3 &generated,
           unsigned int &seed,
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
        vec3 p(hrec[static_cast<vtkm::Id>(HR::Px)], hrec[static_cast<vtkm::Id>(HR::Py)], hrec[static_cast<vtkm::Id>(HR::Pz)]);
        generated = random(p, xorshiftWang::getRandF(seed), xorshiftWang::getRandF(seed),
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
  inline vec3 de_nan(const vec3& c) const {
      vec3 temp = c;
      if (!(temp[0] == temp[0])) temp[0] = 0;
      if (!(temp[1] == temp[1])) temp[1] = 0;
      if (!(temp[2] == temp[2])) temp[2] = 0;
      return temp;
  }

  VTKM_EXEC
  vec3 random_to_sphere(float radius, float distance_squared, float r1, float r2) const {
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
       return de_nan(uvw.local(random_to_sphere(radius, distance_squared, r1,r2)));
  }
  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6);

  VTKM_EXEC
  template<typename PtArrayType, typename HitRecord>
  void operator()(vtkm::Id idx,
          int &which,
           HitRecord &hrec,
           vec3 &generated,
                  unsigned int &seed,
           PtArrayType pt1,
           PtArrayType pt2
           ) const
  {
    if (which == current){
      for (int i = 0; i < pt1.GetNumberOfValues(); i++){
        float radius = pt2.Get(i)[0];
        vec3 p(hrec[static_cast<vtkm::Id>(HR::Px)], hrec[static_cast<vtkm::Id>(HR::Py)], hrec[static_cast<vtkm::Id>(HR::Pz)]);
        generated = random(p, xorshiftWang::getRandF(seed), xorshiftWang::getRandF(seed), pt1.Get(i), radius);
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
    vtkm::Vec<vtkm::Float32, 9> rec;
    vtkm::Vec<vtkm::Int8,2> hid;
    int matId,texId;
    if (surf.hit(o,v, rec, hid, 0.001, std::numeric_limits<float>::max(),x0,x1,z0,z1,k, matId, texId)) {
        float area = (x1-x0)*(z1-z0);
        auto rect = rec[static_cast<vtkm::Id>(HR::T)];
        float distance_squared = rect * rect * vtkm::MagnitudeSquared(v);
        vec3 n(rec[static_cast<vtkm::Id>(HR::Nx)],rec[static_cast<vtkm::Id>(HR::Ny)],rec[static_cast<vtkm::Id>(HR::Nz)]);
        float cosine = fabs(dot(v, n) * vtkm::RMagnitude(v));
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
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9);

  VTKM_EXEC
  template<typename PtArrayType,
          typename FlippedType,
  typename HitRecord,
  int ScatterBitIndex = 3>

  void operator()(vtkm::Id idx,
                  vec3 &origin,
                  vec3 &direction,
                  HitRecord &hrec,
                  FlippedType &scattered,
                  float &sum_value,
                  vec3 &generated,
                  unsigned int &seed,
                  PtArrayType pt1,
                  PtArrayType pt2
                  ) const
  {
    if (scattered & (1UL << ScatterBitIndex)){
      float weight = 1.0/list_size;
      int index = int(xorshiftWang::getRandF(seed) * list_size);
      for (int i = 0; i < pt1.GetNumberOfValues(); i++){
        float x0 = pt1.Get(i)[0];
        float x1 = pt2.Get(i)[0];
        float z0 = pt1.Get(i)[2];
        float z1 = pt2.Get(i)[2];
        float k = pt1.Get(i)[1];
        vec3 p(hrec[static_cast<vtkm::Id>(HR::Px)], hrec[static_cast<vtkm::Id>(HR::Py)], hrec[static_cast<vtkm::Id>(HR::Pz)]);
        sum_value += weight*pdf_value(p, generated, x0,x1,z0,z1,k);
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

    vtkm::Vec<vtkm::Float32, 9> rec;
    vtkm::Vec<vtkm::Id,2> hid;
    int matId,texId;
    if (surf.hit(o, v, rec, hid, 0.001, std::numeric_limits<float>::max(), center, radius, matId, texId )) {
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
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9);

  VTKM_EXEC
  template<typename PtArrayType,
          typename FlippedType,
  typename HitRecord,
  int ScatterBitIndex = 3>
  void operator()(vtkm::Id idx,
                  vec3 &origin,
                  vec3 &direction,
                  HitRecord &hrec,
                  FlippedType &scattered,
                  float &sum_value,
                  vec3 &generated,
                  unsigned int &seed,
                  PtArrayType pt1,
                  PtArrayType pt2
                  ) const
  {
    if (scattered & (1UL << ScatterBitIndex)){
      float weight = 1.0/list_size;
      int index = int(xorshiftWang::getRandF(seed) *list_size);
      for (int i = 0; i < pt1.GetNumberOfValues(); i++){
        float radius = pt2.Get(i)[0];
        vec3 p(hrec[static_cast<vtkm::Id>(HR::Px)], hrec[static_cast<vtkm::Id>(HR::Py)], hrec[static_cast<vtkm::Id>(HR::Pz)]);
        sum_value += weight*pdf_value(p, generated, pt1.Get(i), radius);
      }
    }
  }

  float list_size;
  Sphere surf;
};
#endif
