#ifndef SURFACEWORKLETS_H
#define SURFACEWORKLETS_H
#include <vtkm/worklet/WorkletMapField.h>
#include "Surface.h"
VTKM_EXEC
void rotateY(const vec3 &origin, const vec3 &direction, float angle,
             vec3 &o, vec3 &d){
  float radians = (M_PI / 180.) * angle;
  float sin_theta = sin(radians);
  float cos_theta = cos(radians);
  o = origin;
  d = direction;
  o[0] = cos_theta*origin[0] - sin_theta*origin[2];
  o[2] =  sin_theta*origin[0] + cos_theta*origin[2];
  d[0] = cos_theta*direction[0] - sin_theta*direction[2];
  d[2] = sin_theta*direction[0] + cos_theta*direction[2];
}

VTKM_EXEC
void rotAndTrans(vec3 &origin, vec3 &direction, vec3 offset, float angle,
                 vec3 &o, vec3 &d)
{
  return rotateY(origin - offset, direction, angle, o, d);
}

template<typename HitRecord>
VTKM_EXEC
void applyRotAndTrans(HitRecord &temp_rec, vec3 offset, float angle)
{
  vec3 op(temp_rec[static_cast<vtkm::Id>(HR::Px)], temp_rec[static_cast<vtkm::Id>(HR::Py)], temp_rec[static_cast<vtkm::Id>(HR::Pz)]);
  auto p = op;
  vec3 on(temp_rec[static_cast<vtkm::Id>(HR::Nx)], temp_rec[static_cast<vtkm::Id>(HR::Ny)], temp_rec[static_cast<vtkm::Id>(HR::Nz)]);
  auto normal = on;
  float radians = (M_PI / 180.) * angle;
  float sin_theta = sin(radians);
  float cos_theta = cos(radians);

  p[0] = cos_theta*op[0] + sin_theta*op[2];
  p[2] = -sin_theta*op[0] + cos_theta*op[2];
  normal[0] = cos_theta*on[0] + sin_theta*on[2];
  normal[2] = -sin_theta*on[0] + cos_theta*on[2];

  p += offset;
  temp_rec[static_cast<vtkm::Id>(HR::Px)] = p[0];
  temp_rec[static_cast<vtkm::Id>(HR::Py)] = p[1];
  temp_rec[static_cast<vtkm::Id>(HR::Pz)] = p[2];
  temp_rec[static_cast<vtkm::Id>(HR::Nx)] = normal[0];
  temp_rec[static_cast<vtkm::Id>(HR::Ny)] = normal[1];
  temp_rec[static_cast<vtkm::Id>(HR::Nz)] = normal[2];
}
class XYRectWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  XYRectWorklet(
           vtkm::Id cs,
           vtkm::Id d
      )
    :canvasSize(cs)
    ,depth(d)
  {
  }



  template<typename HitRecord, typename HitId>
  VTKM_EXEC
  bool hit(
          const vec3 &origin,
          const vec3 &direction,
          HitRecord& rec,
          HitId &hid,
          float tmin, float tmax,

          float x0,
          float x1,
          float y0,
          float y1,
          float k,
          int matId,
          int texId) const
  {
    float t = (k-origin[2]) / direction[2];
    if (t < tmin || t > tmax)
        return false;
    float x = origin[0] + t*direction[0];
    float y = origin[1] + t*direction[1];
    if (x < x0 || x > x1 || y < y0 || y > y1)
        return false;
    rec[static_cast<vtkm::Id>(HR::U)] = (x-x0)/(x1-x0);
    rec[static_cast<vtkm::Id>(HR::V)] = (y-y0)/(y1-y0);
    rec[static_cast<vtkm::Id>(HR::T)] = t;
    hid[static_cast<vtkm::Id>(HI::M)] = matId;
    hid[static_cast<vtkm::Id>(HI::T)] = texId;

    vec3 p(origin + direction * (t));
    rec[static_cast<vtkm::Id>(HR::Px)] = p[0];
    rec[static_cast<vtkm::Id>(HR::Py)] = p[1];
    rec[static_cast<vtkm::Id>(HR::Pz)] = p[2];
    rec[static_cast<vtkm::Id>(HR::Nx)] = 0;
    rec[static_cast<vtkm::Id>(HR::Ny)] = 0;
    rec[static_cast<vtkm::Id>(HR::Nz)] = 1;
    return true;
  }

  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14);

  template<typename PtArrayType,
          typename IndexType,
          typename FlippedType,
          typename AngleArray,
          typename HitRecord,
          typename HitId,
  int HitBitIdx = 2,
  int ScatterBitIdx = 3>
  VTKM_EXEC
  void operator()(vtkm::Id idx,
                  vec3 &origin,
                  vec3 &direction,
                  HitRecord &hrec,
                  HitId &hid,
                  float &tmin,
                  float &tmax,
                  vtkm::UInt8 &scattered,
                  PtArrayType pt1,
                  PtArrayType pt2,
                  IndexType matIdx,
                  IndexType texIdx,
                  PtArrayType offsetArray,
                  AngleArray angleArray,
                  FlippedType flipped
                  ) const
  {
    if (scattered & (1UL << ScatterBitIdx)){
      for (int i=0; i<matIdx.GetNumberOfValues(); i++){
        float x0 = pt1.Get(i)[0];
        float x1 = pt2.Get(i)[0];
        float y0 = pt1.Get(i)[1];
        float y1 = pt2.Get(i)[1];
        float k = pt1.Get(i)[2];
        HitRecord  temp_rec;
        HitId temp_hid;
        auto offset = offsetArray.Get(i);
        auto angle = angleArray.Get(i);

        vec3 o,d;
        rotAndTrans(origin, direction, offset, angle, o,d);

        auto h =  hit(o, d, temp_rec, temp_hid, tmin, tmax,
                      x0,x1,y0,y1,k,matIdx.Get(i),texIdx.Get(i));
        if (h){
          if (flipped.Get(i)){
            temp_rec[static_cast<vtkm::Id>(HR::Nx)] = -temp_rec[static_cast<vtkm::Id>(HR::Nx)];
            temp_rec[static_cast<vtkm::Id>(HR::Ny)] = -temp_rec[static_cast<vtkm::Id>(HR::Ny)];
            temp_rec[static_cast<vtkm::Id>(HR::Nz)] = -temp_rec[static_cast<vtkm::Id>(HR::Nz)];
          }
          applyRotAndTrans(temp_rec, offset, angle);
          hrec= temp_rec;
          hid = temp_hid;
          tmax = temp_rec[static_cast<vtkm::Id>(HR::T)];
        }
        scattered |= (h << HitBitIdx);

      }
    }
  }
  const vtkm::Id canvasSize;
  const vtkm::Id depth;
};

class XZRectWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  XZRectWorklet(
           vtkm::Id cs,
           vtkm::Id d
      )
    :canvasSize(cs)
    ,depth(d)
  {
  }
  template<typename HitRecord, typename HitId>
  VTKM_EXEC
  bool hit(const vec3 &origin,
           const vec3 &dir,
                  HitRecord& rec,
           HitId &hid,
                  float tmin, float tmax,
                  float x0, float x1, float z0, float z1,
                  float k,
                  int matId,
                  int texId) const
  {
    float t = (k-origin[1]) / dir[1];
    if (t < tmin || t > tmax)
        return false;
    float x = origin[0] + t*dir[0];
    float z = origin[2] + t*dir[2];
    if (x < x0 || x > x1 || z < z0 || z > z1)
        return false;
    rec[static_cast<vtkm::Id>(HR::U)] = (x-x0)/(x1-x0);
    rec[static_cast<vtkm::Id>(HR::V)] = (z-z0)/(z1-z0);
    rec[static_cast<vtkm::Id>(HR::T)] = t;
    hid[static_cast<vtkm::Id>(HI::T)] = texId;
    hid[static_cast<vtkm::Id>(HI::M)] = matId;
    auto p = origin + dir * t;
    rec[static_cast<vtkm::Id>(HR::Px)] = p[0];
    rec[static_cast<vtkm::Id>(HR::Py)] = p[1];
    rec[static_cast<vtkm::Id>(HR::Pz)] = p[2];

    rec[static_cast<vtkm::Id>(HR::Nx)] = 0;
    rec[static_cast<vtkm::Id>(HR::Ny)] = 1;
    rec[static_cast<vtkm::Id>(HR::Nz)] = 0;
    return true;

  }
  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14);

  template<typename PtArrayType,
          typename IndexType,
          typename FlippedType,
          typename AngleArray,
          typename HitRecord,
          typename HitId,
          int HitBitIdx = 2,
          int ScatterBitIdx= 3>
  VTKM_EXEC
  void operator()(vtkm::Id idx,
                  vec3 &origin,
                  vec3 &direction,
                  HitRecord &hrec,
                  HitId &hid,
                  float &tmin,
                  float &tmax,
                  vtkm::UInt8 &scattered,
                  PtArrayType pt1,
                  PtArrayType pt2,
                  IndexType matIdx,
                  IndexType texIdx,
                  PtArrayType offsetArray,
                  AngleArray angleArray,
                  FlippedType flipped
                  ) const
  {
    if (scattered & (1UL << ScatterBitIdx)){ //scattered
      for (int i=0; i<matIdx.GetNumberOfValues(); i++){
        float x0 = pt1.Get(i)[0];
        float x1 = pt2.Get(i)[0];
        float z0 = pt1.Get(i)[2];
        float z1 = pt2.Get(i)[2];
        float k = pt1.Get(i)[1];
        HitRecord  temp_rec;
        HitId temp_hid;
        auto offset = offsetArray.Get(i);
        auto angle = angleArray.Get(i);
        vec3 o,d;
        rotAndTrans(origin, direction, offset, angle, o,d);

        auto h =  hit(o,d, temp_rec, temp_hid, tmin, tmax,
                      x0,x1,z0,z1,k,matIdx.Get(i),texIdx.Get(i));
        if (h){
          tmax = temp_rec[static_cast<vtkm::Id>(HR::T)];
          if (flipped.Get(i)){
            temp_rec[static_cast<vtkm::Id>(HR::Nx)] = -temp_rec[static_cast<vtkm::Id>(HR::Nx)];
            temp_rec[static_cast<vtkm::Id>(HR::Ny)] = -temp_rec[static_cast<vtkm::Id>(HR::Ny)];
            temp_rec[static_cast<vtkm::Id>(HR::Nz)] = -temp_rec[static_cast<vtkm::Id>(HR::Nz)];
          }
          applyRotAndTrans(temp_rec, offset, angle);
          hrec= temp_rec;
          hid = temp_hid;
        }
        scattered |= (h << HitBitIdx);
      }
    }
  }

  vtkm::Id canvasSize;
  vtkm::Id depth;
};
class YZRectWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  YZRectWorklet(
           vtkm::Id cs,
           vtkm::Id d
      )
    :canvasSize(cs)
    ,depth(d)
  {
  }


  template<typename HitRecord, typename HitId>
  VTKM_EXEC
  bool hit(const vec3 &origin,
           const vec3 &direction,
            HitRecord& rec,
           HitId &hid,
            float tmin, float tmax,
            float y0, float y1, float z0, float z1,
            float k,
            int matId,
            int texId) const
  {
    float t = (k-origin[0]) / direction[0];
    if (t < tmin || t > tmax){
      return false;
    }

    float y = origin[1] + t*direction[1];
    float z = origin[2] + t*direction[2];
    if (y < y0 || y > y1 || z < z0 || z > z1){
      return false;
    }

    rec[static_cast<vtkm::Id>(HR::U)] = (y-y0)/(y1-y0);
    rec[static_cast<vtkm::Id>(HR::V)] = (z-z0)/(z1-z0);
    rec[static_cast<vtkm::Id>(HR::T)] = t;
    hid[static_cast<vtkm::Id>(HI::T)] = texId;
    hid[static_cast<vtkm::Id>(HI::M)] = matId;
    auto p = origin + direction * (t);
    rec[static_cast<vtkm::Id>(HR::Px)] = p[0];
    rec[static_cast<vtkm::Id>(HR::Py)] = p[1];
    rec[static_cast<vtkm::Id>(HR::Pz)] = p[2];

    rec[static_cast<vtkm::Id>(HR::Nx)] = 1;
    rec[static_cast<vtkm::Id>(HR::Ny)] = 0;
    rec[static_cast<vtkm::Id>(HR::Nz)] = 0;

    return true;
  }
  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14);

  template<typename PtArrayType,
          typename IndexType,
          typename FlippedType,
          typename AngleArray,
          typename HitRecord,
          typename HitId,
          int HitBitIdx = 2,
          int ScatterBitIdx= 3>
  VTKM_EXEC
  void operator()(vtkm::Id idx,
                  vec3 &origin,
                  vec3 &direction,
                  HitRecord &hrec,
                  HitId &hid,
                  float &tmin,
                  float &tmax,
                  vtkm::UInt8 &scattered,
                  PtArrayType pt1,
                  PtArrayType pt2,
                  IndexType matIdx,
                  IndexType texIdx,
                  PtArrayType offsetArray,
                  AngleArray angleArray,
                  FlippedType flipped
                  ) const
  {
    if (scattered & (1UL << ScatterBitIdx)){
      for (int i=0; i<pt1.GetNumberOfValues(); i++){
        float y0 = pt1.Get(i)[1];
        float y1 = pt2.Get(i)[1];
        float z0 = pt1.Get(i)[2];
        float z1 = pt2.Get(i)[2];
        float k = pt1.Get(i)[0];

        HitRecord temp_rec;
        HitId temp_id;
        auto offset = offsetArray.Get(i);
        auto angle = angleArray.Get(i);
        vec3 o,d;
        rotAndTrans(origin, direction, offset, angle, o,d);
        auto h =  hit(o,d, temp_rec, temp_id,  tmin, tmax,
                      y0,y1,z0,z1,k,matIdx.Get(i),texIdx.Get(i));
        if (h){
          if (flipped.Get(i)){
            temp_rec[static_cast<vtkm::Id>(HR::Nx)] = -temp_rec[static_cast<vtkm::Id>(HR::Nx)];
            temp_rec[static_cast<vtkm::Id>(HR::Ny)] = -temp_rec[static_cast<vtkm::Id>(HR::Ny)];
            temp_rec[static_cast<vtkm::Id>(HR::Nz)] = -temp_rec[static_cast<vtkm::Id>(HR::Nz)];
          }
          tmax = temp_rec[static_cast<vtkm::Id>(HR::T)];


          applyRotAndTrans(temp_rec, offset, angle);
          hrec= temp_rec;
          hid = temp_id;
        }
        scattered |= (h << HitBitIdx);

      }
    }
  }

  vtkm::Id canvasSize;
  vtkm::Id depth;
};

class SphereIntersecttWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  SphereIntersecttWorklet(
           vtkm::Id cs,
           vtkm::Id d)
    :canvasSize(cs)
    ,depth(d)
  {
  }
  VTKM_EXEC
  void get_sphere_uv(const vec3& p, float& u, float& v) const {
      float phi = atan2(p[2], p[0]);
      float theta = asin(p[1]);
      u = 1-(phi + M_PI) / (2*M_PI);
      v = (theta + M_PI/2) / M_PI;
  }

  template<typename HitRecord, typename HitId>
  VTKM_EXEC
  bool hit(const vec3 &origin, const vec3 &direction,  HitRecord& rec, HitId &hid, float &tmin, float &tmax,
           vec3 center, float radius,
           int matId, int texId) const {
    vec3 oc = origin - center;
    float a = dot(direction, direction);
    float b = dot(oc, direction);
    float c = dot(oc, oc) - radius*radius;
    float discriminant = b*b - a*c;
    if (discriminant > 0) {
        float temp = (-b - sqrt(b*b-a*c))/a;
        if (temp < tmax && temp > tmin) {
            rec[static_cast<vtkm::Id>(HR::T)] = temp;
            auto p = origin + direction * rec[static_cast<vtkm::Id>(HR::T)];
            rec[static_cast<vtkm::Id>(HR::Px)] = p[0];
            rec[static_cast<vtkm::Id>(HR::Py)] = p[1];
            rec[static_cast<vtkm::Id>(HR::Pz)] = p[2];
            get_sphere_uv((p-center)/radius, rec[static_cast<vtkm::Id>(HR::U)], rec[static_cast<vtkm::Id>(HR::V)]);
            auto n = (p - center) / radius;
            rec[static_cast<vtkm::Id>(HR::Nx)] = n[0];
            rec[static_cast<vtkm::Id>(HR::Ny)] = n[1];
            rec[static_cast<vtkm::Id>(HR::Nz)] = n[2];
            hid[static_cast<vtkm::Id>(HI::M)] = matId;
            hid[static_cast<vtkm::Id>(HI::T)] = texId;

            return true;
        }
        temp = (-b + sqrt(b*b-a*c))/a;
        if (temp < tmax && temp > tmin) {
            rec[static_cast<vtkm::Id>(HR::T)] = temp;
            auto p = origin + direction * (rec[static_cast<vtkm::Id>(HR::T)]);
            rec[static_cast<vtkm::Id>(HR::Px)] = p[0];
            rec[static_cast<vtkm::Id>(HR::Py)] = p[1];
            rec[static_cast<vtkm::Id>(HR::Pz)] = p[2];
            get_sphere_uv((p-center)/radius, rec[static_cast<vtkm::Id>(HR::U)], rec[static_cast<vtkm::Id>(HR::V)]);
            auto n = (p - center) / radius;
            rec[static_cast<vtkm::Id>(HR::Nx)] = n[0];
            rec[static_cast<vtkm::Id>(HR::Ny)] = n[1];
            rec[static_cast<vtkm::Id>(HR::Nz)] = n[2];
            hid[static_cast<vtkm::Id>(HI::M)] = matId;
            hid[static_cast<vtkm::Id>(HI::T)] = texId;

            return true;
        }
    }
    return false;
  }


  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>
  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11);

  template<typename PtArrayType,
            typename IndexType,
  typename HitRecord,
  typename HitId,
  int HitBitIdx = 2,
  int ScatterBitIdx= 3>
  VTKM_EXEC
  void operator()(vtkm::Id idx,
                  vec3 &origin,
                  vec3 &direction,
                  HitRecord &hrec,
                  HitId &hid,
                  float &tmin,
                  float &tmax,
                  vtkm::UInt8 &scattered,
                  PtArrayType pt1,
                  PtArrayType pt2,
                  IndexType matIdx,
                  IndexType texIdx
                  ) const
  {
    if (scattered & (1UL << ScatterBitIdx)){ //scattered
      for (int i=0; i<pt1.GetNumberOfValues(); i++){
        HitRecord  temp_rec;
        HitId temp_hid;
        auto h =   hit(origin, direction, temp_rec, temp_hid, tmin, tmax,
                       pt1.Get(i), pt2.Get(i)[0],matIdx.Get(i),texIdx.Get(i));
        if (h){

          tmax = temp_rec[static_cast<vtkm::Id>(HR::T)];
          hrec = temp_rec;
          hid = temp_hid;
        }
        scattered |= (h << HitBitIdx);
      }
    }
  }

  Sphere surf;
  vtkm::Id canvasSize;
  vtkm::Id depth;
};


class CollectIntersecttWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  CollectIntersecttWorklet(
           vtkm::Id cs,
           int d)
    :canvasSize(cs)
    ,depth(d)
  {
  }


  using ControlSignature = void(FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>
  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3);

  template<typename PtArrayType,
  int HitBitIdx = 2,
  int ScatterBitIndex = 3>
  VTKM_EXEC
  void operator()(vtkm::Id idx,
                  vtkm::UInt8 &sctr,
                  PtArrayType &emitted,
                  PtArrayType &attenuation
                  ) const
  {
    if (!((sctr & (1UL << ScatterBitIndex)) && (sctr & (1UL << HitBitIdx)))){ //hitRay
      sctr &= ~(1UL << ScatterBitIndex);
      attenuation.Set(idx + canvasSize * depth, vec3(1.0));
      emitted.Set(idx + canvasSize * depth, vec3(0.0f));
    }
    sctr &= ~(1UL << HitBitIdx);
    //scattered.GetPortalControl().Set(i,sctr);
  }

  int depth;
  vtkm::Id canvasSize;
};
#endif
