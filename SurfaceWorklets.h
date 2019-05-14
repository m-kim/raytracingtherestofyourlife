#ifndef SURFACEWORKLETS_H
#define SURFACEWORKLETS_H
#include <vtkm/worklet/WorkletMapField.h>
#include "Surface.h"

auto rotateY(const vec3 &origin, const vec3 &direction, float angle){
  float radians = (M_PI / 180.) * angle;
  float sin_theta = sin(radians);
  float cos_theta = cos(radians);
  vec3 o, d;
  o = origin;
  d = direction;
  o[0] = cos_theta*origin[0] - sin_theta*origin[2];
  o[2] =  sin_theta*origin[0] + cos_theta*origin[2];
  d[0] = cos_theta*direction[0] - sin_theta*direction[2];
  d[2] = sin_theta*direction[0] + cos_theta*direction[2];
  return std::make_tuple(o, d);
}


auto rotAndTrans(vec3 &origin, vec3 &direction, vec3 offset, float angle)
{
  return rotateY(origin - offset, direction, angle);
}

void applyRotAndTrans(HitRecord &temp_rec, vec3 offset, float angle)
{
  vec3 p = temp_rec.p;
  vec3 normal = temp_rec.normal;
  float radians = (M_PI / 180.) * angle;
  float sin_theta = sin(radians);
  float cos_theta = cos(radians);

  p[0] = cos_theta*temp_rec.p[0] + sin_theta*temp_rec.p[2];
  p[2] = -sin_theta*temp_rec.p[0] + cos_theta*temp_rec.p[2];
  normal[0] = cos_theta*temp_rec.normal[0] + sin_theta*temp_rec.normal[2];
  normal[2] = -sin_theta*temp_rec.normal[0] + cos_theta*temp_rec.normal[2];
  temp_rec.p = p;
  temp_rec.normal = normal;
  temp_rec.p += offset;

}
class XYRectWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  XYRectWorklet(
           int cs,
           int d
      )
    :canvasSize(cs)
    ,depth(d)
  {
  }



  VTKM_EXEC
  bool hit(
          const vec3 &origin,
          const vec3 &direction,
          HitRecord& rec,
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
    rec.u = (x-x0)/(x1-x0);
    rec.v = (y-y0)/(y1-y0);
    rec.t = t;
    rec.matId = matId;
    rec.texId = texId;

    rec.p = origin + direction * (t);
    rec.normal = vec3(0, 0, 1);
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

  VTKM_EXEC
  template<typename PtArrayType,
          typename IndexType,
          typename FlippedType,
          typename AngleArray>
  void operator()(vtkm::Id idx,
                  vec3 &origin,
                  vec3 &direction,
                  HitRecord &hrec,
                  float &tmin,
                  float &tmax,
                  vtkm::Int8 &scattered,
                  vtkm::Int8 &rayHit,
                  PtArrayType pt1,
                  PtArrayType pt2,
                  IndexType matIdx,
                  IndexType texIdx,
                  PtArrayType offsetArray,
                  AngleArray angleArray,
                  FlippedType flipped
                  ) const
  {
    if (scattered){
      for (int i=0; i<matIdx.GetNumberOfValues(); i++){
        float x0 = pt1.Get(i)[0];
        float x1 = pt2.Get(i)[0];
        float y0 = pt1.Get(i)[1];
        float y1 = pt2.Get(i)[1];
        float k = pt1.Get(i)[2];
        HitRecord  temp_rec;

        auto offset = offsetArray.Get(i);
        auto angle = angleArray.Get(i);
        auto moved_r = rotAndTrans(origin, direction, offset, angle);

        auto h =  hit(std::get<0>(moved_r), std::get<1>(moved_r), temp_rec, tmin, tmax,
                      x0,x1,y0,y1,k,matIdx.Get(i),texIdx.Get(i));
        if (h){
          if (flipped.Get(i))
            temp_rec.normal = -temp_rec.normal;
          applyRotAndTrans(temp_rec, offset, angle);
          hrec= temp_rec;
          tmax = temp_rec.t;
        }
        rayHit |= h;

      }
    }
  }
  int canvasSize;
  int depth;
};

class XZRectWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  XZRectWorklet(
           int cs,
           int d
      )
    :canvasSize(cs)
    ,depth(d)
  {
  }
  VTKM_EXEC
  bool hit(const vec3 &origin,
           const vec3 &direction,
                  HitRecord& rec,
                  float tmin, float tmax,
                  float x0, float x1, float z0, float z1,
                  float k,
                  int matId,
                  int texId) const
  {
//    float t = (k-r.origin()[1]) / r.direction()[1];
//    if (t < tmin || t > tmax)
//        return false;
//    float x = r.origin()[0] + t*r.direction()[0];
//    float z = r.origin()[2] + t*r.direction()[2];
//    if (x < x0 || x > x1 || z < z0 || z > z1)
//        return false;
//    rec.u = (x-x0)/(x1-x0);
//    rec.v = (z-z0)/(z1-z0);
//    rec.t = t;
//    rec.texId = texId;
//    rec.matId = matId;

//    rec.p = r.point_at_parameter(t);
//    rec.normal = vec3(0, 1, 0);
//    return true;
    return surf.hit(origin, direction,rec,tmin,tmax,x0,x1,z0,z1,k,matId,texId);


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

  VTKM_EXEC
  template<typename PtArrayType,
          typename IndexType,
          typename FlippedType,
          typename AngleArray>
  void operator()(vtkm::Id idx,
                  vec3 &origin,
                  vec3 &direction,
                  HitRecord &hrec,
                  float &tmin,
                  float &tmax,
                  vtkm::Int8 &scattered,
                  vtkm::Int8 &rayHit,
                  PtArrayType pt1,
                  PtArrayType pt2,
                  IndexType matIdx,
                  IndexType texIdx,
                  PtArrayType offsetArray,
                  AngleArray angleArray,
                  FlippedType flipped
                  ) const
  {
    if (scattered){
      for (int i=0; i<matIdx.GetNumberOfValues(); i++){
        float x0 = pt1.Get(i)[0];
        float x1 = pt2.Get(i)[0];
        float z0 = pt1.Get(i)[2];
        float z1 = pt2.Get(i)[2];
        float k = pt1.Get(i)[1];
        HitRecord  temp_rec;

        auto offset = offsetArray.Get(i);
        auto angle = angleArray.Get(i);
        auto moved_r = rotAndTrans(origin, direction, offset, angle);

        auto h =  hit(std::get<0>(moved_r), std::get<1>(moved_r), temp_rec, tmin, tmax,
                      x0,x1,z0,z1,k,matIdx.Get(i),texIdx.Get(i));
        if (h){
          tmax = temp_rec.t;
          if (flipped.Get(i))
            temp_rec.normal = -temp_rec.normal;

          applyRotAndTrans(temp_rec, offset, angle);
          hrec= temp_rec;
        }
        rayHit |= h;
      }
    }
  }

  int canvasSize;
  int depth;
  XZRect surf;
};
class YZRectWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  YZRectWorklet(
           int cs,
           int d
      )
    :canvasSize(cs)
    ,depth(d)
  {
  }


  VTKM_EXEC
  bool hit(const vec3 &origin,
           const vec3 &direction,
            HitRecord& rec,
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

    rec.u = (y-y0)/(y1-y0);
    rec.v = (z-z0)/(z1-z0);
    rec.t = t;
    rec.texId = texId;
    rec.matId = matId;

    rec.p = origin + direction * (t);
    rec.normal = vec3(1, 0, 0);
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

  VTKM_EXEC
  template<typename PtArrayType,
          typename IndexType,
          typename FlippedType,
          typename AngleArray>
  void operator()(vtkm::Id idx,
                  vec3 &origin,
                  vec3 &direction,
                  HitRecord &hrec,
                  float &tmin,
                  float &tmax,
                  vtkm::Int8 &scattered,
                  vtkm::Int8 &rayHit,
                  PtArrayType pt1,
                  PtArrayType pt2,
                  IndexType matIdx,
                  IndexType texIdx,
                  PtArrayType offsetArray,
                  AngleArray angleArray,
                  FlippedType flipped
                  ) const
  {
    if (scattered){
      for (int i=0; i<pt1.GetNumberOfValues(); i++){
        float y0 = pt1.Get(i)[1];
        float y1 = pt2.Get(i)[1];
        float z0 = pt1.Get(i)[2];
        float z1 = pt2.Get(i)[2];
        float k = pt1.Get(i)[0];

        HitRecord temp_rec;

        auto offset = offsetArray.Get(i);
        auto angle = angleArray.Get(i);
        auto moved_r = rotAndTrans(origin, direction, offset, angle);
        auto h =  hit(std::get<0>(moved_r), std::get<1>(moved_r), temp_rec, tmin, tmax,
                      y0,y1,z0,z1,k,matIdx.Get(i),texIdx.Get(i));
        if (h){
          if (flipped.Get(i))
            temp_rec.normal = -temp_rec.normal;
          tmax = temp_rec.t;


          applyRotAndTrans(temp_rec, offset, angle);
          hrec= temp_rec;
        }
        rayHit |= h;

      }
    }
  }

  int canvasSize;
  int depth;
};

class SphereIntersecttWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  SphereIntersecttWorklet(
           int cs,
           int d)
    :canvasSize(cs)
    ,depth(d)
  {
  }
  VTKM_EXEC
  bool hit(const vec3 &origin, const vec3 &direction,  HitRecord& rec, float tmin, float tmax,
           vec3 center, float radius,
           int matId, int texId) const {
      surf.hit(origin, direction,rec,tmin, tmax, center, radius, matId, texId);
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

  VTKM_EXEC
  template<typename PtArrayType,
            typename IndexType>
  void operator()(vtkm::Id idx,
                  vec3 &origin,
                  vec3 &direction,
                  HitRecord &hrec,
                  float &tmin,
                  float &tmax,
                  vtkm::Int8 &scattered,
                  vtkm::Int8 &rayHit,
                  PtArrayType pt1,
                  PtArrayType pt2,
                  IndexType matIdx,
                  IndexType texIdx
                  ) const
  {
    if (scattered){
      for (int i=0; i<pt1.GetNumberOfValues(); i++){
        HitRecord  temp_rec;
        auto h =   hit(origin, direction, temp_rec, tmin, tmax,
                       pt1.Get(i), pt2.Get(i)[0],matIdx.Get(i),texIdx.Get(i));
        if (h){

          tmax = temp_rec.t;
          hrec = temp_rec;
        }
        rayHit |= h;
      }
    }
  }

  Sphere surf;
  int canvasSize;
  int depth;
};


class CollectIntersecttWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  CollectIntersecttWorklet(
           int cs,
           int d)
    :canvasSize(cs)
    ,depth(d)
  {
  }


  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>
  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4);

  VTKM_EXEC
  template<typename PtArrayType>
  void operator()(vtkm::Id idx,
                  vtkm::Int8 &sctr,
                  vtkm::Int8 &hit,
                  PtArrayType &emitted,
                  PtArrayType &attenuation
                  ) const
  {
    if (!(sctr && hit)){
      sctr = false;
      attenuation.Set(idx + canvasSize * depth, vec3(1.0));
      emitted.Set(idx + canvasSize * depth, vec3(0.0f));
    }

    //scattered.GetPortalControl().Set(i,sctr);
  }

  int depth, canvasSize;
};
#endif
