#ifndef SURFACEWORKLETS_H
#define SURFACEWORKLETS_H
#include <vtkm/worklet/WorkletMapField.h>

class XYRectWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  XYRectWorklet(
           int cs,
           int d,
      bool _f = false,
      float _tmin = 0.001,
      float _tmax = std::numeric_limits<float>::max()
      )
    :canvasSize(cs)
    ,depth(d)
    ,tmin(_tmin)
    ,tmax(_tmax)
    , flip(_f)
  {
  }


  VTKM_EXEC
  bool hit(
          const ray& r,
          hit_record& rec,
          float x0,
          float x1,
          float y0,
          float y1,
          float k,
          int matId,
          int texId) const
  {
      float t = (k-r.origin()[2]) / r.direction()[2];
      if (t < tmin || t > tmax){
          return false;
      }
      float x = r.origin()[0] + t*r.direction()[0];
      float y = r.origin()[1] + t*r.direction()[1];
      if (x < x0 || x > x1 || y < y0 || y > y1){
        return false;
      }
      rec.u = (x-x0)/(x1-x0);
      rec.v = (y-y0)/(y1-y0);
      rec.t = t;
      rec.matId = matId;
      rec.texId = texId;

      rec.p = r.point_at_parameter(t);
      rec.normal = vec3(0, 0, 1);
      if (flip)
        rec.normal = -rec.normal;
      return true;
  }

  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>
  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12);

  VTKM_EXEC
  template<
          typename MatTypeArrayType,
          typename ColTypeArrayType>
  void operator()(vtkm::Id idx,
                  ray &ray_io,
                  hit_record &hrec,
                  vtkm::Int8 &scattered,
                  vtkm::Int8 &rayHit,
                  vec3 pt1,
                  vec3 pt2,
                  int matIdx,
                  int texIdx,
                  MatTypeArrayType matType,
                  ColTypeArrayType texType
                ) const
  {
    if (scattered){
      auto matId = matType.Get(matIdx);
      auto texId = texType.Get(texIdx);
      float x0 = pt1[0];
      float x1 = pt2[0];
      float y0 = pt1[1];
      float y1 = pt2[1];
      float k = pt1[2];

      rayHit |= hit(ray_io, hrec,
                         x0,x1,y0,y1,k,matId,texId);
    }
  }

  int canvasSize;
  int depth;
  float tmin, tmax;
  bool flip;
};

class XZRectWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  XZRectWorklet(
           int cs,
           int d,
      bool _f = false,
      float _tmin = 0.001,
      float _tmax = std::numeric_limits<float>::max()
      )
    :canvasSize(cs)
    ,depth(d)
    ,tmin(_tmin)
    ,tmax(_tmax)
    , flip(_f)
  {
  }
  VTKM_EXEC
  bool hit(const ray& r,
                  hit_record& rec,
                  float x0, float x1, float z0, float z1,
                  float k,
                  int matId,
                  int texId) const
  {
    float t = (k-r.origin()[1]) / r.direction()[1];
    if (t < tmin || t > tmax){
        return false;
    }
    float x = r.origin()[0] + t*r.direction()[0];
    float z = r.origin()[2] + t*r.direction()[2];
    if (x < x0 || x > x1 || z < z0 || z > z1){
        return false;
    }
    rec.u = (x-x0)/(x1-x0);
    rec.v = (z-z0)/(z1-z0);
    rec.t = t;
    rec.texId = texId;
    rec.matId = matId;

    rec.p = r.point_at_parameter(t);
    rec.normal = vec3(0, 1, 0);
    if (flip)
      rec.normal = -rec.normal;

    return true;
  }
  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>
  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12);

  VTKM_EXEC
  template<
          typename MatTypeArrayType,
          typename ColTypeArrayType>
  void operator()(vtkm::Id idx,
                  ray &ray_io,
                  hit_record &hrec,
                  vtkm::Int8 &scattered,
                  vtkm::Int8 &rayHit,
                  vec3 &pt1,
                  vec3 &pt2,
                  int matIdx,
                  int texIdx,
                  MatTypeArrayType matType,
                  ColTypeArrayType texType
                  ) const
  {
    if (scattered){
      auto matId = matType.Get(matIdx);
      auto texId = texType.Get(texIdx);
      float x0 = pt1[0];
      float x1 = pt2[0];
      float z0 = pt1[2];
      float z1 = pt2[2];
      float k = pt1[1];

      rayHit |= hit(ray_io, hrec,
                         x0,x1,z0,z1,k,matId,texId);
    }
  }

  int canvasSize;
  int depth;
  float tmin, tmax;
  bool flip;
};
class YZRectWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  YZRectWorklet(
           int cs,
           int d,
      bool _f = false,
      float _tmin = 0.001,
      float _tmax = std::numeric_limits<float>::max()
      )
    :canvasSize(cs)
    ,depth(d)
    ,tmin(_tmin)
    ,tmax(_tmax)
    , flip(_f)
  {
  }


  VTKM_EXEC
  bool hit(const ray& r,
            hit_record& rec,
            float y0, float y1, float z0, float z1,
            float k,
            int matId,
            int texId) const
  {
    float t = (k-r.origin()[0]) / r.direction()[0];
    if (t < tmin || t > tmax){
      return false;
    }

    float y = r.origin()[1] + t*r.direction()[1];
    float z = r.origin()[2] + t*r.direction()[2];
    if (y < y0 || y > y1 || z < z0 || z > z1){
      return false;
    }

    rec.u = (y-y0)/(y1-y0);
    rec.v = (z-z0)/(z1-z0);
    rec.t = t;
    rec.texId = texId;
    rec.matId = matId;

    rec.p = r.point_at_parameter(t);
    rec.normal = vec3(1, 0, 0);
    if (flip)
      rec.normal = -rec.normal;
    return true;
  }
  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>
  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12);

  VTKM_EXEC
  template<
          typename MatTypeArrayType,
          typename ColTypeArrayType>
  void operator()(vtkm::Id idx,
                  ray &ray_io,
                  hit_record &hrec,
                  vtkm::Int8 &scattered,
                  vtkm::Int8 &rayHit,
                  vec3 &pt1,
                  vec3 &pt2,
                  int matIdx,
                  int texIdx,
                  MatTypeArrayType matType,
                  ColTypeArrayType texType
                  ) const
  {
    auto matId = matType.Get(matIdx);
    auto texId = texType.Get(texIdx);
    if (scattered){
      float y0 = pt1[1];
      float y1 = pt2[1];
      float z0 = pt1[2];
      float z1 = pt2[2];
      float k = pt1[0];
      rayHit |= hit(ray_io, hrec,
                         y0,y1,z0,z1,k,matId,texId);
    }
  }

  int canvasSize;
  int depth;
  float tmin, tmax;
  bool flip;
};

class SphereIntersecttWorklet : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  SphereIntersecttWorklet(
           int cs,
           int d,
      bool _f = false,
      float _tmin = 0.001,
      float _tmax = std::numeric_limits<float>::max()
      )
    :canvasSize(cs)
    ,depth(d)
    ,tmin(_tmin)
    ,tmax(_tmax)
    , flip(_f)
  {
  }
  VTKM_EXEC
  bool hit(const ray& r,  hit_record& rec,
           vec3 center, float radius,
           int matId, int texId) const {
      vec3 oc = r.origin() - center;
      float a = dot(r.direction(), r.direction());
      float b = dot(oc, r.direction());
      float c = dot(oc, oc) - radius*radius;
      float discriminant = b*b - a*c;
      if (discriminant > 0) {
          float temp = (-b - sqrt(b*b-a*c))/a;
          if (temp < tmin && temp > tmin) {
              rec.t = temp;
              rec.p = r.point_at_parameter(rec.t);
              get_sphere_uv((rec.p-center)/radius, rec.u, rec.v);
              rec.normal = (rec.p - center) / radius;
              rec.matId = matId;
              rec.texId = texId;

              return true;
          }
          temp = (-b + sqrt(b*b-a*c))/a;
          if (temp < tmax && temp > tmin) {
              rec.t = temp;
              rec.p = r.point_at_parameter(rec.t);
              get_sphere_uv((rec.p-center)/radius, rec.u, rec.v);
              rec.normal = (rec.p - center) / radius;
              rec.matId = matId;
              rec.texId = texId;

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
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>
  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12);

  VTKM_EXEC
  template<
          typename MatTypeArrayType,
          typename ColTypeArrayType>
  void operator()(vtkm::Id idx,
                  ray &ray_io,
                  hit_record &hrec,
                  vtkm::Int8 &scattered,
                  vtkm::Int8 &rayHit,
                  vec3 &pt1,
                  vec3 &pt2,
                  int matIdx,
                  int texIdx,
                  MatTypeArrayType matType,
                  ColTypeArrayType texType
                  ) const
  {
    if (scattered){
      auto matId = matType.Get(matIdx);
      auto texId = texType.Get(texIdx);
      float y0 = pt1[1];
      float y1 = pt2[1];
      float z0 = pt1[2];
      float z1 = pt2[2];
      float k = pt1[0];
      rayHit |= hit(ray_io, hrec, pt1, pt2[0],matId,texId);
    }
  }
  int canvasSize;
  int depth;
  float tmin, tmax;
  bool flip;
};
#endif
