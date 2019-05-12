#ifndef SURFACEWORKLETS_H
#define SURFACEWORKLETS_H
#include <vtkm/worklet/WorkletMapField.h>

ray rotateY(const ray &r, float angle){
  float radians = (M_PI / 180.) * angle;
  float sin_theta = sin(radians);
  float cos_theta = cos(radians);

  vec3 origin = r.origin();
  vec3 direction = r.direction();
  origin[0] = cos_theta*r.origin()[0] - sin_theta*r.origin()[2];
  origin[2] =  sin_theta*r.origin()[0] + cos_theta*r.origin()[2];
  direction[0] = cos_theta*r.direction()[0] - sin_theta*r.direction()[2];
  direction[2] = sin_theta*r.direction()[0] + cos_theta*r.direction()[2];
  ray rotated_r(origin, direction, r.time());
  return rotated_r;
}


ray rotAndTrans(ray &ray_io, vec3 offset, float angle)
{
  auto rot_r = rotateY(ray_io, angle);
  return ray(rot_r.origin() - offset, rot_r.direction(), rot_r.time());
}

void applyRotAndTrans(hit_record &temp_rec, vec3 offset, float angle)
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
          const ray& r,
          hit_record& rec,
          float &tmin, float &tmax,

          float x0,
          float x1,
          float y0,
          float y1,
          float k,
          int matId,
          int texId,
          vtkm::Int8 flip) const
  {
    float t = (k-r.origin()[2]) / r.direction()[2];
    if (t < tmin || t > tmax)
        return false;
    float x = r.origin()[0] + t*r.direction()[0];
    float y = r.origin()[1] + t*r.direction()[1];
    if (x < x0 || x > x1 || y < y0 || y > y1)
        return false;
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
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13);

  VTKM_EXEC
  template<typename PtArrayType,
          typename IndexType,
          typename FlippedType,
          typename AngleArray>
  void operator()(vtkm::Id idx,
                  ray &ray_io,
                  hit_record &hrec,
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
        hit_record  temp_rec;

        auto offset = offsetArray.Get(i);
        auto angle = angleArray.Get(i);
        auto moved_r = rotAndTrans(ray_io, offset, angle);
        auto h =  hit(moved_r, temp_rec, tmin, tmax,
                      x0,x1,y0,y1,k,matIdx.Get(i),texIdx.Get(i),flipped.Get(i));
        if (h){
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
  bool hit(const ray& r,
                  hit_record& rec,
                  float &tmin, float &tmax,
                  float x0, float x1, float z0, float z1,
                  float k,
                  int matId,
                  int texId,
           vtkm::Int8 flip) const
  {
    float t = (k-r.origin()[1]) / r.direction()[1];
    if (t < tmin || t > tmax)
        return false;
    float x = r.origin()[0] + t*r.direction()[0];
    float z = r.origin()[2] + t*r.direction()[2];
    if (x < x0 || x > x1 || z < z0 || z > z1)
        return false;
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
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>
  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11);

  VTKM_EXEC
  template<typename PtArrayType,
          typename IndexType,
          typename FlippedType>
  void operator()(vtkm::Id idx,
                  ray &ray_io,
                  hit_record &hrec,
                  float &tmin,
                  float &tmax,
                  vtkm::Int8 &scattered,
                  vtkm::Int8 &rayHit,
                  PtArrayType pt1,
                  PtArrayType pt2,
                  IndexType matIdx,
                  IndexType texIdx,
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
        hit_record  temp_rec;
        auto h =  hit(ray_io, temp_rec, tmin, tmax,
                      x0,x1,z0,z1,k,matIdx.Get(i),texIdx.Get(i),flipped.Get(i));
        if (h){

          tmax = temp_rec.t;
          hrec = temp_rec;
        }
        rayHit |= h;

      }
    }
  }

  int canvasSize;
  int depth;
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
  bool hit(const ray& r,
            hit_record& rec,
            float &tmin, float &tmax,
            float y0, float y1, float z0, float z1,
            float k,
            int matId,
            int texId,
           vtkm::Int8 flip) const
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
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>

  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13);

  VTKM_EXEC
  template<typename PtArrayType,
          typename IndexType,
          typename FlippedType,
          typename AngleArray>
  void operator()(vtkm::Id idx,
                  ray &ray_io,
                  hit_record &hrec,
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

        hit_record temp_rec;

        auto offset = offsetArray.Get(i);
        auto angle = angleArray.Get(i);
        auto moved_r = rotAndTrans(ray_io, offset, angle);
        auto h =  hit(moved_r, temp_rec, tmin, tmax,
                      y0,y1,z0,z1,k,matIdx.Get(i),texIdx.Get(i), flipped.Get(i));
        if (h){
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
  bool hit(const ray& r,  hit_record& rec, float& tmin, float &tmax,
           vec3 center, float radius,
           int matId, int texId) const {
      vec3 oc = r.origin() - center;
      float a = dot(r.direction(), r.direction());
      float b = dot(oc, r.direction());
      float c = dot(oc, oc) - radius*radius;
      float discriminant = b*b - a*c;
      if (discriminant > 0) {
          float temp = (-b - sqrt(b*b-a*c))/a;
          if (temp < tmax && temp > tmin) {
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
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>
  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10);

  VTKM_EXEC
  template<typename PtArrayType,
            typename IndexType>
  void operator()(vtkm::Id idx,
                  ray &ray_io,
                  hit_record &hrec,
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
        hit_record  temp_rec;
        auto h =   hit(ray_io, temp_rec, tmin, tmax,
                       pt1.Get(i), pt2.Get(i)[0],matIdx.Get(i),texIdx.Get(i));
        if (h){

          tmax = temp_rec.t;
          hrec = temp_rec;
        }
        rayHit |= h;
      }
    }
  }
  int canvasSize;
  int depth;
};
#endif
