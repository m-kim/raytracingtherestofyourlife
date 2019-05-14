#ifndef SURFACE_H
#define SURFACE_H
#include "Record.h"

class XZRect
{
public:
  VTKM_EXEC_CONT XZRect(){}
  VTKM_EXEC
  bool hit(const vec3& origin,
           const vec3& dir,
                  HitRecord& rec,
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
    rec.u = (x-x0)/(x1-x0);
    rec.v = (z-z0)/(z1-z0);
    rec.t = t;
    rec.texId = texId;
    rec.matId = matId;

    rec.p = origin + dir * t;

    rec.normal = vec3(0, 1, 0);


    return true;
  }
};


class Sphere
{
public:
  VTKM_EXEC_CONT Sphere(){}

  VTKM_EXEC
  void get_sphere_uv(const vec3& p, float& u, float& v) const {
      float phi = atan2(p[2], p[0]);
      float theta = asin(p[1]);
      u = 1-(phi + M_PI) / (2*M_PI);
      v = (theta + M_PI/2) / M_PI;
  }

  VTKM_EXEC
  bool hit(const vec3& origin, const vec3 &direction,
           HitRecord& rec, float tmin, float tmax,
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
              rec.t = temp;
              rec.p = origin + direction * rec.t;
              get_sphere_uv((rec.p-center)/radius, rec.u, rec.v);
              rec.normal = (rec.p - center) / radius;
              rec.matId = matId;
              rec.texId = texId;

              return true;
          }
          temp = (-b + sqrt(b*b-a*c))/a;
          if (temp < tmax && temp > tmin) {
              rec.t = temp;
              rec.p = origin + direction * (rec.t);
              get_sphere_uv((rec.p-center)/radius, rec.u, rec.v);
              rec.normal = (rec.p - center) / radius;
              rec.matId = matId;
              rec.texId = texId;

              return true;
          }
      }
      return false;
  }
};
#endif
