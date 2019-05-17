#ifndef SURFACE_H
#define SURFACE_H
#include "Record.h"

class XZRect
{
public:
  VTKM_EXEC_CONT XZRect(){}
  template<typename HitRecord, typename HitId>
  VTKM_EXEC
  bool hit(const vec3& origin,
           const vec3& dir,
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
};


class Sphere
{
public:
  Sphere(){}

  VTKM_EXEC
  void get_sphere_uv(const vec3& p, float& u, float& v) const {
      float phi = atan2(p[2], p[0]);
      float theta = asin(p[1]);
      u = 1-(phi + M_PI) / (2*M_PI);
      v = (theta + M_PI/2) / M_PI;
  }

  template<typename HitRecord, typename HitId>
  VTKM_EXEC
  bool hit(const vec3& origin, const vec3 &direction,
           HitRecord& rec, HitId &hid, float tmin, float tmax,
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
};


class SphereExecWrapper : public vtkm::cont::ExecutionObjectBase
{

public:
  SphereExecWrapper()
  {
  }

  template <typename Device>
  VTKM_CONT Sphere PrepareForExecution(Device) const
  {
    return Sphere();
  }
};

class XZRectExecWrapper : public vtkm::cont::ExecutionObjectBase
{

public:
  XZRectExecWrapper()
  {
  }

  template <typename Device>
  VTKM_CONT XZRect PrepareForExecution(Device) const
  {
    return XZRect();
  }
};
#endif
