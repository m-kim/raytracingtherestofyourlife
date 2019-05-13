#ifndef SURFACE_H
#define SURFACE_H

class XZRect
{
public:
  VTKM_EXEC_CONT XZRect(){}
  VTKM_EXEC
  bool hit(const ray& r,
                  hit_record& rec,
                  float tmin, float tmax,
                  float x0, float x1, float z0, float z1,
                  float k,
                  int matId,
                  int texId) const
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


    return true;
  }
};


class Sphere
{
public:
  VTKM_EXEC_CONT Sphere(){}
  VTKM_EXEC
  bool hit(const ray& r,  hit_record& rec, float tmin, float tmax,
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
};
#endif
