#ifndef SURFACE_H
#define SURFACE_H
#include "Record.h"
#include "vec3.h"

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

template <typename Device>
class SphereLeafIntersector
{
public:
  using IdHandle = vtkm::cont::ArrayHandle<vtkm::Id>;
  using Id2Handle = vtkm::cont::ArrayHandle<vtkm::Id2>;
  using FloatHandle = vtkm::cont::ArrayHandle<vtkm::Float32>;
  using IdArrayPortal = typename IdHandle::ExecutionTypes<Device>::PortalConst;
  using Id2ArrayPortal = typename Id2Handle::ExecutionTypes<Device>::PortalConst;
  IdArrayPortal PointIds;
  IdArrayPortal MatIdx, TexIdx;

  SphereLeafIntersector(){}
  SphereLeafIntersector(const IdHandle& pointIds,
                        const IdHandle& matIdx,
                        const IdHandle& texIdx)
    : PointIds(pointIds.PrepareForInput(Device()))
    , MatIdx(matIdx.PrepareForInput(Device()))
    , TexIdx(texIdx.PrepareForInput(Device()))
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

  template<typename PtArrayType,
            typename HitRecord,
            typename HitId,
            typename LeafPortalType,
            typename Id2ArrayPortal,
            int HitBitIdx = 2,
            int ScatterBitIdx= 3>
  VTKM_EXEC void LeafIntersect(
                        const vtkm::Int32 &currentNode,
                        vec3 &origin,
                        vec3 &direction,
                        HitRecord &hrec,
                        HitId &hid,
                        float &tmin,
                        float &tmax,
                        vtkm::UInt8 &scattered,
                        Id2ArrayPortal SphereIds,
                        PtArrayType pts,
                        LeafPortalType leafs)
  {
    if (scattered & (1UL << ScatterBitIdx)){
      const vtkm::Id sphereCount = leafs.Get(currentNode);
      for (vtkm::Id i = 1; i <= sphereCount; ++i)
      {
        const vtkm::Id sphereIndex = leafs.Get(currentNode + i);
        if (sphereIndex < SphereIds.GetNumberOfValues())
        {
          auto pointIndex = SphereIds.Get(sphereIndex);
          vec3 pt = pts.Get(pointIndex[1]);
          vec3 rpt = pts.Get(pointIndex[2]);

          HitRecord  temp_rec;
          HitId temp_hid;
          auto h =   hit(origin, direction, temp_rec, temp_hid, tmin, tmax,
                         pt, rpt[0],MatIdx.Get(i-1),TexIdx.Get(i-1));
          if (h){
            tmax = temp_rec[static_cast<vtkm::Id>(HR::T)];
            hrec = temp_rec;
            hid = temp_hid;
          }
          scattered |= (h << HitBitIdx);
        }
      }
    }
  }
};


class SphereExecWrapper : public vtkm::cont::ExecutionObjectBase
{

public:
  vtkm::cont::ArrayHandle<vtkm::Id> &PointIds;
  vtkm::cont::ArrayHandle<vtkm::Id> &MatIdx;
  vtkm::cont::ArrayHandle<vtkm::Id> &TexIdx;

  SphereExecWrapper(vtkm::cont::ArrayHandle<vtkm::Id> &pointIds,
                    vtkm::cont::ArrayHandle<vtkm::Id> &matIdx,
                    vtkm::cont::ArrayHandle<vtkm::Id> &texIdx

                    )
    : PointIds(pointIds)
    , MatIdx(matIdx)
    , TexIdx(texIdx)

  {
  }

  template <typename Device>
  VTKM_CONT SphereLeafIntersector<Device> PrepareForExecution(Device) const
  {
    return SphereLeafIntersector<Device>(PointIds, MatIdx, TexIdx);
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
