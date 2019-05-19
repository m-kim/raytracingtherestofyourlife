#ifndef SURFACEWORKLETS_H
#define SURFACEWORKLETS_H
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/VectorAnalysis.h>
#include "Surface.h"
VTKM_EXEC


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

  using ControlSignature = void(FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  FieldInOut<>,
  ExecObject surf,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>
  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12);

  template<typename PtArrayType,
            typename HitRecord,
            typename HitId,
            typename SphereExec,
            typename LeafPortalType,
            typename Id2ArrayPortal,
            int HitBitIdx = 2,
            int ScatterBitIdx= 3>
  VTKM_EXEC
  void operator()(vtkm::Id idx,
                  const vtkm::Int32 &currentNode,
                  vec3 &origin,
                  vec3 &direction,
                  HitRecord &hrec,
                  HitId &hid,
                  float &tmin,
                  float &tmax,
                  vtkm::UInt8 &scattered,
                  SphereExec surf,
                  Id2ArrayPortal SphereIds,
                  PtArrayType pts,
                  LeafPortalType leafs
                  ) const
  {
    surf.LeafIntersect(currentNode,
                       origin,
                       direction,
                       hrec,
                       hid,
                       tmin,
                       tmax,
                       scattered,
                       SphereIds,
                       pts,
                       leafs);
  }


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


class QuadIntersect : public vtkm::worklet::WorkletMapField
{
public:

  QuadIntersect(){}
  template <typename vec3, typename Precision>
  VTKM_EXEC bool quad(const vec3& ray_origin,
                      const vec3& ray_direction,
                      const vec3& v00,
                      const vec3& v10,
                      const vec3& v11,
                      const vec3& v01,
                      Precision& u,
                      Precision& v,
                      Precision& t) const
  {

    /* An Eﬃcient Ray-Quadrilateral Intersection Test
         Ares Lagae Philip Dutr´e
         http://graphics.cs.kuleuven.be/publications/LD05ERQIT/index.html

      v01 *------------ * v11
          |\           |
          |  \         |
          |    \       |
          |      \     |
          |        \   |
          |          \ |
      v00 *------------* v10
      */
    // Rejects rays that are parallel to Q, and rays that intersect the plane of
    // Q either on the left of the line V00V01 or on the right of the line V00V10.

    vec3 E03 = v01 - v00;
    vec3 P = vtkm::Cross(ray_direction, E03);
    vec3 E01 = v10 - v00;
    Precision det = vtkm::dot(E01, P);

    if (vtkm::Abs(det) < vtkm::Epsilon<Precision>())
      return false;
    Precision inv_det = 1.0f / det;
    vec3 T = ray_origin - v00;
    Precision alpha = vtkm::dot(T, P) * inv_det;
    if (alpha < 0.0)
      return false;
    vec3 Q = vtkm::Cross(T, E01);
    Precision beta = vtkm::dot(ray_direction, Q) * inv_det;
    if (beta < 0.0)
      return false;

    if ((alpha + beta) > 1.0f)
    {

      // Rejects rays that intersect the plane of Q either on the
      // left of the line V11V10 or on the right of the line V11V01.

      vec3 E23 = v01 - v11;
      vec3 E21 = v10 - v11;
      vec3 P_prime = vtkm::Cross(ray_direction, E21);
      Precision det_prime = vtkm::dot(E23, P_prime);
      if (vtkm::Abs(det_prime) < vtkm::Epsilon<Precision>())
        return false;
      Precision inv_det_prime = 1.0f / det_prime;
      vec3 T_prime = ray_origin - v11;
      Precision alpha_prime = vtkm::dot(T_prime, P_prime) * inv_det_prime;
      if (alpha_prime < 0.0f)
        return false;
      vec3 Q_prime = vtkm::Cross(T_prime, E23);
      Precision beta_prime = vtkm::dot(ray_direction, Q_prime) * inv_det_prime;
      if (beta_prime < 0.0f)
        return false;
    }

    // Compute the ray parameter of the intersection point, and
    // reject the ray if it does not hit Q.

    t = vtkm::dot(E03, Q) * inv_det;
    if (t < 0.0)
      return false;


    // Compute the barycentric coordinates of V11
    Precision alpha_11, beta_11;
    vec3 E02 = v11 - v00;
    vec3 n = vtkm::Cross(E01, E02);

    if ((vtkm::Abs(n[0]) >= vtkm::Abs(n[1])) && (vtkm::Abs(n[0]) >= vtkm::Abs(n[2])))
    {

      alpha_11 = ((E02[1] * E03[2]) - (E02[2] * E03[1])) / n[0];
      beta_11 = ((E01[1] * E02[2]) - (E01[2] * E02[1])) / n[0];
    }
    else if ((vtkm::Abs(n[1]) >= vtkm::Abs(n[0])) && (vtkm::Abs(n[1]) >= vtkm::Abs(n[2])))
    {

      alpha_11 = ((E02[2] * E03[0]) - (E02[0] * E03[2])) / n[1];
      beta_11 = ((E01[2] * E02[0]) - (E01[0] * E02[2])) / n[1];
    }
    else
    {

      alpha_11 = ((E02[0] * E03[1]) - (E02[1] * E03[0])) / n[2];
      beta_11 = ((E01[0] * E02[1]) - (E01[1] * E02[0])) / n[2];
    }

    // Compute the bilinear coordinates of the intersection point.
    if (vtkm::Abs(alpha_11 - 1.0f) < vtkm::Epsilon<Precision>())
    {

      u = alpha;
      if (vtkm::Abs(beta_11 - 1.0f) < vtkm::Epsilon<Precision>())
        v = beta;
      else
        v = beta / ((u * (beta_11 - 1.0f)) + 1.0f);
    }
    else if (vtkm::Abs(beta_11 - 1.0) < vtkm::Epsilon<Precision>())
    {

      v = beta;
      u = alpha / ((v * (alpha_11 - 1.0f)) + 1.0f);
    }
    else
    {

      Precision A = 1.0f - beta_11;
      Precision B = (alpha * (beta_11 - 1.0f)) - (beta * (alpha_11 - 1.0f)) - 1.0f;
      Precision C = alpha;
      Precision D = (B * B) - (4.0f * A * C);
      Precision QQ = -0.5f * (B + ((B < 0.0f ? -1.0f : 1.0f) * vtkm::Sqrt(D)));
      u = QQ / A;
      if ((u < 0.0f) || (u > 1.0f))
        u = C / QQ;
      v = beta / ((u * (beta_11 - 1.0f)) + 1.0f);
    }

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
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>,
  WholeArrayInOut<>
  );
  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13);

  template<typename PtArrayType,
           typename IndexType,
           typename HitRecord,
           typename HitId,
           typename LeafPortalType,
           typename IdArrayPortal,
           int HitBitIdx = 2,
           int ScatterBitIdx= 3>
  VTKM_EXEC
  void operator()(vtkm::Id idx,
                  const vtkm::Int32 &currentNode,
                  vec3 &origin,
                  vec3 &direction,
                  HitRecord &hrec,
                  HitId &hid,
                  float &tmin,
                  float &tmax,
                  vtkm::UInt8 &scattered,
                  PtArrayType pts,
                  LeafPortalType leafs,
                  IndexType matIdx,
                  IndexType texIdx,
                  IdArrayPortal QuadIds
                  ) const

  {

    if (scattered & (1UL << ScatterBitIdx)){
      const vtkm::Id quadCount = leafs.Get(currentNode);
      for (vtkm::Id i = 1; i <= quadCount; ++i)
      {
        const vtkm::Id quadIndex = leafs.Get(currentNode + i);
        if (quadIndex < QuadIds.GetNumberOfValues())
        {
          auto pointIndex = QuadIds.Get(quadIndex);

          vec3 q, r, s, t;
          q = pts.Get(pointIndex[1]);
          r = pts.Get(pointIndex[2]);
          s = pts.Get(pointIndex[3]);
          t = pts.Get(pointIndex[4]);

          HitRecord temp_rec;

          auto h = quad(origin, direction, q,r,s,t,
                        temp_rec[static_cast<vtkm::Id>(HR::U)],
                        temp_rec[static_cast<vtkm::Id>(HR::V)],
                        temp_rec[static_cast<vtkm::Id>(HR::T)]);
          h = h && (temp_rec[static_cast<vtkm::Id>(HR::T)] < tmax) &&
              (temp_rec[static_cast<vtkm::Id>(HR::T)] > tmin);
          if (h){
            tmax = temp_rec[static_cast<vtkm::Id>(HR::T)];
            hrec = temp_rec;

            vec3 normal = vtkm::TriangleNormal(q,r,s);
            vtkm::Normalize(normal);
            if (vtkm::dot(normal, direction) > 0.f)
              normal = -normal;

            vec3 p(origin + direction * tmax);
            hrec[static_cast<vtkm::Id>(HR::Px)] = p[0];
            hrec[static_cast<vtkm::Id>(HR::Py)] = p[1];
            hrec[static_cast<vtkm::Id>(HR::Pz)] = p[2];
            hrec[static_cast<vtkm::Id>(HR::Nx)] = normal[0];
            hrec[static_cast<vtkm::Id>(HR::Ny)] = normal[1];
            hrec[static_cast<vtkm::Id>(HR::Nz)] = normal[2];
            hid[static_cast<vtkm::Id>(HI::M)] = matIdx.Get(i-1);
            hid[static_cast<vtkm::Id>(HI::T)] = texIdx.Get(i-1);
          }
          scattered |= (h << HitBitIdx);
        }
      }
    }
  }
};
#endif
