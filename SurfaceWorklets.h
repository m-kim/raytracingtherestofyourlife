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
           typename SurfExec,
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
                  SurfExec &surf,
                  PtArrayType pts,
                  LeafPortalType leafs,
                  IdArrayPortal QuadIds
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
             pts,
             leafs,
             QuadIds);

  }
};
#endif
