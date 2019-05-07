#ifndef WORKLETS_H
#define WORKLETS_H

#include <vtkm/worklet/AutoDispatcherMapField.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/ArrayHandleCounting.h>

class RayGen : public vtkm::worklet::WorkletMapField
{
//  vtkm::Int32 w;
//  vtkm::Int32 h;
//  vtkm::Int32 Minx;
//  vtkm::Int32 Miny;
//  vtkm::Int32 SubsetWidth;
//  vtkm::Vec<vtkm::Float32, 3> nlook; // normalized look
//  vtkm::Vec<vtkm::Float32, 3> delta_x;
//  vtkm::Vec<vtkm::Float32, 3> delta_y;
//  vtkm::UInt32 RayCount;
  VTKM_CONT

//  RayGen(vtkm::Int32 width,
//                    vtkm::Int32 height,
//                    vtkm::Float32 fovX,
//                    vtkm::Float32 fovY,
//                    vtkm::Vec<vtkm::Float32, 3> look,
//                    vtkm::Vec<vtkm::Float32, 3> up,
//                    vtkm::Float32 _zoom,
//                    vtkm::Int32 subsetWidth,
//                    vtkm::Int32 minx,
//                    vtkm::Int32 miny,
//                    vtkm::UInt32 rayCount)
//    : w(width)
//    , h(height)
//    , Minx(minx)
//    , Miny(miny)
//    , SubsetWidth(subsetWidth)
//    , RayCount(rayCount)
  RayGen()
  {

  }

//  using ControlSignature = void(FieldOut<>, FieldOut<>, FieldOut<>, FieldOut<>);

//  using ExecutionSignature = void(WorkIndex, _1, _2, _3, _4);
//  template <typename Precision>
//  VTKM_EXEC void operator()(vtkm::Id idx,
//                            Precision& rayDirX,
//                            Precision& rayDirY,
//                            Precision& rayDirZ,
//                            vtkm::Id& pixelIndex) const
//  {
//  }

//  VTKM_CONT void run()
//  {
//    vtkm::worklet::AutoDispatcherMapField<RayGen>(
//          RayGen()).Invoke();

//  }
};

#endif
