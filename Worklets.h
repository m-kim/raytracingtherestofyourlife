#ifndef WORKLETS_H
#define WORKLETS_H

#include <vtkm/worklet/AutoDispatcherMapField.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/rendering/xorShift.h>

class UVGen : public vtkm::worklet::WorkletMapField
{
public:
  vtkm::UInt32 RayCount;
  VTKM_CONT

  UVGen(vtkm::Id _nx,
        vtkm::Id _ny,
        vtkm::UInt32 rayCount)
    : numx(_nx),
      numy(_ny),
      RayCount(rayCount)
  {
  }

  using ControlSignature = void(FieldOut<>);

  using ExecutionSignature = void(WorkIndex, _1);
  template <typename Precision>
  VTKM_EXEC void operator()(vtkm::Id &idx,
                            Precision& uv) const
  {
    vtkm::Vec<vtkm::UInt32, 4> randState;
    randState[0] = vtkm::random::xorshift::getRand32(RayCount*1) + 1;
    randState[1] = vtkm::random::xorshift::getRand32(RayCount*2) + 2;
    randState[2] = vtkm::random::xorshift::getRand32(RayCount*3) + 3;
    randState[3] = vtkm::random::xorshift::getRand32(RayCount*4) + 4; //arbitrary random state based off number of rays being shot through

    vtkm::random::xorshift::getRandF(randState);
    vtkm::random::xorshift::getRandF(randState);
    vtkm::random::xorshift::getRandF(randState);
    vtkm::random::xorshift::getRandF(randState);

    vtkm::Id i = idx % numx;
    vtkm::Id j = idx / numx;
    uv[0] = vtkm::Float32(i + vtkm::random::xorshift::getRandF(randState)) / numx;
    uv[1] = vtkm::Float32(j + vtkm::random::xorshift::getRandF(randState)) / numy;

  }

  const vtkm::Id numx, numy;
};

#endif
