#include "CornellBox.h"

#include <vtkm/Transform3D.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include "vtkm/cont/DataSetBuilderExplicit.h"

void CornellBox::invert(vtkm::Vec<vec3,4> &pts)
{
  using vec4 = vtkm::Vec<vtkm::Float32, 4>;
  auto vec3ToVec4 = [&](vec3 in) mutable ->vec4{vec4 ret; ret[0] = in[0]; ret[1] = in[1]; ret[2] = in[2]; ret[3] = 1.f; return ret;};
  auto vec4ToVec3 = [&](vec4 in) mutable ->vec3{vec3 ret; ret[0] = in[0]; ret[1] = in[1]; ret[2] = in[2]; return ret;};

  vec3 offset(265,0,295);
  vtkm::Float32 angle = -15;
  auto translationMatrix = vtkm::Transform3DTranslate(offset[0], offset[1], offset[2]);
  auto rotationMatrix = vtkm::Transform3DRotate(angle, 0.f,1.f,0.f);
  rotationMatrix = vtkm::MatrixTranspose(rotationMatrix);
  auto mat = MatrixMultiply(translationMatrix,rotationMatrix);


  auto pt = vec3ToVec4(pts[0]);
  pts[0] = vec4ToVec3(MatrixMultiply(mat, pt));

  pt = vec3ToVec4(pts[1]);
  pts[1] = vec4ToVec3(MatrixMultiply(mat, pt));

  pt = vec3ToVec4(pts[2]);
  pts[2] = vec4ToVec3(MatrixMultiply(mat, pt));

  pt = vec3ToVec4(pts[3]);
  pts[3] = vec4ToVec3(MatrixMultiply(mat, pt));
}
void CornellBox::build()
{
  tex.Allocate(4);
  tex.GetPortalControl().Set(0, vec3(0.65, 0.05, 0.05));
  tex.GetPortalControl().Set(1, vec3(0.73, 0.73, 0.73));
  tex.GetPortalControl().Set(2, vec3(0.12, 0.45, 0.15));
  tex.GetPortalControl().Set(3, vec3(15, 15, 15));

  matType.Allocate(5);
  matType.GetPortalControl().Set(0, 0); //lambertian
  matType.GetPortalControl().Set(1, 0); //lambertian
  matType.GetPortalControl().Set(2, 0); //lambertian
  matType.GetPortalControl().Set(3, 1); //light
  matType.GetPortalControl().Set(4, 2); //dielectric

  texType.Allocate(5);
  texType.GetPortalControl().Set(0, 0); //red
  texType.GetPortalControl().Set(1, 1); //white
  texType.GetPortalControl().Set(2, 2); //green
  texType.GetPortalControl().Set(3, 3); //super bright
  texType.GetPortalControl().Set(4, 0); //dielectric


  vtkm::cont::ArrayHandle<vec3> pts1;
  pts1.Allocate(12 * 4 + 2);
  ptsIdx[0].Allocate(12*4);

  numindices[0].Allocate(12);
  shapes[0].Allocate(12);
  conn[0].Allocate(12);
  QuadIds.Allocate(12);
  vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleCounting<vtkm::Id>(0,1, 12*4), ptsIdx[0]);
  vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleCounting<vtkm::Id>(0,1, 12*4), conn[0]);

  numindices[1].Allocate(1);
  shapes[1].Allocate(1);
  conn[1].Allocate(1);
  ptsIdx[1].Allocate(1);
  SphereIds.Allocate(1);
  vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleCounting<vtkm::Id>(12*4,1, 1), ptsIdx[1]);
  vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleCounting<vtkm::Id>(0,1, 1), conn[1]);

  matIdx[0].Allocate(12);
  texIdx[0].Allocate(12);
  matIdx[1].Allocate(1);
  texIdx[1].Allocate(1);
  int cell_cnt = 0;
  int pt_idx = 0;
  auto close = [&](){  pt_idx += 4; cell_cnt++; };
  //yz_rect
  QuadIds.GetPortalControl().Set(cell_cnt, vtkm::Vec<vtkm::Id, 5>(0,pt_idx,pt_idx+1,pt_idx+2,pt_idx+3));
  matIdx[0].GetPortalControl().Set(cell_cnt, 2);
  texIdx[0].GetPortalControl().Set(cell_cnt, 2);
  shapes[0].GetPortalControl().Set(cell_cnt, vtkm::CELL_SHAPE_QUAD);
  numindices[0].GetPortalControl().Set(cell_cnt, 4);
  pts1.GetPortalControl().Set(pt_idx, vec3(555,0,0));
  pts1.GetPortalControl().Set(pt_idx+1, vec3(555,555,0));
  pts1.GetPortalControl().Set(pt_idx+2, vec3(555,555,555));
  pts1.GetPortalControl().Set(pt_idx+3, vec3(555,0,555));
  close();

  QuadIds.GetPortalControl().Set(cell_cnt, vtkm::Vec<vtkm::Id, 5>(0,pt_idx,pt_idx+1,pt_idx+2,pt_idx+3));
  matIdx[0].GetPortalControl().Set(cell_cnt, 0);
  texIdx[0].GetPortalControl().Set(cell_cnt, 0);
  shapes[0].GetPortalControl().Set(cell_cnt, vtkm::CELL_SHAPE_QUAD);
  numindices[0].GetPortalControl().Set(cell_cnt, 4);
  pts1.GetPortalControl().Set(pt_idx, vec3(0,0,0));
  pts1.GetPortalControl().Set(pt_idx+1, vec3(0,555,0));
  pts1.GetPortalControl().Set(pt_idx+2, vec3(0,555,555));
  pts1.GetPortalControl().Set(pt_idx+3, vec3(0,0,555));
  close();

  //xz_rect
  QuadIds.GetPortalControl().Set(cell_cnt, vtkm::Vec<vtkm::Id, 5>(0,pt_idx,pt_idx+1,pt_idx+2,pt_idx+3));
  matIdx[0].GetPortalControl().Set(cell_cnt, 3);
  texIdx[0].GetPortalControl().Set(cell_cnt, 3);
  shapes[0].GetPortalControl().Set(cell_cnt, vtkm::CELL_SHAPE_QUAD);
  numindices[0].GetPortalControl().Set(cell_cnt, 4);
  pts1.GetPortalControl().Set(pt_idx, vec3(213,554,227));
  pts1.GetPortalControl().Set(pt_idx+1, vec3(343,554,227));
  pts1.GetPortalControl().Set(pt_idx+2, vec3(343,554,332));
  pts1.GetPortalControl().Set(pt_idx+3, vec3(213,554,332));
  close();

  QuadIds.GetPortalControl().Set(cell_cnt, vtkm::Vec<vtkm::Id, 5>(0,pt_idx,pt_idx+1,pt_idx+2,pt_idx+3));
  matIdx[0].GetPortalControl().Set(cell_cnt, 1);
  texIdx[0].GetPortalControl().Set(cell_cnt, 1);
  shapes[0].GetPortalControl().Set(cell_cnt, vtkm::CELL_SHAPE_QUAD);
  numindices[0].GetPortalControl().Set(cell_cnt, 4);
  pts1.GetPortalControl().Set(pt_idx, vec3(0,555,0));
  pts1.GetPortalControl().Set(pt_idx+1, vec3(555,555,0));
  pts1.GetPortalControl().Set(pt_idx+2, vec3(555,555,555));
  pts1.GetPortalControl().Set(pt_idx+3, vec3(0,555,555));
  close();

  QuadIds.GetPortalControl().Set(cell_cnt, vtkm::Vec<vtkm::Id, 5>(0,pt_idx,pt_idx+1,pt_idx+2,pt_idx+3));
  matIdx[0].GetPortalControl().Set(cell_cnt, 1);
  texIdx[0].GetPortalControl().Set(cell_cnt, 1);
  shapes[0].GetPortalControl().Set(cell_cnt, vtkm::CELL_SHAPE_QUAD);
  numindices[0].GetPortalControl().Set(cell_cnt, 4);
  pts1.GetPortalControl().Set(pt_idx, vec3(0,0,0));
  pts1.GetPortalControl().Set(pt_idx+1, vec3(555,0,0));
  pts1.GetPortalControl().Set(pt_idx+2, vec3(555,0,555));
  pts1.GetPortalControl().Set(pt_idx+3, vec3(0,0,555));
  close();

  //xy_rect
  QuadIds.GetPortalControl().Set(cell_cnt, vtkm::Vec<vtkm::Id, 5>(0,pt_idx,pt_idx+1,pt_idx+2,pt_idx+3));
  matIdx[0].GetPortalControl().Set(cell_cnt, 1);
  texIdx[0].GetPortalControl().Set(cell_cnt, 1);
  shapes[0].GetPortalControl().Set(cell_cnt, vtkm::CELL_SHAPE_QUAD);
  numindices[0].GetPortalControl().Set(cell_cnt, 4);
  pts1.GetPortalControl().Set(pt_idx, vec3(0,0,555));
  pts1.GetPortalControl().Set(pt_idx+1, vec3(555,0,555));
  pts1.GetPortalControl().Set(pt_idx+2, vec3(555,555,555));
  pts1.GetPortalControl().Set(pt_idx+3, vec3(0,555,555));
  close();


//  //small box
  //xy
  QuadIds.GetPortalControl().Set(cell_cnt, vtkm::Vec<vtkm::Id, 5>(0,pt_idx,pt_idx+1,pt_idx+2,pt_idx+3));
  vtkm::Vec<vec3,4> pts;
  pts[0] = vec3(0,0,165);
  pts[1] = vec3(165,0,165);
  pts[2] = vec3(165,330,165);
  pts[3] = vec3(0,330,165);
  invert(pts);
  matIdx[0].GetPortalControl().Set(cell_cnt, 1);
  texIdx[0].GetPortalControl().Set(cell_cnt, 1);
  shapes[0].GetPortalControl().Set(cell_cnt, vtkm::CELL_SHAPE_QUAD);
  numindices[0].GetPortalControl().Set(cell_cnt, 4);
  pts1.GetPortalControl().Set(pt_idx, pts[0]);
  pts1.GetPortalControl().Set(pt_idx+1, pts[1]);
  pts1.GetPortalControl().Set(pt_idx+2, pts[2]);
  pts1.GetPortalControl().Set(pt_idx+3, pts[3]);
  close();

  QuadIds.GetPortalControl().Set(cell_cnt, vtkm::Vec<vtkm::Id, 5>(0,pt_idx,pt_idx+1,pt_idx+2,pt_idx+3));
  pts[0] = vec3(0,0,0);
  pts[1] = vec3(165,0,0);
  pts[2] = vec3(165,330,0);
  pts[3] = vec3(0,330,0);
  invert(pts);
  matIdx[0].GetPortalControl().Set(cell_cnt, 1);
  texIdx[0].GetPortalControl().Set(cell_cnt, 1);
  shapes[0].GetPortalControl().Set(cell_cnt, vtkm::CELL_SHAPE_QUAD);
  numindices[0].GetPortalControl().Set(cell_cnt, 4);
  pts1.GetPortalControl().Set(pt_idx, pts[0]);
  pts1.GetPortalControl().Set(pt_idx+1, pts[1]);
  pts1.GetPortalControl().Set(pt_idx+2, pts[2]);
  pts1.GetPortalControl().Set(pt_idx+3, pts[3]);
  close();

  //yz
  QuadIds.GetPortalControl().Set(cell_cnt, vtkm::Vec<vtkm::Id, 5>(0,pt_idx,pt_idx+1,pt_idx+2,pt_idx+3));
  pts[0] = vec3(165,0,0);
  pts[1] = vec3(165,330,0);
  pts[2] = vec3(165,330,165);
  pts[3] = vec3(165, 0, 165);
  invert(pts);
  matIdx[0].GetPortalControl().Set(cell_cnt, 1);
  texIdx[0].GetPortalControl().Set(cell_cnt, 1);
  shapes[0].GetPortalControl().Set(cell_cnt, vtkm::CELL_SHAPE_QUAD);
  numindices[0].GetPortalControl().Set(cell_cnt, 4);
  pts1.GetPortalControl().Set(pt_idx, pts[0]);
  pts1.GetPortalControl().Set(pt_idx+1, pts[1]);
  pts1.GetPortalControl().Set(pt_idx+2, pts[2]);
  pts1.GetPortalControl().Set(pt_idx+3, pts[3]);
  close();

  QuadIds.GetPortalControl().Set(cell_cnt, vtkm::Vec<vtkm::Id, 5>(0,pt_idx,pt_idx+1,pt_idx+2,pt_idx+3));
  pts[0] = vec3(0,0,0);
  pts[1] = vec3(0,330,0);
  pts[2] = vec3(0,330,165);
  pts[3] = vec3(0, 0, 165);
  invert(pts);
  matIdx[0].GetPortalControl().Set(cell_cnt, 1);
  texIdx[0].GetPortalControl().Set(cell_cnt, 1);
  shapes[0].GetPortalControl().Set(cell_cnt, vtkm::CELL_SHAPE_QUAD);
  numindices[0].GetPortalControl().Set(cell_cnt, 4);
  pts1.GetPortalControl().Set(pt_idx, pts[0]);
  pts1.GetPortalControl().Set(pt_idx+1, pts[1]);
  pts1.GetPortalControl().Set(pt_idx+2, pts[2]);
  pts1.GetPortalControl().Set(pt_idx+3, pts[3]);
  close();


  //xz_rect
  QuadIds.GetPortalControl().Set(cell_cnt, vtkm::Vec<vtkm::Id, 5>(0,pt_idx,pt_idx+1,pt_idx+2,pt_idx+3));
  pts[0] = vec3(0,333,0);
  pts[1] = vec3(165,330,0);
  pts[2] = vec3(165,330,165);
  pts[3] = vec3(0, 330, 165);
  invert(pts);
  matIdx[0].GetPortalControl().Set(cell_cnt, 1);
  texIdx[0].GetPortalControl().Set(cell_cnt, 1);
  shapes[0].GetPortalControl().Set(cell_cnt, vtkm::CELL_SHAPE_QUAD);
  numindices[0].GetPortalControl().Set(cell_cnt, 4);
  pts1.GetPortalControl().Set(pt_idx, pts[0]);
  pts1.GetPortalControl().Set(pt_idx+1, pts[1]);
  pts1.GetPortalControl().Set(pt_idx+2, pts[2]);
  pts1.GetPortalControl().Set(pt_idx+3, pts[3]);
  close();

  QuadIds.GetPortalControl().Set(cell_cnt, vtkm::Vec<vtkm::Id, 5>(0,pt_idx,pt_idx+1,pt_idx+2,pt_idx+3));
  pts[0] = vec3(0,0,0);
  pts[1] = vec3(165,0,0);
  pts[2] = vec3(165,0,165);
  pts[3] = vec3(0, 0, 165);
  invert(pts);
  matIdx[0].GetPortalControl().Set(cell_cnt, 1);
  texIdx[0].GetPortalControl().Set(cell_cnt, 1);
  shapes[0].GetPortalControl().Set(cell_cnt, vtkm::CELL_SHAPE_QUAD);
  numindices[0].GetPortalControl().Set(cell_cnt, 4);
  pts1.GetPortalControl().Set(pt_idx, pts[0]);
  pts1.GetPortalControl().Set(pt_idx+1, pts[1]);
  pts1.GetPortalControl().Set(pt_idx+2, pts[2]);
  pts1.GetPortalControl().Set(pt_idx+3, pts[3]);
  close();

  //    //sphere
    SphereIds.GetPortalControl().Set(0, pt_idx);

    matIdx[1].GetPortalControl().Set(0, 4);
    texIdx[1].GetPortalControl().Set(0, 0);
    shapes[1].GetPortalControl().Set(0, vtkm::CELL_SHAPE_VERTEX);
    numindices[1].GetPortalControl().Set(0, 2);
    pts1.GetPortalControl().Set(pt_idx, vec3(190,90,190));

    SphereRadii.Allocate(1);
    SphereRadii.GetPortalControl().Set(0, 90);
    coord.SetData(pts1);
}

vtkm::cont::DataSet CornellBox::buildDataSet()
{
  if (coord.GetData().GetNumberOfValues() == 0){
    build();

  }

  vtkm::cont::DataSetBuilderExplicit dsb;

  vtkm::cont::DataSet ds;
  ds.AddCoordinateSystem(coord);
  vtkm::cont::CellSetExplicit<> cellsetQuads;
  cellsetQuads.Fill(12*4, shapes[0], numindices[0], conn[0], ptsIdx[0]);
  ds.AddCellSet(cellsetQuads);


  vtkm::cont::CellSetExplicit<> cellsetSphere;
  cellsetSphere.Fill(1,
                     shapes[1],
                    numindices[1],
                    conn[1],
                    ptsIdx[1]);

  ds.AddCellSet(cellsetSphere);


  return ds;
}
