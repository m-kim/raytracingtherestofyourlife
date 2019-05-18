#include "CornellBox.h"

#include <vtkm/Transform3D.h>

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


  cellTypeArray.resize(5);

  //yz_rect
    matIdx[0].Allocate(4);
    texIdx[0].Allocate(4);
    pts1[0].Allocate(16);
    matIdx[0].GetPortalControl().Set(0, 2);
    texIdx[0].GetPortalControl().Set(0, 2);
    pts1[0].GetPortalControl().Set(0, vec3(555,0,0));
    pts1[0].GetPortalControl().Set(1, vec3(555,555,0));
    pts1[0].GetPortalControl().Set(2, vec3(555,555,555));
    pts1[0].GetPortalControl().Set(3, vec3(555,0,555));


    matIdx[0].GetPortalControl().Set(1, 0);
    texIdx[0].GetPortalControl().Set(1, 0);
    pts1[0].GetPortalControl().Set(4, vec3(0,0,0));
    pts1[0].GetPortalControl().Set(5, vec3(0,555,0));
    pts1[0].GetPortalControl().Set(6, vec3(0,555,555));
    pts1[0].GetPortalControl().Set(7, vec3(0,0,555));

    //xz_rect
    matIdx[1].Allocate(5);
    texIdx[1].Allocate(5);
    pts1[1].Allocate(20);
    cellTypeArray[1] =  1;
    matIdx[1].GetPortalControl().Set(0, 3);
    texIdx[1].GetPortalControl().Set(0, 3);
    pts1[1].GetPortalControl().Set(0, vec3(213,554,227));
    pts1[1].GetPortalControl().Set(1, vec3(343,554,227));
    pts1[1].GetPortalControl().Set(2, vec3(343,554,332));
    pts1[1].GetPortalControl().Set(3, vec3(213,554,332));

    matIdx[1].GetPortalControl().Set(1, 1);
    texIdx[1].GetPortalControl().Set(1, 1);
    pts1[1].GetPortalControl().Set(4, vec3(0,555,0));
    pts1[1].GetPortalControl().Set(5, vec3(555,555,0));
    pts1[1].GetPortalControl().Set(6, vec3(555,555,555));
    pts1[1].GetPortalControl().Set(7, vec3(0,555,555));

    matIdx[1].GetPortalControl().Set(2, 1);
    texIdx[1].GetPortalControl().Set(2, 1);
    pts1[1].GetPortalControl().Set(8, vec3(0,0,0));
    pts1[1].GetPortalControl().Set(9, vec3(555,0,0));
    pts1[1].GetPortalControl().Set(10, vec3(555,0,555));
    pts1[1].GetPortalControl().Set(11, vec3(0,0,555));

    //xy_rect
    matIdx[2].Allocate(3);
    texIdx[2].Allocate(3);
    pts1[2].Allocate(12);
    cellTypeArray[2] =  2;
    matIdx[2].GetPortalControl().Set(0, 1);
    texIdx[2].GetPortalControl().Set(0, 1);
    pts1[2].GetPortalControl().Set(0, vec3(0,0,555));
    pts1[2].GetPortalControl().Set(1, vec3(555,0,555));
    pts1[2].GetPortalControl().Set(2, vec3(555,555,555));
    pts1[2].GetPortalControl().Set(3, vec3(0,555,555));

//    //sphere
//    flipped[3].Allocate(1);
    matIdx[3].Allocate(1);
    texIdx[3].Allocate(1);
    pts1[3].Allocate(2);
//    pts2[3].Allocate(1);
    cellTypeArray[3] =  3;
    matIdx[3].GetPortalControl().Set(0, 4);
    texIdx[3].GetPortalControl().Set(0, 0);
    pts1[3].GetPortalControl().Set(0, vec3(190,90,190));
    pts1[3].GetPortalControl().Set(1, vec3(90,0,0));
//    flipped[3].GetPortalControl().Set(0, 0);

//    flipped[4].Allocate(1);
//    matIdx[4].Allocate(1);
//    texIdx[4].Allocate(1);
//    pts1[4].Allocate(1);
//    pts2[4].Allocate(1);
    cellTypeArray[4] =  4;
//    matIdx[4].GetPortalControl().Set(0, 1);
//    texIdx[4].GetPortalControl().Set(0, 1);
//    flipped[4].GetPortalControl().Set(0, 0);

//  //small box
//  const vec3 constp1(0,0,0);
//  const vec3 constp2(165,330, 165);

//  vec3 p1 = constp1;
//  vec3 p2 = constp2;

  //xy

  vtkm::Vec<vec3,4> pts;
  pts[0] = vec3(0,0,165);
  pts[1] = vec3(165,0,165);
  pts[2] = vec3(165,330,165);
  pts[3] = vec3(0,330,165);

  invert(pts);

  int ct = 2;
  matIdx[ct].GetPortalControl().Set(1, 1);
  texIdx[ct].GetPortalControl().Set(1, 1);
  pts1[ct].GetPortalControl().Set(4, pts[0]);
  pts1[ct].GetPortalControl().Set(5, pts[1]);
  pts1[ct].GetPortalControl().Set(6, pts[2]);
  pts1[ct].GetPortalControl().Set(7, pts[3]);

  //pts2[ct].GetPortalControl().Set(1, p2);
//  translateOffset[ct].GetPortalControl().Set(1, offset);
//  angleArray[ct].GetPortalControl().Set(1, 15);
//  flipped[ct].GetPortalControl().Set(1, 0);

  pts[0] = vec3(0,0,0);
  pts[1] = vec3(165,0,0);
  pts[2] = vec3(165,330,0);
  pts[3] = vec3(0,330,0);

  invert(pts);


  matIdx[ct].GetPortalControl().Set(2, 1);
  texIdx[ct].GetPortalControl().Set(2, 1);
  pts1[ct].GetPortalControl().Set(8, pts[0]);
  pts1[ct].GetPortalControl().Set(9, pts[1]);
  pts1[ct].GetPortalControl().Set(10, pts[2]);
  pts1[ct].GetPortalControl().Set(11, pts[3]);

  pts[0] = vec3(165,0,0);
  pts[1] = vec3(165,330,0);
  pts[2] = vec3(165,330,165);
  pts[3] = vec3(165, 0, 165);

  invert(pts);
  //yz
  ct = 0;
  matIdx[ct].GetPortalControl().Set(2, 1);
  texIdx[ct].GetPortalControl().Set(2, 1);
  pts1[ct].GetPortalControl().Set(8, pts[0]);
  pts1[ct].GetPortalControl().Set(9, pts[1]);
  pts1[ct].GetPortalControl().Set(10, pts[2]);
  pts1[ct].GetPortalControl().Set(11, pts[3]);

  pts[0] = vec3(0,0,0);
  pts[1] = vec3(0,330,0);
  pts[2] = vec3(0,330,165);
  pts[3] = vec3(0, 0, 165);

  invert(pts);
//  p1 = vec3(0,0,0);
//  p2 = vec3(165,330,165);
  matIdx[ct].GetPortalControl().Set(3, 1);
  texIdx[ct].GetPortalControl().Set(3, 1);
  pts1[ct].GetPortalControl().Set(12, pts[0]);
  pts1[ct].GetPortalControl().Set(13, pts[1]);
  pts1[ct].GetPortalControl().Set(14, pts[2]);
  pts1[ct].GetPortalControl().Set(15, pts[3]);


  pts[0] = vec3(0,333,0);
  pts[1] = vec3(165,330,0);
  pts[2] = vec3(165,330,165);
  pts[3] = vec3(0, 330, 165);

  invert(pts);

  //xz_rect
  ct = 1;
  matIdx[ct].GetPortalControl().Set(3, 1);
  texIdx[ct].GetPortalControl().Set(3, 1);
  pts1[ct].GetPortalControl().Set(12, pts[0]);
  pts1[ct].GetPortalControl().Set(13, pts[1]);
  pts1[ct].GetPortalControl().Set(14, pts[2]);
  pts1[ct].GetPortalControl().Set(15, pts[3]);

  pts[0] = vec3(0,0,0);
  pts[1] = vec3(165,0,0);
  pts[2] = vec3(165,0,165);
  pts[3] = vec3(0, 0, 165);

  invert(pts);

  matIdx[ct].GetPortalControl().Set(4, 1);
  texIdx[ct].GetPortalControl().Set(4, 1);
  pts1[ct].GetPortalControl().Set(16, pts[0]);
  pts1[ct].GetPortalControl().Set(17, pts[1]);
  pts1[ct].GetPortalControl().Set(18, pts[2]);
  pts1[ct].GetPortalControl().Set(19, pts[3]);

}
