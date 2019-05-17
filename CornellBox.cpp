#include "CornellBox.h"
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
  norms.Allocate(16);

  //yz_rect
  translateOffset[0].Allocate(4);
  angleArray[0].Allocate(4);
    flipped[0].Allocate(4);
    matIdx[0].Allocate(4);
    texIdx[0].Allocate(4);
    pts1[0].Allocate(4);
    pts2[0].Allocate(4);
    cellTypeArray[0] = 0;
    matIdx[0].GetPortalControl().Set(0, 2);
    texIdx[0].GetPortalControl().Set(0, 2);
    pts1[0].GetPortalControl().Set(0, vec3(555,0,0));
    pts2[0].GetPortalControl().Set(0, vec3(555,555,555));
    flipped[0].GetPortalControl().Set(0, 1);
    translateOffset[0].GetPortalControl().Set(0, vec3(0.0f));
    angleArray[0].GetPortalControl().Set(0, 0);


    matIdx[0].GetPortalControl().Set(1, 0);
    texIdx[0].GetPortalControl().Set(1, 0);
    pts1[0].GetPortalControl().Set(1, vec3(0,0,0));
    pts2[0].GetPortalControl().Set(1, vec3(0,555,555));
    flipped[0].GetPortalControl().Set(1, 0);
    translateOffset[0].GetPortalControl().Set(1, vec3(0.0f));
    angleArray[0].GetPortalControl().Set(1, 0);

    //xz_rect
    translateOffset[1].Allocate(5);
    angleArray[1].Allocate(5);
    flipped[1].Allocate(5);
    matIdx[1].Allocate(5);
    texIdx[1].Allocate(5);
    pts1[1].Allocate(5);
    pts2[1].Allocate(5);
    cellTypeArray[1] =  1;
    matIdx[1].GetPortalControl().Set(0, 3);
    texIdx[1].GetPortalControl().Set(0, 3);
    pts1[1].GetPortalControl().Set(0, vec3(213,554,227));
    pts2[1].GetPortalControl().Set(0, vec3(343,554,332));
    flipped[1].GetPortalControl().Set(0, 1);

    matIdx[1].GetPortalControl().Set(1, 1);
    texIdx[1].GetPortalControl().Set(1, 1);
    pts1[1].GetPortalControl().Set(1, vec3(0,555,0));
    pts2[1].GetPortalControl().Set(1, vec3(555,555,555));
    flipped[1].GetPortalControl().Set(1, 1);

    matIdx[1].GetPortalControl().Set(2, 1);
    texIdx[1].GetPortalControl().Set(2, 1);
    pts1[1].GetPortalControl().Set(2, vec3(0,0,0));
    pts2[1].GetPortalControl().Set(2, vec3(555,0,555));
    flipped[1].GetPortalControl().Set(2, 0);

    //xy_rect
    translateOffset[2].Allocate(3);
    angleArray[2].Allocate(3);
    flipped[2].Allocate(3);
    matIdx[2].Allocate(3);
    texIdx[2].Allocate(3);
    pts1[2].Allocate(3);
    pts2[2].Allocate(3);
    cellTypeArray[2] =  2;
    matIdx[2].GetPortalControl().Set(0, 1);
    texIdx[2].GetPortalControl().Set(0, 1);
    pts1[2].GetPortalControl().Set(0, vec3(0,0,555));
    pts2[2].GetPortalControl().Set(0, vec3(555,555,555));
    translateOffset[2].GetPortalControl().Set(0, vec3(0,0,0));
    angleArray[2].GetPortalControl().Set(0,0);
    flipped[2].GetPortalControl().Set(0, 1);

    //sphere
    flipped[3].Allocate(1);
    matIdx[3].Allocate(1);
    texIdx[3].Allocate(1);
    pts1[3].Allocate(1);
    pts2[3].Allocate(1);
    cellTypeArray[3] =  3;
    matIdx[3].GetPortalControl().Set(0, 4);
    texIdx[3].GetPortalControl().Set(0, 0);
    pts1[3].GetPortalControl().Set(0, vec3(190,90,190));
    pts2[3].GetPortalControl().Set(0, vec3(90,0,0));
    flipped[3].GetPortalControl().Set(0, 0);

    flipped[4].Allocate(1);
    matIdx[4].Allocate(1);
    texIdx[4].Allocate(1);
    pts1[4].Allocate(1);
    pts2[4].Allocate(1);
    cellTypeArray[4] =  4;
    matIdx[4].GetPortalControl().Set(0, 1);
    texIdx[4].GetPortalControl().Set(0, 1);
    flipped[4].GetPortalControl().Set(0, 0);

  //small box
  const vec3 constp1(0,0,0);
  const vec3 constp2(165,330, 165);

  vec3 p1 = constp1;
  vec3 p2 = constp2;

  //xy
  p1 = vec3(0,0,165);
  p2 = vec3(165,330,165);
  vec3 offset(265,0,295);
  int ct = 2;
  matIdx[ct].GetPortalControl().Set(1, 1);
  texIdx[ct].GetPortalControl().Set(1, 1);
  pts1[ct].GetPortalControl().Set(1, p1);
  pts2[ct].GetPortalControl().Set(1, p2);
  translateOffset[ct].GetPortalControl().Set(1, offset);
  angleArray[ct].GetPortalControl().Set(1, 15);
  flipped[ct].GetPortalControl().Set(1, 0);


  matIdx[ct].GetPortalControl().Set(2, 1);
  texIdx[ct].GetPortalControl().Set(2, 1);
  pts1[ct].GetPortalControl().Set(2, vec3(0,0,0));
  pts2[ct].GetPortalControl().Set(2,  vec3(165,330,165));
  translateOffset[ct].GetPortalControl().Set(2, offset);
  angleArray[ct].GetPortalControl().Set(2, 15);
  flipped[ct].GetPortalControl().Set(2, 1);

  //yz
  ct = 0;
  p1 = vec3(165,0,0);
  p2 = vec3(165,330,165);
  matIdx[ct].GetPortalControl().Set(2, 1);
  texIdx[ct].GetPortalControl().Set(2, 1);
  pts1[ct].GetPortalControl().Set(2, p1);
  pts2[ct].GetPortalControl().Set(2, p2);
  translateOffset[ct].GetPortalControl().Set(2, offset);
  angleArray[ct].GetPortalControl().Set(2, 15);
  flipped[ct].GetPortalControl().Set(2, 0);

  p1 = vec3(0,0,0);
  p2 = vec3(165,330,165);
  matIdx[ct].GetPortalControl().Set(3, 1);
  texIdx[ct].GetPortalControl().Set(3, 1);
  pts1[ct].GetPortalControl().Set(3, p1);
  pts2[ct].GetPortalControl().Set(3, p2);
  translateOffset[ct].GetPortalControl().Set(3, offset);
  angleArray[ct].GetPortalControl().Set(3, 15);
  flipped[ct].GetPortalControl().Set(3, 1);

  //xz_rect
  ct = 1;
  p1 = vec3(0,330,0);
  p2 = vec3(165,330,165);
  matIdx[ct].GetPortalControl().Set(3, 1);
  texIdx[ct].GetPortalControl().Set(3, 1);
  pts1[ct].GetPortalControl().Set(3, p1);
  pts2[ct].GetPortalControl().Set(3, p2);
  translateOffset[ct].GetPortalControl().Set(3, offset);
  angleArray[ct].GetPortalControl().Set(3, 15);
  flipped[ct].GetPortalControl().Set(3, 0);

  p1 = vec3(0,0,0);
  p2 = vec3(165,330,165);
  matIdx[ct].GetPortalControl().Set(4, 1);
  texIdx[ct].GetPortalControl().Set(4, 1);
  pts1[ct].GetPortalControl().Set(4, p1);
  pts2[ct].GetPortalControl().Set(4, p2);
  translateOffset[ct].GetPortalControl().Set(4, offset);
  angleArray[ct].GetPortalControl().Set(4, 15);
  flipped[ct].GetPortalControl().Set(4, 1);

}
