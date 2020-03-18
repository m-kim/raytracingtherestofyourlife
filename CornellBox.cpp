#include "CornellBox.h"

#include <vtkm/Transform3D.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include "vtkm/cont/DataSetBuilderExplicit.h"
#include "pathtracing/SphereExtractor.h"
#include <vtkm/rendering/raytracing/QuadExtractor.h>

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

void CornellBox::buildQuad(vtkm::Vec<vec3,4> &pts,
    std::vector<vtkm::UInt8> &shapes,
    std::vector<vtkm::IdComponent> &numindices,
    std::vector<vec3> &pts1,
    std::vector<vtkm::Float32> &vecField,
    std::vector<vtkm::Id> &vecMatIdx,
    std::vector<vtkm::Id> &vecTexIdx,
                           int midx,
                           int tidx)

{
  vecMatIdx.push_back(midx);
  vecTexIdx.push_back(tidx);
  shapes.push_back( vtkm::CELL_SHAPE_QUAD);
  numindices.push_back( 4);
  pts1.push_back( pts[0]/555.0);
  pts1.push_back( pts[1]/555.0);
  pts1.push_back( pts[2]/555.0);
  pts1.push_back( pts[3]/555.0);
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));

}

void CornellBox::buildBox(const vec3 &near,
              const vec3 &far,
              std::vector<vtkm::UInt8> &shapes,
              std::vector<vtkm::IdComponent> &numindices,
              std::vector<vec3> &pts1,
              std::vector<vtkm::Float32> &vecField,
              std::vector<vtkm::Id> &vecMatIdx,
              std::vector<vtkm::Id> &vecTexIdx,
                                     int midx,
                                     int tidx)
{
  vtkm::Vec<vec3,4> pts;
  pts[0] = vec3(near[0], near[1], near[2]);
  pts[1] = vec3(far[0], near[1], near[2]);
  pts[2] = vec3(far[0], far[1], near[2]);
  pts[3] = vec3(near[0], far[1], near[2]);
  buildQuad(pts,
            shapes,
            numindices,
            pts1,
            vecField,
            vecMatIdx,
            vecTexIdx,
      midx, tidx);

  pts[0] = vec3(near[0], near[1], far[2]);
  pts[1] = vec3(far[0], near[1], far[2]);
  pts[2] = vec3(far[0], far[1], far[2]);
  pts[3] = vec3(near[0], far[1], far[2]);
  buildQuad(pts,
            shapes,
            numindices,
            pts1,
            vecField,
            vecMatIdx,
            vecTexIdx,
      midx,tidx);


  pts[0] = vec3(near[0], far[1], near[2]);
  pts[1] = vec3(far[0], far[1], near[2]);
  pts[2] = vec3(far[0], far[1], far[2]);
  pts[3] = vec3(near[0], far[1], far[2]);
  buildQuad(pts,
            shapes,
            numindices,
            pts1,
            vecField,
            vecMatIdx,
            vecTexIdx,
      midx,tidx);
  pts[0] = vec3(near[0], near[1], near[2]);
  pts[1] = vec3(far[0], near[1], near[2]);
  pts[2] = vec3(far[0], far[1], near[2]);
  pts[3] = vec3(near[0], far[1], near[2]);
  buildQuad(pts,
            shapes,
            numindices,
            pts1,
            vecField,
            vecMatIdx,
            vecTexIdx,
      midx,tidx);
  pts[0] = vec3(near[0], near[1], far[2]);
  pts[1] = vec3(far[0], near[1], far[2]);
  pts[2] = vec3(far[0], far[1], far[2]);
  pts[3] = vec3(near[0], far[1], far[2]);
  buildQuad(pts,
            shapes,
            numindices,
            pts1,
            vecField,
            vecMatIdx,
            vecTexIdx,
      midx,tidx);

}

vtkm::cont::DataSet CornellBox::buildDataSet()
{
  tex.Allocate(4);
  tex.WritePortal().Set(0, vec3(0.65, 0.05, 0.05));
  tex.WritePortal().Set(1, vec3(0.73, 0.73, 0.73));
  tex.WritePortal().Set(2, vec3(0.12, 0.45, 0.15));
  tex.WritePortal().Set(3, vec3(15, 15, 15));

  matType.Allocate(5);
  matType.WritePortal().Set(0, 0); //lambertian
  matType.WritePortal().Set(1, 0); //lambertian
  matType.WritePortal().Set(2, 0); //lambertian
  matType.WritePortal().Set(3, 1); //light
  matType.WritePortal().Set(4, 2); //dielectric

  texType.Allocate(5);
  texType.WritePortal().Set(0, 0); //red
  texType.WritePortal().Set(1, 1); //white
  texType.WritePortal().Set(2, 2); //green
  texType.WritePortal().Set(3, 3); //super bright
  texType.WritePortal().Set(4, 0); //dielectric


  std::vector<vtkm::UInt8> shapes;
  std::vector<vtkm::IdComponent> numindices;
  vtkm::cont::ArrayHandle<vtkm::Id> conn;
  std::vector<vec3> pts1;

  std::vector<vtkm::Float32> vecField;
  std::vector<vtkm::Id> vecMatIdx[2], vecTexIdx[2];



  //yz_rect
  vecMatIdx[0].push_back(2);
  vecTexIdx[0].push_back(2);
  shapes.push_back( vtkm::CELL_SHAPE_QUAD);
  numindices.push_back( 4);
  pts1.push_back( vec3(555,0,0)/555.0);
  pts1.push_back( vec3(555,555,0)/555.0);
  pts1.push_back( vec3(555,555,555)/555.0);
  pts1.push_back( vec3(555,0,555)/555.0);
  vecField.push_back( 0);
  vecField.push_back( 0);
  vecField.push_back( 0);
  vecField.push_back( 0);


  vecMatIdx[0].push_back(0);
  vecTexIdx[0].push_back(0);
  shapes.push_back( vtkm::CELL_SHAPE_QUAD);
  numindices.push_back( 4);
  pts1.push_back( vec3(0,0,0));
  pts1.push_back( vec3(0,555,0)/555.0);
  pts1.push_back( vec3(0,555,555)/555.0);
  pts1.push_back( vec3(0,0,555)/555.0);
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));


  //xz_rect
  vecMatIdx[0].push_back(3);
  vecTexIdx[0].push_back(3);
  shapes.push_back( vtkm::CELL_SHAPE_QUAD);
  numindices.push_back( 4);
  pts1.push_back( vec3(213,554,227)/555.0);
  pts1.push_back( vec3(343,554,227)/555.0);
  pts1.push_back( vec3(343,554,332)/555.0);
  pts1.push_back( vec3(213,554,332)/555.0);
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));


  vecMatIdx[0].push_back(1);
  vecTexIdx[0].push_back(1);
  shapes.push_back( vtkm::CELL_SHAPE_QUAD);
  numindices.push_back( 4);
  pts1.push_back( vec3(0,555,0)/555.0);
  pts1.push_back( vec3(555,555,0)/555.0);
  pts1.push_back( vec3(555,555,555)/555.0);
  pts1.push_back( vec3(0,555,555)/555.0);
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));


  vecMatIdx[0].push_back(1);
  vecTexIdx[0].push_back(1);
  shapes.push_back( vtkm::CELL_SHAPE_QUAD);
  numindices.push_back( 4);
  pts1.push_back( vec3(0,0,0));
  pts1.push_back( vec3(555,0,0)/555.0);
  pts1.push_back( vec3(555,0,555)/555.0);
  pts1.push_back( vec3(0,0,555)/555.0);
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));


  //xy_rect
  vecMatIdx[0].push_back(1);
  vecTexIdx[0].push_back(1);
  shapes.push_back( vtkm::CELL_SHAPE_QUAD);
  numindices.push_back( 4);
  pts1.push_back( vec3(0,0,555)/555.0);
  pts1.push_back( vec3(555,0,555)/555.0);
  pts1.push_back( vec3(555,555,555)/555.0);
  pts1.push_back( vec3(0,555,555)/555.0);
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));
  vecField.push_back( float(numindices.size()-1));



//  //small box
  //xy
  vtkm::Vec<vec3,4> pts;
  pts[0] = vec3(0,0,165);
  pts[1] = vec3(165,0,165);
  pts[2] = vec3(165,330,165);
  pts[3] = vec3(0,330,165);
  invert(pts);
  buildQuad(pts,
            shapes,
            numindices,
            pts1,
            vecField,
            vecMatIdx[0],
            vecTexIdx[0],
      1,1);

  pts[0] = vec3(0,0,0);
  pts[1] = vec3(165,0,0);
  pts[2] = vec3(165,330,0);
  pts[3] = vec3(0,330,0);
  invert(pts);
  buildQuad(pts,
            shapes,
            numindices,
            pts1,
            vecField,
            vecMatIdx[0],
            vecTexIdx[0],
      1,1);


  //yz
  pts[0] = vec3(165,0,0);
  pts[1] = vec3(165,330,0);
  pts[2] = vec3(165,330,165);
  pts[3] = vec3(165, 0, 165);
  invert(pts);
  buildQuad(pts,
            shapes,
            numindices,
            pts1,
            vecField,
            vecMatIdx[0],
            vecTexIdx[0],
      1,1);


  pts[0] = vec3(0,0,0);
  pts[1] = vec3(0,330,0);
  pts[2] = vec3(0,330,165);
  pts[3] = vec3(0, 0, 165);
  invert(pts);
  buildQuad(pts,
            shapes,
            numindices,
            pts1,
            vecField,
            vecMatIdx[0],
            vecTexIdx[0],
      1,1);



  //xz_rect
  pts[0] = vec3(0,333,0);
  pts[1] = vec3(165,330,0);
  pts[2] = vec3(165,330,165);
  pts[3] = vec3(0, 330, 165);
  invert(pts);
  buildQuad(pts,
            shapes,
            numindices,
            pts1,
            vecField,
            vecMatIdx[0],
            vecTexIdx[0],
      1,1);

  pts[0] = vec3(0,0,0);
  pts[1] = vec3(165,0,0);
  pts[2] = vec3(165,0,165);
  pts[3] = vec3(0, 0, 165);
  invert(pts);
  buildQuad(pts,
            shapes,
            numindices,
            pts1,
            vecField,
            vecMatIdx[0],
            vecTexIdx[0],
      1,1);


  //    //sphere
  vec3 sphere_center(-335,90, 290);
  vtkm::Float32 sphere_radii = 90;

  vecMatIdx[1].push_back(4);
  vecTexIdx[1].push_back(0);
  shapes.push_back( vtkm::CELL_SHAPE_VERTEX);
  numindices.push_back( 1);
  pts1.push_back( sphere_center/555.0);
  vecField.push_back( float(numindices.size()-1));


  vec3 box_center(135,90, 290);
  vec3 near(box_center[0] - sphere_radii, 0, box_center[2] - sphere_radii);
  vec3 far(box_center[0] + sphere_radii, 180,box_center[2] + sphere_radii);
  buildBox(near,far,
           shapes, numindices,
           pts1, vecField,
           vecMatIdx[0], vecTexIdx[0],
           1,1);

  //------------------------------------------------------------------------

  near = vec3(50, 0, 50);
  far = vec3(450, 100, 100);

  buildBox(near,far,
           shapes, numindices,
           pts1, vecField,
           vecMatIdx[0], vecTexIdx[0],
           1,1);



  for (auto &val: vecField)
    val /= float(vecField.size());

  matIdx[0] = vtkm::cont::make_ArrayHandle(vecMatIdx[0], vtkm::CopyFlag::On);
  matIdx[1] = vtkm::cont::make_ArrayHandle(vecMatIdx[1], vtkm::CopyFlag::On);
  texIdx[0] = vtkm::cont::make_ArrayHandle(vecTexIdx[0], vtkm::CopyFlag::On);
  texIdx[1] = vtkm::cont::make_ArrayHandle(vecTexIdx[1], vtkm::CopyFlag::On);

  coord.SetData(vtkm::cont::make_ArrayHandle(pts1, vtkm::CopyFlag::On));

  conn.Allocate(pts1.size());
  vtkm::cont::ArrayCopy(vtkm::cont::ArrayHandleCounting<vtkm::Id>(0,1, pts1.size()), conn);

  vtkm::cont::DataSetBuilderExplicit dsb;
  auto arr = coord.GetData().Cast<vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,3>>>();
  ds = dsb.Create(arr,
                  vtkm::cont::make_ArrayHandle(shapes, vtkm::CopyFlag::On),
                  vtkm::cont::make_ArrayHandle(numindices, vtkm::CopyFlag::On),
                  conn, "coords");

  field = vtkm::cont::make_ArrayHandle(vecField, vtkm::CopyFlag::On);
  vtkm::cont::Field pfield(
    "point_var",
    vtkm::cont::Field::Association::POINTS,
    field);

  ds.AddField(pfield);
  return ds;
}

void CornellBox::extract()
{

  vtkm::rendering::pathtracing::SphereExtractor sphereExtractor;
  sphereExtractor.ExtractCells(ds.GetCellSet(), 90.0/555.0);
  SphereIds = sphereExtractor.GetPointIds();
  SphereRadii = sphereExtractor.GetRadii();
  ShapeOffset = ds.GetCellSet().Cast<vtkm::cont::CellSetExplicit<>>().GetOffsetsArray(vtkm::TopologyElementTagPoint(), vtkm::TopologyElementTagCell());
  for (int i=0; i<SphereIds.GetNumberOfValues(); i++){
    std::cout << SphereIds.ReadPortal().Get(i) << std::endl;
    std::cout << SphereRadii.ReadPortal().Get(i) << std::endl;

  }
  vtkm::rendering::raytracing::QuadExtractor quadExtractor;
  quadExtractor.ExtractCells(ds.GetCellSet());
  QuadIds = quadExtractor.GetQuadIds();

}
