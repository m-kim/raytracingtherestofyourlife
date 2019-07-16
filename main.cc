//==================================================================================================
// Written in 2019 by Mark Kim
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is distributed
// without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication along
// with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==================================================================================================

#include <iostream>
#include <limits>
#include <vector>
#include <tuple>
#include <sstream>
#include <vtkm/cont/ArrayHandleConstant.h>
#include "pathtracing/Camera.h"
#include <vtkm/cont/internal/DeviceAdapterAlgorithmGeneral.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleCompositeVector.h>
#include <omp.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/worklet/Invoker.h>
#include <vtkm/rendering/raytracing/BoundingVolumeHierarchy.h>
#include <vtkm/cont/ArrayHandleExtractComponent.h>
#include <adios2.h>

#include <raytracing/RayTracerNormals.h>
#include "MapperPathTracer.h"

#include <fstream>
#include "CornellBox.h"

#include <vtkm/rendering/MapperRayTracer.h>
#include "MapperQuad.h"
#include "MapperQuadNormals.h"
#include <vtkm/rendering/Scene.h>
#include "View3D.h"

using ArrayType = vtkm::cont::ArrayHandle<vec3>;

template<typename VecType>
inline VecType de_nan(const VecType& c) {
    auto temp = c;
    if (!(temp[0] == temp[0])) temp[0] = 0;
    if (!(temp[1] == temp[1])) temp[1] = 0;
    if (!(temp[2] == temp[2])) temp[2] = 0;
    return temp;
}


const auto
parse(int argc, char **argv){
  int x = 128;
  int y = 128;
  int s = 10;
  int depth = 5;

  bool hemi = false;
  bool direct = false;
  for (int i=1; i<argc; i++){
    if (!strcmp(argv[i], "-x")){
      if (i+1 < argc){
        x = atoi(argv[i+1]);
        i += 1;
      }

    }
    else if (!strcmp(argv[i], "-y")){
      if (i+1 < argc){
        y = atoi(argv[i+1]);
        i += 1;
      }
    }
    else if (!strcmp(argv[i], "-samplecount")){
      if (i+1 < argc){
        s = atoi(argv[i+1]);
        i += 1;
      }
    }
    else if (!strcmp(argv[i], "-raydepth")){
      if (i+1 < argc){
        depth = atoi(argv[i+1]);
        i += 1;
      }
    }
    else if (!strcmp(argv[i], "-hemisphere"))
    {
      hemi = true;
    }
    else if(!strcmp(argv[i], "-direct"))
      direct = true;
  }

  return std::make_tuple(x,y, s, depth, hemi, direct);
}
void runRay(int nx, int ny, int samplecount, int depthcount,
              vtkm::rendering::Canvas &canvas, vtkm::rendering::Camera &cam)
{
  CornellBox cb;
  path::rendering::MapperQuad mapper;
  auto ds = cb.buildDataSet();
  vtkm::rendering::Scene scene;

  scene.AddActor(vtkm::rendering::Actor(
                   ds.GetCellSet(),
                   ds.GetCoordinateSystem(),
                   ds.GetField("point_var"),
                   vtkm::cont::ColorTable{vtkm::cont::ColorTable::Preset::COOL_TO_WARM_EXTENDED}));
  vtkm::rendering::Color background(0,0,0, 1.0f);
  vtkm::rendering::Color foreground(1,1,1, 1.0f);
  vtkm::rendering::View3D view(scene, mapper, canvas, cam, background, foreground);

  view.Initialize();
  view.Paint();

}
void runNorms(int nx, int ny, int samplecount, int depthcount,
              vtkm::rendering::Canvas &canvas, vtkm::rendering::Camera &cam)
{
  CornellBox cb;
  path::rendering::MapperQuadNormals mapper;
  auto ds = cb.buildDataSet();
  vtkm::rendering::Scene scene;

  scene.AddActor(vtkm::rendering::Actor(
                   ds.GetCellSet(),
                   ds.GetCoordinateSystem(),
                   ds.GetField("point_var"),
                   vtkm::cont::ColorTable{vtkm::cont::ColorTable::Preset::COOL_TO_WARM_EXTENDED}));
  vtkm::rendering::Color background(0,0,0, 1.0f);
  vtkm::rendering::Color foreground(1,1,1, 1.0f);
  vtkm::rendering::View3D view(scene, mapper, canvas, cam, background, foreground);

  view.Initialize();
  view.Paint();

}
void runPath(int nx, int ny, int samplecount, int depthcount,
              vtkm::rendering::Canvas &canvas, vtkm::rendering::Camera &cam)
{
  using MyAlgos = details::PathAlgorithms<vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>, VTKM_DEFAULT_DEVICE_ADAPTER_TAG>;
  using StorageTag = vtkm::cont::StorageTagBasic;
  using Device = VTKM_DEFAULT_DEVICE_ADAPTER_TAG;

  CornellBox cb;

  auto ds = cb.buildDataSet();

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 2>> uvs;
  uvs.Allocate(nx*ny);

  vtkm::rendering::MapperPathTracer mapper(samplecount,
                                           depthcount,
                                           cb.matIdx,
                                           cb.texIdx,
                                           cb.matType,
                                           cb.texType,
                                           cb.tex);

  mapper.SetCanvas(&canvas);



  vtkm::cont::Field field;
  vtkm::cont::ColorTable ct;
  vtkm::Range sr;
  mapper.RenderCells(ds.GetCellSet(0),
                     cb.coord,
                     field,
                     ct,
                     cam,
                     sr);


}

template<typename ValueType>
void save(std::fstream &fs,
          int samplecount,
          ValueType &col);

template<>
void save(std::fstream &fs,
          int samplecount,
          vtkm::Vec<vtkm::Float32,4> &col)
{
  col = de_nan(col);
  col = col / float(samplecount);
  col[0] = sqrt(col[0]);
  col[1] = sqrt(col[1]);
  col[2] = sqrt(col[2]);
  int ir = int(255.99*col[0]);
  int ig = int(255.99*col[1]);
  int ib = int(255.99*col[2]);
  fs << ir << " " << ig << " " << ib << std::endl;

}
template<>
void save(std::fstream &fs,
          int samplecount,
          vtkm::Float32 &col)
{
  col = sqrt(col);
  int ir = int(255.99*col);
  int ig = int(255.99*col);
  int ib = int(255.99*col);
  fs << ir << " " << ig << " " << ib << std::endl;

}

template<typename ArrayType>
void save(std::string fn,
          int nx, int ny, int samplecount,
          ArrayType &cols)
{
  std::fstream fs;
  fs.open(fn.c_str(), std::fstream::out);
  if (fs.is_open()){
    fs << "P3\n" << nx << " "  << ny << " 255" << std::endl;
    for (int i=0; i<cols.GetNumberOfValues(); i++){
      auto col = cols.GetPortalConstControl().Get(i);
      if (col != col)
        col = 0.0f;
      save(fs, samplecount, col);
    }
    fs.close();
  }
  else
    std::cout << "Couldn't save pnm." << std::endl;
//  std::vector<std::uint8_t> PngBuffer;
}

template<typename ArrayType>
void saveADIOS(std::string fn,
               int nx, int ny, int samplecount,
               ArrayType &cols);

template<>
void saveADIOS(std::string fn,
               int nx, int ny, int samplecount,
               vtkm::cont::ArrayHandle<vtkm::Float32> &cols)
{
  using ArrayType = vtkm::Float32;
  adios2::ADIOS adios(adios2::DebugON);
  adios2::IO bpIO = adios.DeclareIO("BPFile_N2N");

  adios2::Variable<ArrayType> bpOut = bpIO.DefineVariable<ArrayType>(
        "pnms", {}, {}, {static_cast<std::size_t>(nx*ny)}, adios2::ConstantDims);

  adios2::Engine writer = bpIO.Open(fn, adios2::Mode::Write);

  auto *ptr = cols.GetStorage().GetArray();
  writer.Put<vtkm::Float32>(bpOut, ptr );
  writer.Close();

}
template<>
void saveADIOS(std::string fn,
               int nx, int ny, int samplecount,
               vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,4>> &cols)
{
  using OutArrayType = vtkm::Float32;
  adios2::ADIOS adios(adios2::DebugON);
  adios2::IO bpIO = adios.DeclareIO("BPFile_N2N");

  adios2::Variable<OutArrayType> bpOut = bpIO.DefineVariable<OutArrayType>(
        "pnms", {}, {}, {static_cast<std::size_t>(nx*ny*4)}, adios2::ConstantDims);

  adios2::Engine writer = bpIO.Open(fn, adios2::Mode::Write);

  std::vector<vtkm::Float32> arrayOut(cols.GetNumberOfValues()*4);
  for (int i=0; i<cols.GetNumberOfValues(); i++){
    auto col = cols.GetPortalConstControl().Get(i);
    arrayOut[i*4] = col[0];
    arrayOut[i*4 + 1] = col[1];
    arrayOut[i*4 + 2] = col[2];
    arrayOut[i*4 + 3] = col[3];
  }
  writer.Put<OutArrayType>(bpOut, arrayOut.data() );
  writer.Close();

}

void generateHemisphere(int nx, int ny, int samplecount, int depthcount, bool direct)
{
  vtkm::rendering::CanvasRayTracer canvas(nx,ny);
  vtkm::rendering::Camera cam;
  cam.SetClippingRange(01.f, 5.f);
  cam.SetPosition(vec3(278/555.0,278/555.0,-800/555.0));
  cam.SetFieldOfView(40.f);
  cam.SetViewUp(vec3(0,1,0));
  cam.SetLookAt(vec3(278/555.0,278/555.0,278/555.0));

  int numPhi = 30;
  int numTheta = 30;

  float rTheta = (2.0*M_PI)/float(numTheta);
  float rPhi = (M_PI/2.0)/float(numPhi);

  float r = -1078/555.0;
  for (float phi=M_PI*0.25; phi<M_PI*0.5; phi += rPhi){
	std::cout << "Phi: " << phi << std::endl;
    for (float theta=0; theta<2*M_PI; theta+=rTheta){
      auto x = r * cos(theta) * sin(phi);
      auto y = r * sin(theta) * sin(phi);
      auto z = r * cos(phi);

      vec3 pos(x+278/555.0, y+278/555.0, z+278/555.0 );
      cam.SetPosition(pos);
      std::stringstream sstr;
      if (direct){
        sstr << "direct-" << phi << "-" << theta << ".pnm";
        runRay(nx,ny,samplecount, depthcount, canvas, cam);
        save(sstr.str(), nx, ny, samplecount, canvas.GetColorBuffer());
        sstr.str("");
        sstr << "depth-" << phi << "-" << theta << ".pnm";
        save(sstr.str(), nx, ny, samplecount, canvas.GetDepthBuffer());
        runNorms(nx,ny,samplecount,depthcount, canvas, cam);
        sstr.str("");
        sstr << "normals-" << phi << "-" << theta << ".pnm";
        save(sstr.str(), nx, ny, samplecount, canvas.GetColorBuffer());

      }
      else{
        sstr << "output-" << phi << "-" << theta << ".pnm";
        runPath(nx,ny, samplecount, depthcount, canvas, cam);
      }
      save(sstr.str(), nx, ny, samplecount, canvas.GetColorBuffer());
    }
  }
}
int main(int argc, char *argv[]) {

  const auto tup = parse(argc, argv);
  const int nx = std::get<0>(tup);
  const int ny = std::get<1>(tup);
  const int samplecount = std::get<2>(tup);
  const int depthcount = std::get<3>(tup);
  const bool hemi = std::get<4>(tup);
  const bool direct = std::get<5>(tup);

  if (hemi)
    generateHemisphere(nx,ny, samplecount, depthcount, direct);
  else
  {
    vtkm::rendering::CanvasRayTracer canvas(nx,ny);
    vtkm::rendering::Camera cam;
    cam.SetClippingRange(0.1f, 5.f);
    cam.SetPosition(vec3(278/555.0,278/555.0,-800/555.0));
    cam.SetFieldOfView(40.);
    cam.SetViewUp(vec3(0,1,0));
    cam.SetLookAt(vec3(278/555.0,278/555.0,278/555.0));
    if (direct){
      runRay(nx,ny,samplecount,depthcount, canvas, cam);
      std::stringstream sstr;
      sstr << "direct.pnm";
      save(sstr.str(), nx, ny, samplecount, canvas.GetColorBuffer());
      saveADIOS("direct.bp", nx,ny, samplecount, canvas.GetColorBuffer());
      sstr.str("depth.pnm");
      save(sstr.str(), nx, ny, samplecount, canvas.GetDepthBuffer());
      runNorms(nx,ny,samplecount,depthcount, canvas, cam);
      sstr << "normals.pnm";
      save(sstr.str(), nx, ny, samplecount, canvas.GetColorBuffer());

      saveADIOS("depth.bp", nx,ny, samplecount, canvas.GetDepthBuffer());
    }
    else{
      runPath(nx,ny, samplecount, depthcount, canvas, cam);
      std::stringstream sstr;
      sstr << "output.pnm";
      save(sstr.str(), nx, ny, samplecount, canvas.GetColorBuffer());

    }
  }
}

