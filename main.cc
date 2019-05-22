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
#include "Camera.h"
#include <vtkm/cont/internal/DeviceAdapterAlgorithmGeneral.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleCompositeVector.h>
#include <vtkm/rendering/raytracing/Ray.h>
#include <vtkm/rendering/raytracing/RayOperations.h>
#include <omp.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/worklet/Invoker.h>
#include <vtkm/rendering/raytracing/BoundingVolumeHierarchy.h>
#include <vtkm/cont/ArrayHandleExtractComponent.h>

#include "MapperPathTracer.h"

#include <fstream>
#include "CornellBox.h"

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
  }

  return std::make_tuple(x,y, s, depth, hemi);
}

void run(int nx, int ny, int samplecount, int depthcount,
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
    fs.close();
  }
  else
    std::cout << "Couldn't save pnm." << std::endl;
//  std::vector<std::uint8_t> PngBuffer;
}

void generateHemisphere(int nx, int ny, int samplecount, int depthcount)
{
  vtkm::rendering::CanvasRayTracer canvas(nx,ny);
  vtkm::rendering::Camera cam;
  cam.SetPosition(vec3(278,278,-800));
  cam.SetFieldOfView(40.f);
  cam.SetViewUp(vec3(0,1,0));
  cam.SetLookAt(vec3(278,278,278));

  int numPhi = 30;
  int numTheta = 30;

  float rTheta = (2.0*M_PI)/float(numTheta);
  float rPhi = (M_PI/2.0)/float(numPhi);

  float r = -1078;
  for (int i=0; i<numTheta; i++){
    for (int j=0; j<numPhi; j++){
      auto x = r * cos(i * rTheta) * sin(j * rPhi);
      auto y = r * sin(i * rTheta) * sin(j * rPhi);
      auto z = r * cos(i*rPhi);

      vec3 pos(x+278, y+278, z+278 );
      cam.SetPosition(pos);
      run(nx,ny, samplecount, depthcount, canvas, cam);
      std::stringstream sstr;
      sstr << "output-" << i*rTheta << "-" << j*rPhi << ".pnm";
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
  if (hemi)
    generateHemisphere(nx,ny, samplecount, depthcount);
  else
  {
    vtkm::rendering::CanvasRayTracer canvas(nx,ny);
    vtkm::rendering::Camera cam;
    cam.SetPosition(vec3(278,278,-800));
    cam.SetFieldOfView(40.);
    cam.SetViewUp(vec3(0,1,0));
    cam.SetLookAt(vec3(278,278,278));
    run(nx,ny, samplecount, depthcount, canvas, cam);
    std::stringstream sstr;
    sstr << "output.pnm";
    save(sstr.str(), nx, ny, samplecount, canvas.GetColorBuffer());

  }
}

