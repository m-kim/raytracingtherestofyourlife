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

#include <raytracing/RayTracerNormals.h>
#include "MapperPathTracer.h"

#include <fstream>
#include "CornellBox.h"

#include <vtkm/rendering/MapperRayTracer.h>
#include "MapperQuad.h"
#include "MapperQuadNormals.h"
#include "MapperQuadAlbedo.h"
#include <vtkm/rendering/Scene.h>
#include "View3D.h"
#include "PAVE.h"

// EXAMPLE CALL
// For buffers:
// ./CornellBox_CUDA -x 256 -y 256 -samplecount 300 -raydepth 20 -hemisphere -direct
// For path traced images:
//./CornellBox_CUDA -x 256 -y 256 -samplecount 300 -raydepth 20 -hemisphere
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
std::vector<double> norm_color_range( std::vector<double> color_vals){
    std::vector<double> normalized_colors = color_vals;
    for(int i=0; i<color_vals.size(); i++)
        normalized_colors[i] = color_vals[i]/255.0;
    return normalized_colors;
}
void runRay(int nx, int ny, int samplecount, int depthcount,
              vtkm::rendering::Canvas &canvas, vtkm::rendering::Camera &cam)
{
  CornellBox cb;
  path::rendering::MapperQuad mapper;
  auto ds = cb.buildDataSet();
  vtkm::rendering::Scene scene;

  std::vector<double> colors ={113, 31 ,30,
                               84, 23, 23,
                              132 ,128 ,126,
                               48, 93 ,53,
                              69 ,63 ,59,
                              54 ,47 ,43,
                              109,104, 102,
                              251, 251, 251,
                              1, 1, 1,
                              35, 66 ,38,
                              91 ,86, 84,
                              25, 32, 21};
  colors = norm_color_range(colors);


  std::vector<double> red = {1, 0, 0};
  std::vector<double> white = {1, 1, 1};
  std::vector<double> green = {0, 0, 1};

  std::vector<double> blue = {0, 0, 1};
  //0,1,1 0,1,0 1,1,0
  std::vector<double> lamb1 = {0, 1, 1};
  std::vector<double> lamb2 = {0, 1, 0};
  std::vector<double> lamb3 = {1, 1, 0};

  std::vector<double> c1 = {0.65, 0.05, 0.05}; //red
  std::vector<double> c2 = {0.73, 0.73, 0.73}; //white
  std::vector<double> c3 = {0.12, 0.45, 0.15}; //green
  std::vector<double> fill1 = c2;//{.15, .15, .15};

  std::vector<double> pallet;
  int num_quads = 12+6+6;
  int num_colors = 3;
  pallet.reserve(num_quads*num_colors);


  pallet.insert(pallet.end(), c3.begin(), c3.end()); //green
  pallet.insert(pallet.end(), c1.begin(), c1.end()); //red
  pallet.insert(pallet.end(), fill1.begin(), fill1.end()); //light
  for (int i=0; i<num_quads - 3; i++)
    pallet.insert(pallet.end(), c2.begin(), c2.end()); //white


  std::vector<double> alpha(num_quads);
  for (int i=0; i<alpha.size(); i++)
      alpha[i] = 1.0;

  vtkm::cont::ColorTable ct_12_quad("pallet_color_table",
                            vtkm::cont::ColorSpace::RGB,
                            vtkm::Vec<double,3>(0,0,0),
                            pallet, alpha);
  scene.AddActor(vtkm::rendering::Actor(
                   ds.GetCellSet(),
                   ds.GetCoordinateSystem(),
                   ds.GetField("point_var"),
                   ct_12_quad));
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

void runAlbedo(int nx, int ny, int samplecount, int depthcount,
              vtkm::rendering::Canvas &canvas, vtkm::rendering::Camera &cam)
{
  CornellBox cb;
  path::rendering::MapperQuadAlbedo mapper;
  auto ds = cb.buildDataSet();
  vtkm::rendering::Scene scene;

  std::vector<double> c1 = {0.65, 0.05, 0.05}; //red
  std::vector<double> c2 = {0.73, 0.73, 0.73}; //white
  std::vector<double> c3 = {0.12, 0.45, 0.15}; //green
  std::vector<double> fill1 = c2;//{.15, .15, .15};

  std::vector<double> pallet;
  int num_quads = 12+6+6;
  int num_colors = 3;
  pallet.reserve(num_quads*num_colors);
  pallet.insert(pallet.end(), c3.begin(), c3.end()); //green
  pallet.insert(pallet.end(), c1.begin(), c1.end()); //red
  pallet.insert(pallet.end(), fill1.begin(), fill1.end()); //light
  for (int i=0; i<num_quads - 3; i++)
    pallet.insert(pallet.end(), c2.begin(), c2.end()); //white
  std::vector<double> alpha(num_quads); //alpha
  for (int i=0; i<alpha.size(); i++)
      alpha[i] = 1.0;

  vtkm::cont::ColorTable ct_12_quad("pallet_color_table",
                            vtkm::cont::ColorSpace::RGB,
                            vtkm::Vec<double,3>(0,0,0),
                            pallet, alpha);

  scene.AddActor(vtkm::rendering::Actor(
                   ds.GetCellSet(),
                   ds.GetCoordinateSystem(),
                   ds.GetField("point_var"),
                   ct_12_quad));
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
  //col = col / float(samplecount);
  if(col != col)
      col = 0.f;
  col = sqrt(col);
  int ir = int(255.99*col);
  int ig = int(255.99*col);
  int ib = int(255.99*col);
  fs << ir << " " << ig << " " << ib << std::endl;

}

template<typename ArrayType>
void save(std::stringstream &fnstream,
          int nx, int ny, int samplecount,
          ArrayType &cols)
{
  auto orig = fnstream.str();
  fnstream << ".pnm";
  auto fn = fnstream.str();
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
  fnstream.str(orig);
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

  std::unique_ptr<PAVE> paver[5];
  paver[0] = std::unique_ptr<PAVE>(new PAVE("direct.bp"));
  paver[1] = std::unique_ptr<PAVE>(new PAVE("depth.bp"));
  paver[2] = std::unique_ptr<PAVE>(new PAVE("normals.bp"));
  paver[3] = std::unique_ptr<PAVE>(new PAVE("albedo.bp"));
  paver[4] = std::unique_ptr<PAVE>(new PAVE("outputs.bp"));

  for (float phi=M_PI*0.25; phi<M_PI*0.5; phi += rPhi){
	std::cout << "Phi: " << phi << std::endl;
    for (float theta=0; theta<2*M_PI; theta+=rTheta){
      auto x = r * cos(theta) * sin(phi);
      auto y = r * sin(theta) * sin(phi);
      auto z = r * cos(phi);

      vec3 pos(x+278/555.0, y+278/555.0, z+278/555.0 );
      cam.SetPosition(pos);
      std::stringstream phiTheta;
      phiTheta << phi << "-" << theta;
      std::stringstream sstr;
      if (direct){
        sstr << "direct-" << phiTheta.str();
        runRay(nx,ny,samplecount, depthcount, canvas, cam);
        save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());
        paver[0]->save(phiTheta,nx,ny, samplecount, canvas.GetColorBuffer());

        sstr.str("");
        sstr << "depth-" << phiTheta.str();
        save(sstr, nx, ny, samplecount, canvas.GetDepthBuffer());

        paver[1]->save(phiTheta,nx,ny, samplecount, canvas.GetDepthBuffer());

        runNorms(nx,ny,samplecount,depthcount, canvas, cam);
        sstr.str("");
        sstr << "normals-" << phiTheta.str();
        save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());
        paver[2]->save(phiTheta,nx,ny, samplecount, canvas.GetColorBuffer());


        sstr.str("");
        runAlbedo(nx,ny,samplecount,depthcount, canvas, cam);
        sstr << "albedo-" << phiTheta.str();
        save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());
        paver[3]->save(phiTheta,nx,ny, samplecount, canvas.GetColorBuffer());
      }
      else{
        sstr << "output-" << phiTheta.str();
        runPath(nx,ny, samplecount, depthcount, canvas, cam);
        save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());
        paver[4]->save(phiTheta,nx,ny, samplecount, canvas.GetColorBuffer());
      }
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
    std::unique_ptr<PAVE> paver[3];
    if (direct){

      paver[0] = std::unique_ptr<PAVE>(new PAVE("direct.bp"));
      runRay(nx,ny,samplecount,depthcount, canvas, cam);
      std::stringstream sstr;
      sstr << "direct";
      save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());
      sstr.str("");
      sstr << "direct";
      paver[0]->save(sstr,nx,ny, samplecount, canvas.GetColorBuffer());

      sstr.str("");
      sstr << "depth";
      save(sstr, nx, ny, samplecount, canvas.GetDepthBuffer());
      sstr.str("");
      sstr << "depth";
      paver[1]  = std::unique_ptr<PAVE>(new PAVE("depth.bp"));
      paver[1]->save(sstr, nx,ny, samplecount, canvas.GetDepthBuffer());
      runNorms(nx,ny,samplecount,depthcount, canvas, cam);
      sstr.str("");
      sstr << "normals";
      save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());
      sstr.str("");
      sstr << "normals";
      paver[2] = std::unique_ptr<PAVE>(new PAVE("normals.bp"));
      paver[2]->save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());

      sstr.str("");
      runAlbedo(nx,ny,samplecount,depthcount, canvas, cam);
      sstr << "albedo";
      save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());

    }
    else{
      runPath(nx,ny, samplecount, depthcount, canvas, cam);
      std::stringstream sstr;
      sstr << "output";
      save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());

      PAVE paver("output.bp");
      paver.save(sstr,  nx, ny, samplecount, canvas.GetColorBuffer());

    }
  }
}

