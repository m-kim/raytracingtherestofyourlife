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
#include <vtkm/cont/Initialize.h>
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
#include <vtkm/cont/Invoker.h>
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
#include <iomanip>
#include <signal.h>
#include <vtkm/cont/Timer.h>

// EXAMPLE CALL
// For buffers:
// ./CornellBox_CUDA -x 256 -y 256 -samplecount 300 -raydepth 20 -hemisphere -direct
// For path traced images:
//./CornellBox_CUDA -x 256 -y 256 -samplecount 300 -raydepth 20 -hemisphere
using ArrayType = vtkm::cont::ArrayHandle<vec3>;
std::unique_ptr<CornellBox> cb;



const auto
parse(int argc, char **argv){
  int x = 128;
  int y = 128;
  int s = 10;
  int phiCount = 15;
  int thetaCount = 15;

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
    else if(!strcmp(argv[i], "-phicount")){
      if (i+1 < argc){
        phiCount = atoi(argv[i+1]);
        i += 1;
      }
    }
    else if(!strcmp(argv[i], "-thetacount")){
      if (i+1 < argc){
        thetaCount = atoi(argv[i+1]);
        i += 1;
      }
    }
  }

  return std::make_tuple(x,y, s, depth, hemi, direct, phiCount, thetaCount);
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
  path::rendering::MapperQuad mapper;
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
                   cb->ds.GetCellSet(),
                   cb->ds.GetCoordinateSystem(),
                   cb->ds.GetField("point_var"),
                   ct_12_quad));
  vtkm::rendering::Color background(0,0,0, 1.0f);
  vtkm::rendering::Color foreground(1,1,1, 1.0f);
  vtkm::rendering::pathtracing::View3D view(scene, mapper, canvas, cam, background, foreground);

  view.Initialize();
  view.Paint();

}
void runNorms(int nx, int ny, int samplecount, int depthcount,
              vtkm::rendering::Canvas &canvas, vtkm::rendering::Camera &cam)
{
  path::rendering::MapperQuadNormals mapper;
  vtkm::rendering::Scene scene;

  scene.AddActor(vtkm::rendering::Actor(
                   cb->ds.GetCellSet(),
                   cb->ds.GetCoordinateSystem(),
                   cb->ds.GetField("point_var"),
                   vtkm::cont::ColorTable{vtkm::cont::ColorTable::Preset::COOL_TO_WARM_EXTENDED}));
  vtkm::rendering::Color background(0,0,0, 1.0f);
  vtkm::rendering::Color foreground(1,1,1, 1.0f);
  vtkm::rendering::pathtracing::View3D view(scene, mapper, canvas, cam, background, foreground);

  view.Initialize();
  view.Paint();

}

void runAlbedo(int nx, int ny, int samplecount, int depthcount,
              vtkm::rendering::Canvas &canvas, vtkm::rendering::Camera &cam)
{
  path::rendering::MapperQuadAlbedo mapper;
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
                   cb->ds.GetCellSet(),
                   cb->ds.GetCoordinateSystem(),
                   cb->ds.GetField("point_var"),
                   ct_12_quad));
  vtkm::rendering::Color background(0,0,0, 1.0f);
  vtkm::rendering::Color foreground(1,1,1, 1.0f);
  vtkm::rendering::pathtracing::View3D view(scene, mapper, canvas, cam, background, foreground);

  view.Initialize();
  view.Paint();

}

struct NormalizeFunctor
{
  NormalizeFunctor(int sc)
  {
    samplecount = static_cast<float>(sc);
  }

  VTKM_EXEC_CONT
  vtkm::Vec<vtkm::Float32,4> de_nan(const vtkm::Vec<vtkm::Float32,4>& c) const
  {
      auto temp = c;
      if (!(temp[0] == temp[0])) temp[0] = 0;
      if (!(temp[1] == temp[1])) temp[1] = 0;
      if (!(temp[2] == temp[2])) temp[2] = 0;
      return temp;
  }
  VTKM_EXEC_CONT
  vtkm::Vec<vtkm::Float32,4> de_nan(const vtkm::Float32 & c) const
  {
      auto temp = c;
      if (!(temp == temp)) temp = 0;
      return temp;
  }

  template <typename T>
  VTKM_EXEC_CONT T operator()(const T& x, const T& y) const
  {
    auto tmp = de_nan(x);
    tmp = vtkm::Sqrt(tmp/samplecount);

    return tmp;
  }

  float samplecount;
};

void runPath(int nx, int ny, int samplecount, int depthcount,
              vtkm::rendering::Canvas &canvas, vtkm::rendering::Camera &cam)
{


  vtkm::rendering::MapperPathTracer mapper(samplecount,
                                           depthcount,
                                           cb->matIdx,
                                           cb->texIdx,
                                           cb->matType,
                                           cb->texType,
                                           cb->tex);

  mapper.SetCanvas(&canvas);



  vtkm::cont::Field field;
  vtkm::cont::ColorTable ct;
  vtkm::Range sr;
  mapper.RenderCells(cb->ds.GetCellSet(),
                     cb->coord,
                     field,
                     ct,
                     cam,
                     sr);


  vtkm::cont::Algorithm::Transform(
              canvas.GetColorBuffer(),
              canvas.GetColorBuffer(),
              canvas.GetColorBuffer(),
              NormalizeFunctor(samplecount));

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
  if ((col[0] != col[0]) || (col[1] != col[1]) || (col[2] != col[2])){
    col = 0.0f;
  }

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
      auto col = cols.ReadPortal().Get(i);

      save(fs, samplecount, col);
    }
    fs.close();
  }
  else
    std::cout << "Couldn't save pnm." << std::endl;
//  std::vector<std::uint8_t> PngBuffer;
  fnstream.str(orig);
}


void generate(vtkm::rendering::Camera &cam,
              vtkm::rendering::CanvasRayTracer &canvas,
              int nx, int ny,
              int samplecount,
              int depthcount,
              bool direct,
              float phi, float theta)
{
  std::stringstream phiTheta;
//  phiTheta << std::fixed << std::setw(4) << std::setprecision(4) << x << "-";
//  phiTheta << std::fixed << std::setw(4) << std::setprecision(4) <<  y <<"-";
//  phiTheta << std::fixed << std::setw(4) << std::setprecision(4) <<  z;
    phiTheta << std::fixed << std::setw(4) << std::setprecision(4) << phi << "-";
    phiTheta << std::fixed << std::setw(4) << std::setprecision(4) <<  theta;
  std::stringstream sstr;
  if (direct){
    sstr << "direct-" << phiTheta.str();
    runRay(nx,ny,samplecount, depthcount, canvas, cam);
    save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());
    
    sstr.str("");
    sstr << "depth-" << phiTheta.str();
    save(sstr, nx, ny, samplecount, canvas.GetDepthBuffer());

    
    runNorms(nx,ny,samplecount,depthcount, canvas, cam);
    sstr.str("");
    sstr << "normals-" << phiTheta.str();
    save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());
    

    sstr.str("");
    runAlbedo(nx,ny,samplecount,depthcount, canvas, cam);
    sstr << "albedo-" << phiTheta.str();
    save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());
  }
  else{
    sstr << "output-" << phiTheta.str();
    runPath(nx,ny, samplecount, depthcount, canvas, cam);
    save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());
    
  }
}
void fibonacciHemisphere(int sampleCount,
                         int nx, int ny,
                         int samplecount,
                         int depthcount,
                         bool direct)
{
//  def fibonacci_sphere(samples=1,randomize=True):
//      rnd = 1.
//      if randomize:
//          rnd = random.random() * samples

//      points = []
//      offset = 2./samples
//      increment = math.pi * (3. - math.sqrt(5.));

//      #for i in range(samples):
//      #    y = ((i * offset) - 1) + (offset / 2);
//      for i in range(samples / 2, samples):
//          y = ((i * offset) - 1) + (offset / 2)
//          r = math.sqrt(1 - pow(y,2))

//          phi = ((i + rnd) % samples) * increment

//          x = math.cos(phi) * r
//          z = math.sin(phi) * r

//          points.append([x,y,z])

//      return points

  vtkm::rendering::CanvasRayTracer canvas(nx,ny);
  vtkm::rendering::Camera cam;
  cam.SetClippingRange(01.f, 5.f);
  cam.SetPosition(vec3(278/555.0,278/555.0,-800/555.0));
  cam.SetFieldOfView(40.f);
  cam.SetViewUp(vec3(0,1,0));
  cam.SetLookAt(vec3(278/555.0,278/555.0,278/555.0));

  int rnd = rand() % sampleCount;
  vtkm::Float32 offset = 2./sampleCount;
  auto increment = vtkm::Pi() * (3. - vtkm::Sqrt(5.0));
  float zdepth = -1078/555.0;

  //for (int i=0; i<sampleCount ;i++){
  for (int i=0; i<sampleCount ;i++){
    vtkm::Float32  z = ((i * offset) - 1) + (offset /2);
    vtkm::Float32 r = vtkm::Sqrt(1 - pow(z,2));
    auto phi = ((i + rnd) %sampleCount) * increment;
    vtkm::Float32 x = vtkm::Cos(phi) * r;
    vtkm::Float32 y = vtkm::Sin(phi) * r;
    std::cout << x << " " << y << " " << z << std::endl;
    if (z < 0){
      vec3 pos(x+278/555.0, y+278/555.0, z+278/555.0 );
      cam.SetPosition(pos);
      generate(cam,
               canvas,
               nx, ny,
               samplecount,
               depthcount,
               direct,
               0,i);

    }


  }

}

volatile sig_atomic_t flag = 0;
void my_function(int sig){ // can be called asynchronously
  flag = 1; // set flag
}

bool generateHemisphere(int nx, int ny,
                        int samplecount,
                        int depthcount,
                        bool direct,
                        float phiBegin = 0,
                        float phiEnd = 1.0, //(M_PI/2.0)
                        int numPhi = 15,
                        float rPhi = 1.0/15.0,
                        float thetaBegin = 0,
                        float thetaEnd = 2*vtkm::Pi(),
                        float rTheta = 1.0/15.0
                        )
{
  vtkm::rendering::CanvasRayTracer canvas(nx,ny);
  vtkm::rendering::Camera cam;
  cam.SetClippingRange(01.f, 5.f);
  cam.SetPosition(vec3(278/555.0,278/555.0,-800/555.0));
  cam.SetFieldOfView(40.f);
  cam.SetViewUp(vec3(0,1,0));
  cam.SetLookAt(vec3(278/555.0,278/555.0,278/555.0));



  float r = -1078/555.0;


  //if (fabs(phiBegin) < 1e-6)
  //  phiBegin += rPhi;

  for (float phi=phiBegin; phi<(phiEnd - 0.5*rPhi); phi += rPhi){
    for (float theta=thetaBegin; theta<thetaEnd; theta+=rTheta){
      if (flag){
        flag = 0;
        std::cout << "release" << std::endl;
        return 1;
      }
      auto x = r * cos(theta) * sin(phi);
      auto y = r * sin(theta) * sin(phi);
      auto z = r * cos(phi);

      vec3 pos(x+278/555.0, y+278/555.0, z+278/555.0 );
      cam.SetPosition(pos);
      generate(cam,
               canvas,
               nx, ny,
               samplecount,
               depthcount,
               direct,
               phi,theta);
    }


  }



  return 0;
}
int main(int argc, char *argv[]) {

  const auto tup = parse(argc, argv);
  const int nx = std::get<0>(tup);
  const int ny = std::get<1>(tup);
  const int samplecount = std::get<2>(tup);
  const int depthcount = std::get<3>(tup);
  const bool hemi = std::get<4>(tup);
  const bool direct = std::get<5>(tup);

  vtkm::cont::Initialize();

  

  cb = std::make_unique<CornellBox>(CornellBox());
  const int phiCount = std::get<6>(tup);
  const int thetaCount = std::get<7>(tup);

  vtkm::Float64 computationTime = 0.0;
  vtkm::Float64 elapsedTime1, elapsedTime2, elapsedTime3;

  // Decompose
  vtkm::cont::Timer timer;
  timer.Start();

  cb->buildDataSet();
  if (hemi){
    float phiBegin, phiEnd;
    phiBegin = 0.0;//rPhi * static_cast<float>(rank);
    phiEnd = 1.0;//rPhi * static_cast<float>(rank + 1);

    const float thetaEnd = 2*vtkm::Pi();//rTheta * static_cast<float>(rank+1) * 2*M_PI;
    float rTheta = thetaEnd/(static_cast<float>(thetaCount));
	float thetaBegin = 0;

    float rPhi = (phiEnd - phiBegin)/float(phiCount);

    std::cout << " rPhi: " << rPhi;
    std::cout << " rTheta: " << rTheta;
    std::cout << " countTheta: " << thetaCount;
    std::cout << " countPhi: " << phiCount;
    std:: cout << " phiBegin " << phiBegin << " phiEnd: " << phiEnd;
    std::cout << " thetaBegin " << thetaBegin;
    std::cout << " thetaEnd " << thetaEnd << std::endl;

    generateHemisphere(nx,ny, samplecount, depthcount,
                       direct, phiBegin, phiEnd, phiCount, rPhi,
                       thetaBegin,
                       thetaEnd,
                       rTheta);
//    fibonacciHemisphere(10000, nx,ny, samplecount, depthcount,direct);
  }
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
      sstr << "direct";
      save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());
      sstr.str("");
      sstr << "direct";

      sstr.str("");
      sstr << "depth";
      save(sstr, nx, ny, samplecount, canvas.GetDepthBuffer());
      sstr.str("");
      sstr << "depth";
      runNorms(nx,ny,samplecount,depthcount, canvas, cam);
      sstr.str("");
      sstr << "normals";
      save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());
      sstr.str("");
      sstr << "normals";

      sstr.str("");
      runAlbedo(nx,ny,samplecount,depthcount, canvas, cam);
      sstr << "albedo";
      save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());
      sstr.str("");
      sstr << "albedo";

    }
    else{
      runPath(nx,ny, samplecount, depthcount, canvas, cam);
      std::stringstream sstr;
      sstr << "output";
      save(sstr, nx, ny, samplecount, canvas.GetColorBuffer());

    }
  }

  timer.Stop();
  elapsedTime1 = timer.GetElapsedTime();
  std::cout << " Elapsed time         = " << elapsedTime1 << std::endl;
}

