cmake_minimum_required(VERSION 3.8)

project(PathTracing VERSION 0.1 LANGUAGES CXX)

find_package(VTKm REQUIRED
  OPTIONAL_COMPONENTS CUDA OpenMP)

set (PathTracingSrc main.cc CornellBox.cpp Camera.cxx )
set (PathTracingHdr Worklets.h AABBSurface.h BVHTraverser.h
  EmitWorklet.h ScatterWorklet.h SurfaceWorklets.h Surface.h PdfWorklet.h Record.h onb.h
  CornellBox.h vec3.h wangXor.h)

add_executable(CornellBox ${PathTracingSrc} ${PathTracingHdr})

target_include_directories(CornellBox PUBLIC /home/mark/local/anaconda3/envs/cling/include/)
target_include_directories(CornellBox PUBLIC /home/mark/local/vtk-m-1.3.0/include/vtkm-1.3/)

set_target_properties(CornellBox PROPERTIES
  CXX_STANDARD 17)


target_link_libraries(CornellBox PUBLIC
  "/home/mark/local/vtk-m-1.3.0/lib/libvtkm_rendering-1.3.so"
  "/home/mark/local/vtk-m-1.3.0/lib/libvtkm_cont-1.3.so"
  "/home/mark/local/anaconda3/envs/cling/lib/libgomp.so"
  "/home/mark/local/anaconda3/envs/cling/lib/liblodepng.so")


if(VTKm_ENABLE_CUDA)
  vtkm_compile_as_cuda(cudaSource ${PathTracingSrc})
  add_executable(CornellBox_CUDA ${cudaSource} )
  target_compile_definitions(CornellBox_CUDA PRIVATE "VTKM_DEVICE_ADAPTER=VTKM_DEVICE_ADAPTER_CUDA")

  target_include_directories(CornellBox_CUDA PUBLIC /home/mark/local/anaconda3/envs/cling/include/)
  target_include_directories(CornellBox_CUDA PUBLIC /home/mark/local/vtk-m-1.3.0/include/vtkm-1.3/)
  set_target_properties(CornellBox_CUDA PROPERTIES
    CUDA_STANDARD 14)
  target_link_libraries(CornellBox_CUDA PUBLIC
    "/home/mark/local/vtk-m-1.3.0/lib/libvtkm_rendering-1.3.so"
    "/home/mark/local/vtk-m-1.3.0/lib/libvtkm_cont-1.3.so"
    "/home/mark/local/anaconda3/envs/cling/lib/libgomp.so"
    "/home/mark/local/anaconda3/envs/cling/lib/liblodepng.so")

endif()
