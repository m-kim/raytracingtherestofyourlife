cmake_minimum_required(VERSION 3.8)

project(PathTracing VERSION 0.1 LANGUAGES CXX)

find_package(VTKm REQUIRED
  OPTIONAL_COMPONENTS CUDA OpenMP)

set (PathTracingSrc main.cc CornellBox.cpp Camera.cxx PathTracer.cxx
  SphereExtractor.cxx QuadExtractor.cxx QuadIntersector.cxx
  SphereIntersector.cxx QuadPdf.cxx SpherePdf.cxx QuadGenerateDir.cxx SphereGenerateDir.cxx
  CosineGenerateDir.cxx WhichGenerateDir.cxx MapperPathTracer.cxx)
set (PathTracingHdr Worklets.h AABBSurface.h BVHTraverser.h
  EmitWorklet.h ScatterWorklet.h SurfaceWorklets.h Surface.h PdfWorklet.h Record.h onb.h
  CornellBox.h vec3.h wangXor.h GenerateDir.h)

add_executable(CornellBox ${PathTracingSrc} ${PathTracingHdr})


set_target_properties(CornellBox PROPERTIES
  CXX_STANDARD 17)


target_link_libraries(CornellBox PUBLIC vtkm_cont vtkm_rendering)


if(VTKm_ENABLE_CUDA)
  vtkm_compile_as_cuda(cudaSource ${PathTracingSrc})
  add_executable(CornellBox_CUDA ${cudaSource} )
  target_compile_definitions(CornellBox_CUDA PRIVATE "VTKM_DEVICE_ADAPTER=VTKM_DEVICE_ADAPTER_CUDA")

  set_target_properties(CornellBox_CUDA PROPERTIES
    CUDA_STANDARD 14)
  target_link_libraries(CornellBox_CUDA PUBLIC vtkm_cont vtkm_rendering)

endif()
