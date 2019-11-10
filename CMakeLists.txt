cmake_minimum_required(VERSION 3.8)

project(PathTracing VERSION 0.1 LANGUAGES CXX)

find_package(VTKm REQUIRED
  OPTIONAL_COMPONENTS CUDA OpenMP)

find_package(adios2 REQUIRED)
set(mpilib) #intentionally set mpilib blank
#VTKm requires OpenMP because pthreads is required in ColorTable(!!)
  find_package(OpenMP REQUIRED)
  find_package(MPI REQUIRED)

get_target_property(MPI_pthread MPI::MPI_CXX INTERFACE_COMPILE_OPTIONS)
list(REMOVE_ITEM MPI_pthread "-pthread")
get_target_property(OMP_pthread OpenMP::OpenMP_CXX INTERFACE_COMPILE_OPTIONS)
list(REMOVE_ITEM OMP_pthread "-pthread")

set (PathTracingSrc 
main.cc 
CornellBox.cpp 
pathtracing/Camera.cxx 
pathtracing/PathTracer.cxx
pathtracing/SphereExtractor.cxx 
pathtracing/QuadIntersector.cxx
pathtracing/SphereIntersector.cxx 
pathtracing/WhichGenerateDir.cxx 
MapperPathTracer.cxx 
MapperQuad.cxx
MapperQuadAlbedo.cxx
raytracing/ChannelBuffer.cxx
raytracing/RayTracerNormals.cxx
raytracing/RayTracerAlbedo.cxx
MapperQuadNormals.cxx
View3D.cxx )
set (PathTracingHdr
pathtracing/Camera.h

pathtracing/CosineGenerateDir.h
pathtracing/PathAlgorithms.h
pathtracing/PathTracer.h
pathtracing/Pdf.h
pathtracing/QuadGenerateDir.h
pathtracing/SphereExtractor.h
pathtracing/SphereIntersector.h
pathtracing/SphereGenerateDir.h
pathtracing/SpherePdf.h
pathtracing/WhichGenerateDir.h
pathtracing/Worklets.h 
pathtracing/QuadIntersector.h
pathtracing/AABBSurface.h 
pathtracing/BVHTraverser.h
pathtracing/EmitWorklet.h 
pathtracing/ScatterWorklet.h 
pathtracing/SurfaceWorklets.h 
pathtracing/Surface.h 
pathtracing/PdfWorklet.h 
pathtracing/Record.h 
pathtracing/onb.h
MapperPathTracer.h
MapperQuad.h
raytracing/ChannelBuffer.h
raytracing/ChannelBufferOperations.h
raytracing/Ray.h
raytracing/RayTracerNormals.h
MapperQuadNormals.h
CornellBox.h 
pathtracing/vec3.h 
pathtracing/wangXor.h 
pathtracing/GenerateDir.h 
pathtracing/QuadPdf.h 
pathtracing/SpherePdf.h
View3D.h
pathtracing/RayOperations.h)

add_executable(CornellBox ${PathTracingSrc} ${PathTracingHdr})


set_target_properties(CornellBox PROPERTIES
  CXX_STANDARD 17)
set_target_properties(CornellBox PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_compile_features(CornellBox PRIVATE cxx_std_14)

target_include_directories(CornellBox
  PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  )
target_link_libraries(CornellBox PUBLIC vtkm_cont vtkm_rendering OpenMP::OpenMP_CXX MPI::MPI_CXX adios2::adios2 ${mpilib} )


if(VTKm_ENABLE_CUDA)
  set_source_files_properties(${PathTracingSrc} PROPERTIES LANGUAGE "CUDA")
  set_target_properties(CornellBox PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  target_compile_options(CornellBox PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=compute_50,code=sm_50>)

  target_include_directories(CornellBox
    PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
    )
  set_target_properties(CornellBox PROPERTIES
    CUDA_STANDARD 14)
  target_link_libraries(CornellBox PUBLIC vtkm::cuda)
endif()
