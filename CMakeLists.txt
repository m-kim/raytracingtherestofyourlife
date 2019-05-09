cmake_minimum_required(VERSION 3.8)

project(PathTracing VERSION 0.1 LANGUAGES CXX)

add_executable(CornellBox main.cc aabb.h      camera.h           material.h       ray.h              texture.h
aarect.h    constant_medium.h  moving_sphere.h  sphere.h           vec3.h
box.h       foomaterial.h      onb.h            stb_image.h
bucamera.h  hitable.h          pdf.h            stb_image_write.h
bvh.h       hitable_list.h     perlin.h         surface_texture.h Worklets.h
EmitWorklet.h ScatterWorklet.h)

target_include_directories(CornellBox PUBLIC /home/mark/local/anaconda3/envs/cling/include/)
target_include_directories(CornellBox PUBLIC /home/mark/local/mark-vtk-m/include/vtkm-1.0/)

set_target_properties(CornellBox PROPERTIES
  CXX_STANDARD 14)


target_link_libraries(CornellBox PUBLIC
  "/home/mark/local/mark-vtk-m/lib/libvtkm_rendering-1.0.so.1"
  "/home/mark/local/mark-vtk-m/lib/libvtkm_cont-1.0.so.1"
  "/home/mark/local/anaconda3/envs/cling/lib/libgomp.so"
  "/home/mark/local/anaconda3/envs/cling/lib/liblodepng.so")
