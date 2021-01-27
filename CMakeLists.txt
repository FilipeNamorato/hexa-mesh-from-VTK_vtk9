cmake_minimum_required(VERSION 2.6)
 
project(HexaMeshFromVTK)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(HexaMeshFromVTK src/main.cpp)


if(VTK_LIBRARIES)
  target_link_libraries(HexaMeshFromVTK ${VTK_LIBRARIES})
else()
  target_link_libraries(HexaMeshFromVTK vtkHybrid)
endif()
