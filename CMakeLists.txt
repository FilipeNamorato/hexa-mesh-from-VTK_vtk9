cmake_minimum_required(VERSION 3.12)  # Recomendo pelo menos 3.12 para bom suporte ao VTK 9

project(HexaMeshFromVTK)

find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    FiltersCore
    FiltersGeneral
    FiltersPoints
    IOXML
    IOXMLParser
    FiltersSources
    FiltersGeometry
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# compiling from python environment
set(CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX})

add_executable(HexaMeshFromVTK src/main.cpp src/config.cpp)

target_link_libraries(HexaMeshFromVTK PRIVATE ${VTK_LIBRARIES})

target_include_directories(HexaMeshFromVTK PRIVATE src)