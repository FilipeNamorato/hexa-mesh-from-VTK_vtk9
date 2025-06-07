cmake_minimum_required(VERSION 3.12)  # Recomendo pelo menos 3.12 para bom suporte ao VTK 9

project(HexaMeshFromVTK)

# Especifica explicitamente os módulos que você está usando (com base no seu main.cpp)
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

# Definir saída dos binários
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Criar o executável
add_executable(HexaMeshFromVTK src/main.cpp src/config.cpp)

# Linkar com VTK
target_link_libraries(HexaMeshFromVTK PRIVATE ${VTK_LIBRARIES})

target_include_directories(HexaMeshFromVTK PRIVATE src)