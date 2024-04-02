cmake_minimum_required(VERSION 3.25)
project(DIProcess)
set(INC_DIR E:\\Arib\\handout\\code\\rsdiplib\\include)
set(LINK_DIR E:\\Arib\\handout\\code\\rsdiplib\\lib)
include_directories(${INC_DIR})
link_directories(${LINK_DIR})
link_libraries(RSDIPLib)

set(CMAKE_CXX_STANDARD 98)
set(CMAKE_MFC_FLAG 2)

add_executable(Ref ImageProcessingEx.cpp ImageProcessingEx.h)
add_executable(DIProcess CImageProcessingEx.cpp CImageProcessingEx.h CalFunction.cpp CalFunction.h PCAMatrix.cpp PCAMatrix.h)
add_executable(TEST Test.cpp)

target_link_libraries(DIProcess RSDIPLib)
target_link_libraries(TEST RSDIPLib)

