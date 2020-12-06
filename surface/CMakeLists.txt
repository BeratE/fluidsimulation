cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)

### Project
set(TARGET_NAME surface)

set(LIBRARIES "")

set(INCLUDE_DIRS
  PUBLIC
  ${PROJECT_SOURCE_DIR}/extern
  ${PROJECT_SOURCE_DIR}/common
  ${PROJECT_SOURCE_DIR}/extern/Eigen)

### Sources
set(SOURCE_FILES
  mc_lut.hpp
  surface.h
  surface.cpp)

### Target Library
add_library(${TARGET_NAME} SHARED ${SOURCE_FILES})
target_link_libraries(${TARGET_NAME} ${LIBRARIES})
target_include_directories(${TARGET_NAME} ${INCLUDE_DIRS})
