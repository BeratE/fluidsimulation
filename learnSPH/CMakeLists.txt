cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

### Project
set(TARGET_NAME learnSPH)

set(LIBRARIES 
  CompactNSearch
  surface)

set(INCLUDE_DIRS
  PUBLIC
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/extern
  ${PROJECT_SOURCE_DIR}/common/inc)


### External dependencies
## Eigen
set(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/extern/Eigen)
if(NOT EIGEN3_INCLUDE_DIR)
  message(FATAL_ERROR
    "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
endif()
list(APPEND INCLUDE_DIRS "${EIGEN3_INCLUDE_DIR}")

### Sources
set(SOURCE_FILES
  ${COMMON_SOURCES}
  system/particlesystem.h
  system/particlesystem.cpp
  system/fluidsystem.h
  system/fluidsystem.cpp
  system/boundarysystem.h        
  system/boundarysystem.cpp
  system/emitter.h
  system/emitter.cpp
  solver.h
  solver.cpp
  solver_sph.h
  solver_sph.cpp
  kernel.h
  kernel.cpp)

### Target Library
add_library(${TARGET_NAME} ${SOURCE_FILES})
target_link_libraries(${TARGET_NAME} ${LIBRARIES})
target_include_directories(${TARGET_NAME} ${INCLUDE_DIRS})
