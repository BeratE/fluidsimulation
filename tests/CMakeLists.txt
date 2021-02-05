cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)

### Project
set(TARGET_NAME tests)

set(LIBRARIES
  CompactNSearch
  learnSPH)

set(INCLUDE_DIRS
  PUBLIC
  ${PROJECT_SOURCE_DIR}/extern/Eigen
  ${PROJECT_SOURCE_DIR}/extern
  ${PROJECT_SOURCE_DIR}/common/inc)

### Sources
set(SOURCE_FILES
  util.h
  main.cpp  
  comparison.cpp
  mc_surface.cpp
  pbf_simulation.cpp
  sph_kernel.cpp
  sph_sampling.cpp
  sph_simulation.cpp
  tension_test.cpp
  )

### Target  
add_executable(${TARGET_NAME} ${SOURCE_FILES})
add_dependencies(${TARGET_NAME} ${LIBRARIES})
target_link_libraries(${TARGET_NAME} ${LIBRARIES})
target_include_directories(${TARGET_NAME} ${INCLUDE_DIRS})
