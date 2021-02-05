cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)

### Project
set(TARGET_NAME surface)

set(LIBRARIES 
    CompactNSearch
    learnSPH)

set(INCLUDE_DIRS
  PUBLIC
  ${PROJECT_SOURCE_DIR}/extern
  ${PROJECT_SOURCE_DIR}/common
  ${PROJECT_SOURCE_DIR}/extern/Eigen)

### Sources
set(SOURCE_FILES
  ${COMMON_SOURCES}
  main.cpp
  mc_lut.h
  surface.h
  surface.cpp)

if(USE_CUDA)
  list(APPEND SOURCE_FILES
    marching_cubes.cu)
else()
  list(APPEND SOURCE_FILES
    marching_cubes.cpp)
endif()


### Target Library
add_executable(${TARGET_NAME} ${SOURCE_FILES})
add_dependencies(${TARGET_NAME} ${LIBRARIES})
target_link_libraries(${TARGET_NAME} ${LIBRARIES})
target_include_directories(${TARGET_NAME} ${INCLUDE_DIRS})
