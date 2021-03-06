cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

project(learnSPH)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

configure_file (
  "${PROJECT_SOURCE_DIR}/cmake/config.h.in"
  "${PROJECT_SOURCE_DIR}/common/inc/config.h")

option(USE_CUDA "Use CUDA" OFF)


### CUDA
if (USE_CUDA)
  enable_language(CUDA)
  add_definitions(-DUSE_CUDA)
  
  if(NOT DEFINED ${CMAKE_CUDA_ARCHITECTURES})
    set(CMAKE_CUDA_ARCHITECTURES 50)
  endif()

  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)

  set(CMAKE_CUDA_FLAGS "-Wall -Wno-deprecated-gpu-targets")
  set(CMAKE_CUDA_FLAGS_DEBUG "-g")

  set(LIBRARIES ${CUDA_LIBRARIES})

  set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
  set(CUDA_RESOLVE_DEVICE_SYMBOLSON ON)  
endif()

if (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # Parallel compilation in Visual Studio
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else () 
   # Disable warnings for all the unknown OpenMP pragmas. 
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas") 
endif () 


# Dependencies
add_subdirectory(extern/CompactNSearch)

# Common/shared source files
file(GLOB COMMON_SOURCES ${PROJECT_SOURCE_DIR}/common/src/*.cpp)

# Modules
add_subdirectory(learnSPH)
add_subdirectory(surface)
add_subdirectory(app)
add_subdirectory(tests)

# Copy res folder
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/res DESTINATION ${CMAKE_BINARY_DIR})
