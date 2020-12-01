cmake_minimum_required(VERSION 3.1)
project(learnSPH)
set(CMAKE_CXX_STANDARD 17)

configure_file (
    "${PROJECT_SOURCE_DIR}/cmake/config.h.in"
    "${PROJECT_SOURCE_DIR}/learnSPH/util/config.h"
)

if (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # Parallel compilation in Visual Studio
endif()

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Dependencies
add_subdirectory(extern/CompactNSearch)

# Modules
add_subdirectory(learnSPH)
add_subdirectory(tests)
add_subdirectory(app)

# Copy res folder
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/res DESTINATION ${CMAKE_BINARY_DIR})
