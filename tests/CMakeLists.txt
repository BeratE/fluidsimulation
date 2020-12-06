cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

### Project
set(TARGET_NAME tests)

set(LIBRARIES
  learnSPH)

set(INCLUDE_DIRS
  PUBLIC
  ${PROJECT_SOURCE_DIR})

### Sources
set(SOURCE_FILES
    main.cpp
    sph_kernel.cpp
    sph_sampling.cpp
    sph_runtime.cpp
    sph_simulation.cpp)

### Target  
add_executable(${TARGET_NAME} ${SOURCE_FILES})
add_dependencies(${TARGET_NAME} learnSPH)
target_link_libraries(${TARGET_NAME} ${LIBRARIES})
target_include_directories(${TARGET_NAME} ${INCLUDE_DIRS})
