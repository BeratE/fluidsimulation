cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)

### Project
set(TARGET_NAME tension)

set(LIBRARIES
	CompartNSearch
	learnSPH)

set(INCLUDE_DIRS
  PUBLIC
  ${PROJECT_SOURCE_DIR}/extern/Eigen
  ${PROJECT_SOURCE_DIR}/extern
  ${PROJECT_SOURCE_DIR}/common
  )

### Sources
set(SOURCE_FILES
  tension.h
  tension.cpp)



### Target Library
add_library(${TARGET_NAME} ${SOURCE_FILES})
add_dependencies(${TARGET_NAME} learnSPH)
target_link_libraries(${TARGET_NAME} ${LIBRARIES})
target_include_directories(${TARGET_NAME} ${INCLUDE_DIRS})
