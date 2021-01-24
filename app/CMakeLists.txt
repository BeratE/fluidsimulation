cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

### Project
set(TARGET_NAME learnSPH_app)



set(INCLUDE_DIRS
  PUBLIC
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/extern
  ${PROJECT_SOURCE_DIR}/extern/Eigen
  ${PROJECT_SOURCE_DIR}/common/inc)

### External dependencies

### Sources
set(SOURCE_FILES
    main.cpp)

### Target  
add_executable(${TARGET_NAME} ${SOURCE_FILES})
target_link_libraries(${TARGET_NAME} ${LIBRARIES})
target_include_directories(${TARGET_NAME} ${INCLUDE_DIRS})
