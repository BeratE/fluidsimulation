cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

### Project
set(TARGET_NAME learnSPH_app)

set(LIBRARIES
  learnSPH)

set(INCLUDE_DIRS
  PUBLIC
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/common/inc)

### External dependencies

### Sources
set(SOURCE_FILES
    main.cpp)

### Target  
add_executable(${TARGET_NAME} ${SOURCE_FILES})
add_dependencies(${TARGET_NAME} learnSPH)
target_link_libraries(${TARGET_NAME} ${LIBRARIES})
target_include_directories(${TARGET_NAME} ${INCLUDE_DIRS})
