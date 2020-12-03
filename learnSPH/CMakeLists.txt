# External dependencies

set(LIBRARIES 
  CompactNSearch)

set(INCLUDE_DIRS
  PUBLIC
  ${PROJECT_SOURCE_DIR}/extern)

## Eigen
set(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/extern/Eigen)
if(NOT EIGEN3_INCLUDE_DIR)
  message(FATAL_ERROR
    "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
endif()
list(APPEND INCLUDE_DIRS "${EIGEN3_INCLUDE_DIR}")

## OpenCL
find_package(OpenCL)
if (OpenCL_FOUND)
  list(APPEND LIBRARIES "${OpenCL_LIBRARIES}")
  list(APPEND INCLUDE_DIRS "${OpenCL_INCLUDE_DIRS}")
endif()

set(SOURCE_FILES
	util/obj_reader.cpp
	util/obj_reader.h
	util/vtk_writer.cpp
	util/vtk_writer.h
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
        kernel.h
	kernel.cpp)
    
# Create a library version so it can be linked with other modules
add_library(learnSPH ${SOURCE_FILES})
target_link_libraries(learnSPH ${LIBRARIES})
target_include_directories(learnSPH ${INCLUDE_DIRS} )
