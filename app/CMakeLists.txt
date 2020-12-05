cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Local
set(SOURCE_FILES
    main.cpp)
    
add_executable(learnSPH_app ${SOURCE_FILES})
add_dependencies(learnSPH_app learnSPH)
target_link_libraries(learnSPH_app learnSPH)
target_include_directories(learnSPH_app PUBLIC ${PROJECT_SOURCE_DIR})
