cmake_minimum_required(VERSION 3.0.0)
project(QCD_Analysis VERSION 0.1.0)

include(CTest)
enable_testing()

find_package(Eigen3 REQUIRED)
MESSAGE( [Main] " EIGEN3_INCLUDE_DIRS = ${EIGEN3_INCLUDE_DIRS}") 
# EIGEN3_INCLUDE_DIRS = /usr/include/eigen3

include_directories(${EIGEN3_INCLUDE_DIRS})

# SET(CMAKE_CXX_FLAGS "-lgsl -lgslcblas -lm -O3")

find_package(GSL REQUIRED)
link_libraries(GSL::gsl GSL::gslcblas)

if (MSVC)
  string(REGEX REPLACE "/W[0-9]" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
endif (MSVC)

# list(APPEND CMAKE_PREFIX_PATH /home/pesznyak/GSL)

# 2D function fit for two correlated sets
add_executable(Fit Fit.cpp)
add_executable(Extrapolate Extrapolate.cpp)
add_executable(ContinuumLimes1D ContinuumLimes1D.cpp)
add_executable(ContinuumLimes2D ContinuumLimes2D.cpp)
add_executable(ContinuumLimes2DEx ContinuumLimes2DEx.cpp)
add_executable(FitExtra FitExtra.cpp)
add_executable(HRGCalculations HRGCalculations.cpp)

target_link_libraries(ContinuumLimes2D PRIVATE Eigen3::Eigen)

set_target_properties(Fit PROPERTIES CXX_STANDARD 17
                                                 CXX_STANDARD_REQUIRED ON
                                                 CXX_EXTENSIONS OFF)

target_compile_options(Fit PRIVATE $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-Wall -Wextra -pedantic>
                                               $<$<CXX_COMPILER_ID:MSVC>:/W4 /permissive->)
                                               
#if(BUILD_TESTING)
#  add_executable(test1 test1.cpp)
#  add_test(NAME test_op COMMAND test_op)
#endif(BUILD_TESTING)
                                               
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)