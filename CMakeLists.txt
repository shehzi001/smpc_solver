# General parameters
cmake_minimum_required (VERSION 2.6)
project (smpc_solver CXX C)

include (CheckFunctionExists)
include (CheckIncludeFile)


####################################
# Options
####################################
option (BUILD_TESTS         "Build tests" ON)


####################################
# Settings
####################################
set (smpc_solver_SOURCE_DIR "${PROJECT_SOURCE_DIR}/solver/")
set (wmg_SOURCE_DIR "${PROJECT_SOURCE_DIR}/WMG/")
set (CMAKE_SKIP_ASSEMBLY_SOURCE_RULES ON)
set (CMAKE_SKIP_PREPROCESSED_SOURCE_RULES ON)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/test/")
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib/")
set (CMAKE_VERBOSE_MAKEFILE ON)
set (CMAKE_RULE_MESSAGES OFF CACHE BOOL "")



####################################
# Other parameters
####################################


find_path(EIGEN_INCLUDE_DIRS "Eigen/Core" PATHS "${SDK_DIR}/include" PATH_SUFFIXES "eigen3")
if (EIGEN_INCLUDE_DIRS)
    include_directories ("${EIGEN_INCLUDE_DIRS}")
endif (EIGEN_INCLUDE_DIRS)
include_directories ("${PROJECT_SOURCE_DIR}/include/")


set (CMAKE_REQUIRED_LIBRARIES "m")
check_function_exists (feenableexcept HAVE_FEENABLEEXCEPT)
configure_file ("${smpc_solver_SOURCE_DIR}/solver_config.h.in" "${smpc_solver_SOURCE_DIR}/solver_config.h" )


file (GLOB SMPC_SRC "${smpc_solver_SOURCE_DIR}/*.cpp")
add_library (smpc_solver STATIC ${SMPC_SRC})

file (GLOB WMG_SRC "${wmg_SOURCE_DIR}/*.cpp")
add_library (wmg STATIC ${WMG_SRC})


if (BUILD_TESTS)
    set (test_DIR "${PROJECT_SOURCE_DIR}/test/")
    file (GLOB TESTS RELATIVE "${test_DIR}" "${test_DIR}test_*.cpp")
    foreach (testname ${TESTS})
        string(REPLACE ".cpp" ".a" targetname "${testname}")
        add_executable (${targetname} "${test_DIR}/${testname}")
        target_link_libraries (${targetname} smpc_solver wmg)
    endforeach (testname ${TESTS})
endif (BUILD_TESTS)
