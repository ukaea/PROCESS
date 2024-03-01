# CMake file for Pythonic PROCESS
# Author  :   PROCESS Team (UKAEA)
# Date    :   last modified 2023-07-24

# Specify the minimum version for CMake
# 3.12 is required due to use of list TRANSFORM commands
# 3.15 for cmake policy CMP0094
CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

# Set project name
PROJECT(process LANGUAGES Fortran)

cmake_policy(SET CMP0094 NEW)

# Ensure python3 interpreter is used
if(CMAKE_HOST_APPLE)
    SET(CMAKE_FIND_FRAMEWORK NEVER)
    set(CMAKE_MACOSX_RPATH ON)
endif()

# Python 3.8 or greater
find_package(Python3 3.8 COMPONENTS Interpreter Development)

# Read in external CMake scripts for preprocessing, f2py, ford, etc.
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/pip.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/process_preprocessing.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/preprocess.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/f2py.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/ford.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/gfortranlibs.cmake)

# Retrieve and set information for preprocessor
FindPreprocessingVars()

SET(CMAKE_Fortran_FLAGS "-cpp -ffree-line-length-none -fPIC -ftest-coverage -fprofile-arcs -finit-local-zero")

IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # space before -g is very necessary
    STRING(CONCAT CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} " -g -O0 -fcheck=all,no-bounds,no-array-temps")
ENDIF()

# Get System information, log as much information as we can
# for ease later on when debugging system specific issues
CMAKE_HOST_SYSTEM_INFORMATION(RESULT OS_LABEL QUERY OS_NAME)
CMAKE_HOST_SYSTEM_INFORMATION(RESULT OS_VER QUERY OS_RELEASE)
CMAKE_HOST_SYSTEM_INFORMATION(RESULT OS_PLAT QUERY OS_PLATFORM)
MESSAGE(STATUS "[System Information]: ")
MESSAGE(STATUS "\tOperating System: ${OS_LABEL} ")
MESSAGE(STATUS "\tVersion: ${OS_VER}")
MESSAGE(STATUS "\tPlatform: ${OS_PLAT}")
MESSAGE(STATUS "[PROCESS Build Information]: ")
MESSAGE(STATUS "\tPython Binary: ${Python3_EXECUTABLE}")
MESSAGE(STATUS "\tBinary Location: ${CMAKE_BINARY_DIR}")
MESSAGE(STATUS "\tFortran Compiler: ${CMAKE_Fortran_COMPILER}")
MESSAGE(STATUS "\tFortran Compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
MESSAGE(STATUS "\tFortran Compiler Version: ${CMAKE_Fortran_COMPILER_VERSION}")
MESSAGE(STATUS "\tPROCESS Module Installation Directory: ${PROCESS_MODULE_INSTALL_LOCATION}")
MESSAGE(STATUS "\tBuild type: ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "\tFortran arguments: ${CMAKE_Fortran_FLAGS}")

# Define various addresses for locations of files
# relative to the CMake root directory
SET(PROCESS_SRC_DIR ${CMAKE_SOURCE_DIR}/source/fortran)
SET(PROCESS_PYTHON_SRC_DIR ${CMAKE_SOURCE_DIR}/process)

SET(PYTHON_MODULE_DIR ${CMAKE_SOURCE_DIR}/process)
SET(PYTHON_LIBS_DIR ${PYTHON_MODULE_DIR}/lib)
SET(PYTHON_SOURCE_IO_DIR ${PYTHON_MODULE_DIR}/io)

FILE(GLOB PROCESS_SRC_PATHS ${PROCESS_SRC_DIR}/*.f90 ${PROCESS_SRC_DIR}/*.f)
FILE(GLOB PROCESS_PYTHON_SRC_PATHS ${PROCESS_PYTHON_SRC_DIR}/*.py)

# Define interface source filenames to wrap
LIST(APPEND PROCESS_SRCS
    numerics.f90
    scan.f90
    fw.f90
    hcpb.f90
    pfcoil.f90
    reinke_module.f90
    sctfcoil.f90
    plasma_profiles.f90
    final_module.f90
    cost_variables.f90
    divertor_variables.f90
    fwbs_variables.f90
    physics.f90
    physics_functions.f90
    physics_variables.f90
    tfcoil_variables.f90
    times_variables.f90
    ife_variables.f90
    ife.f90
    heat_transport_variables.f90
    buildings_variables.f90
    maths_library.f90
    iteration_variables.f90
    evaluators.f90
    water_usage_variables.f90
    constraint_equations.f90
    constants.f90
    build_variables.f90
    current_drive_variables.f90
    pfcoil_variables.f90
    structure_variables.f90
    output.f90
    init_module.f90
    main_module.f90
    error_handling.f90
    global_variables.f90
    constraint_variables.f90
    vacuum_variables.f90
    CS_fatigue_variables.f90
    impurity_radiation.f90
    pulse_variables.f90
    pf_power_variables.f90
    dcll.f90
    rebco_variables.f90
    primary_pumping_variables.f90
    reinke_variables.f90
    neoclassics_module.f90
    blanket_library.f90
    stellarator_variables.f90
    stellarator.f90
    stellarator_configuration.f90
    stellarator_fwbs.f90
    startup_variables.f90
)

PREPROCESS()

# Install all requirements in the current Python environment
PIP_INSTALL()

# Build targets
# Stage 1: Compile Fortran sources to a shared library
# (libprocess.so on Linux)
ADD_LIBRARY(${PROJECT_NAME} SHARED ${PROCESS_SRC_PATHS})
SET_TARGET_PROPERTIES(${PROJECT_NAME}
    PROPERTIES LINKER_LANGUAGE Fortran
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/process/lib
)

IF(NOT CMAKE_HOST_APPLE)
    SET_TARGET_PROPERTIES(${PROJECT_NAME}
        PROPERTIES BUILD_RPATH $ORIGIN
    )
ENDIF()

ADD_DEPENDENCIES(${PROJECT_NAME} ${PIP_NAME})

IF(CMAKE_HOST_APPLE)
    SET(LIBRARY_OUTPUT_SUFFIX ".dylib")
ELSEIF(CMAKE_HOST_UNIX)
    SET(LIBRARY_OUTPUT_SUFFIX ".so")
ENDIF()

# Stage 3: Run f2py
# Builds every time to make sure the libprocess.so file is up-to-date
F2PY()

# Run FORD on a "fast" md file, which doesn't make graphs or a search index
# These are wanted for making docs, but not when making dicts quickly
# This produces the pickled project file
FORD()
DICTS()

# Copy libgfortran from system into Python module folder
GET_GFORTRANLIBS()

SET(PROCESS_MODULE_DIST ${CMAKE_SOURCE_DIR}/process_dist)

# Delete any existing .gcda files
ADD_CUSTOM_TARGET(coverage_cleanup
    COMMAND echo "-- Deleting any existing .gcda coverage files..."
    COMMAND find ${CMAKE_SOURCE_DIR}/build/ -name '*.gcda' -exec rm -rf {} \\\;
)

# Run the Pip Install
ADD_CUSTOM_TARGET(install_process ALL
    COMMAND ${Python3_EXECUTABLE} -m pip install --no-cache-dir -e ${CMAKE_SOURCE_DIR}
)

ADD_DEPENDENCIES(install_process ${PIP_NAME} f2py ${FORD_NAME} ${GFORTLIB_NAME} ${DICTS_NAME} coverage_cleanup)

if(EXISTS ${CMAKE_SOURCE_DIR}/lcov_results)
    execute_process(COMMAND echo "-- Delete existing lcov_coverage directory...")
    execute_process(COMMAND cmake -E remove_directory ${CMAKE_SOURCE_DIR}/lcov_results)
    execute_process(COMMAND mkdir ${CMAKE_SOURCE_DIR}/lcov_results)
elseif(NOT EXISTS ${CMAKE_SOURCE_DIR}/lcov_results)
    execute_process(COMMAND echo "-- Create lcov_coverage directory...")
    execute_process(COMMAND mkdir ${CMAKE_SOURCE_DIR}/lcov_results)
endif()

# Custom target to run the regression tests and generate coverage report
ADD_CUSTOM_TARGET(regression
    COMMAND pytest -v tests/regression/ -s | tee tests/regression/pytest.log
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

ADD_CUSTOM_TARGET(coverage
    COMMAND lcov --gcov-tool gcov --capture --directory build/CMakeFiles/process.dir/source/fortran/ --output-file lcov_results/coverage.info
    COMMAND genhtml --output-directory lcov_results/html lcov_results/coverage.info
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
