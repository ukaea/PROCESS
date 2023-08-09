# PROCESS Preprocessing Variables for CMake
# Author    :   K. Zarebski (UKAEA)
# Date      :   last modified 2020-11-05
#
# Retrieves all required information for setting the preprocessor
# variables during the PROCESS build

MACRO(FindPreprocessingVars)
    INCLUDE(${CMAKE_SOURCE_DIR}/cmake/utilities.cmake)

    EXECUTE_PROCESS(
        COMMAND bash -c "git -C ${CMAKE_SOURCE_DIR} show -s --format=format:%s"
        OUTPUT_VARIABLE COMMIT_MSG
    )
    STRING(STRIP ${COMMIT_MSG} COMMIT_MSG)
    EXECUTE_PROCESS(
        COMMAND bash -c "echo \"$(git -C ${CMAKE_SOURCE_DIR} diff | wc -l)\"|tr '\n' ' '"
        OUTPUT_VARIABLE GIT_DIFF
    )
    STRING(STRIP ${GIT_DIFF} GIT_DIFF)
    STRING(REPLACE "'" "|" COMMIT_MSG ${COMMIT_MSG})
    STRING(REPLACE "\"" "|" COMMIT_MSG ${COMMIT_MSG})
    STRING(REPLACE "#" "[hash]" COMMIT_MSG ${COMMIT_MSG})

    EXECUTE_PROCESS(
        COMMAND bash -c "echo \"$(git -C ${CMAKE_SOURCE_DIR} describe --tags)\"|tr '\n' ' '"
        OUTPUT_VARIABLE GIT_TAG
    )
    STRING(STRIP ${GIT_TAG} GIT_TAG)

    EXECUTE_PROCESS(
        COMMAND bash -c "echo \"$(git -C ${CMAKE_SOURCE_DIR} rev-parse --abbrev-ref HEAD)\"|tr '\n' ' '"
        OUTPUT_VARIABLE GIT_BRANCH
    )
    STRING(STRIP ${GIT_BRANCH} GIT_BRANCH)

    # CMake uses f95 (rather than gfortran) by default. However, this is usually
    # a symbolic link to gfortran. Set CMAKE_Fortran_COMPILER to the actual
    # location of gfortran for clarity
    EXECUTE_PROCESS(
        COMMAND bash -c "which gfortran | tr -d '[:space:]'"
        OUTPUT_VARIABLE CMAKE_Fortran_COMPILER
    )
    STRING(STRIP ${CMAKE_Fortran_COMPILER} CMAKE_Fortran_COMPILER)

    # The CMAKE_Fortran_COMPILER_VERSION can differ from "gfortran --version"
    # if gfortran has been updated, perhaps due to the f95 linking. Set
    # CMAKE_Fortran_COMPILER_VERSION to the actual gfortran version
    EXECUTE_PROCESS(
        COMMAND bash -c "gfortran -dumpversion"
        OUTPUT_VARIABLE CMAKE_Fortran_COMPILER_VERSION
    )
    STRING(STRIP ${CMAKE_Fortran_COMPILER_VERSION} CMAKE_Fortran_COMPILER_VERSION)

    # gfortran >= 9 required: differing regression test results with lower
    # versions. Reason unknown
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 9)
        MESSAGE(STATUS "Fortran Compiler: ${CMAKE_Fortran_COMPILER}")
        MESSAGE(STATUS "Fortran Compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
        MESSAGE(STATUS "Fortran Compiler Version: ${CMAKE_Fortran_COMPILER_VERSION}")
        MESSAGE(FATAL_ERROR "cmake detected gfortran version is " ${CMAKE_Fortran_COMPILER_VERSION} ". 9 or above required")
    ENDIF()

    FOREACH(VAR COMMIT_MSG GIT_DIFF GIT_TAG GIT_BRANCH CMAKE_Fortran_COMPILER CMAKE_Fortran_VERSION)
        IF(NOT VAR)
            MESSAGE(FATAL_ERROR "Failed to obtain value for '${VAR}'")
        ENDIF()
    ENDFOREACH()

    EXECUTE_PROCESS(
        COMMAND bash -c "${Python3_EXECUTABLE} -c \"import site, os;print(os.path.join(site.getsitepackages()[0], '${PROJECT_NAME}'))\""
        OUTPUT_VARIABLE PROCESS_MODULE_INSTALL_LOCATION
    )
    STRING(STRIP ${PROCESS_MODULE_INSTALL_LOCATION} PROCESS_MODULE_INSTALL_LOCATION)

    # ---------- Summarise Preprocessor Flags in Output ---------- #
    MESSAGE(STATUS "[Preprocessor Variables]: ")
    MESSAGE(STATUS "\tINSTALLDIR : ${CMAKE_SOURCE_DIR}")
    MESSAGE(STATUS "\tCOMMSG : ${COMMIT_MSG}")
    MESSAGE(STATUS "\tbranch_name : ${GIT_BRANCH}")
    MESSAGE(STATUS "\ttagno : ${GIT_TAG}")
    MESSAGE(STATUS "\tuntracked : ${GIT_DIFF}")
    MESSAGE(STATUS "\tFortran compiler version : ${CMAKE_Fortran_COMPILER_VERSION}")

    # ------------------------------------------------------------ #
    ADD_DEFINITIONS(-DINSTALLDIR="${CMAKE_SOURCE_DIR}")
    ADD_DEFINITIONS(-Dtagno="${GIT_TAG}")
    ADD_DEFINITIONS(-Dbranch_name="${GIT_BRANCH}")
    ADD_DEFINITIONS(-Duntracked=${GIT_DIFF})

    ensure_string_length(${COMMIT_MSG} 145 COMMIT_MSG) # 1502 line truncation error fix

    ADD_DEFINITIONS(-DCOMMSG="${COMMIT_MSG}")
    MESSAGE(STATUS "\tTruncated commit message length : ${COMMIT_MSG}")
ENDMACRO()
