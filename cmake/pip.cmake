# Check all Pip modules exist
# Author  :   K. Zarebski
# Date    : last modified 2020-11-05
#
# Checks to see if all Python module requirements are satisfied else
# runs python pip on the requirements file

MACRO(PIP_INSTALL)
    SET(PIP_NAME "pip_installs")

    EXECUTE_PROCESS(
        COMMAND bash -c "${Python3_EXECUTABLE} -c \"import ford\""
        OUTPUT_VARIABLE FORD_CHECK
        ERROR_QUIET
    )

    IF(NOT FORD_CHECK)
        MESSAGE(STATUS "\tFORD Install: Python module 'ford' not found, will install.")
        include(ExternalProject)
        ExternalProject_Add(
            ford_git
            GIT_REPOSITORY https://github.com/jonmaddock/ford.git
            UPDATE_COMMAND ""
            INSTALL_COMMAND bash -c "${Python3_EXECUTABLE} -m pip install ."
            BUILD_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_IN_SOURCE 1
            BYPRODUCTS ${CMAKE_BINARY_DIR}/ford_git-prefix/
        )
    ELSE()
        ADD_CUSTOM_TARGET(
            ford_git
        )
        MESSAGE(STATUS "\tFORD Install: Python module 'ford' found, skipping install.")
    ENDIF()

    SET(MODULE_REQUIREMENTS_FILE ${CMAKE_SOURCE_DIR}/requirements.txt)
    STRING(REPLACE "/" "_" PIP_OUT_PREFIX ${Python3_EXECUTABLE})
    SET(PIP_COMPLETE_FILE ${CMAKE_BINARY_DIR}/${PIP_OUT_PREFIX}.touch)
    ADD_CUSTOM_TARGET(
        ${PIP_NAME}
        DEPENDS ${PIP_COMPLETE_FILE}
    )

    # Manually install numpy as including it in requirements doesnt install it
    # It is a pre-requisite to f90wrap install
    ADD_CUSTOM_COMMAND(
        OUTPUT ${PIP_COMPLETE_FILE}
        COMMAND ${Python3_EXECUTABLE} -m pip install numpy
        COMMAND ${Python3_EXECUTABLE} -m pip install -r ${MODULE_REQUIREMENTS_FILE}
        COMMAND touch ${PIP_COMPLETE_FILE}
    )

    ADD_DEPENDENCIES(${PIP_NAME} ford_git)
ENDMACRO()
