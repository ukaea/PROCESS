# Check all Pip modules exist
# Author  :   PROCESS Team (UKAEA)
# Date    : last modified 2023-07-24
#
# Checks to see if all Python module requirements are satisfied else
# runs python pip on the requirements file

MACRO(PIP_INSTALL)
    SET(PIP_NAME "pip_installs")

    if(RELEASE)
        SET(RELEASE TRUE)
    else()
        SET(RELEASE FALSE)
    endif()

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
    SET(DEVELOP_REQUIREMENTS_FILE ${CMAKE_SOURCE_DIR}/requirements_dev.txt)
    STRING(REPLACE "/" "_" PIP_OUT_PREFIX ${Python3_EXECUTABLE})
    SET(PIP_COMPLETE_FILE ${CMAKE_BINARY_DIR}/${PIP_OUT_PREFIX}.touch)
    ADD_CUSTOM_TARGET(
        ${PIP_NAME}
        DEPENDS ${PIP_COMPLETE_FILE}
    )

    ADD_CUSTOM_COMMAND(
        OUTPUT ${PIP_COMPLETE_FILE}
        COMMAND ${Python3_EXECUTABLE} -m pip install -r ${MODULE_REQUIREMENTS_FILE}
        COMMAND touch ${PIP_COMPLETE_FILE}
    )

    ADD_DEPENDENCIES(${PIP_NAME} ford_git)
ENDMACRO()
