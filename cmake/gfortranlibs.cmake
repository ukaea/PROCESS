# GFortranLibs Finder
# Author  :   K. Zarebski
# Date    :   last modified 2020-11-06
#
# Get the location of the GFortran Libraries so they can
# be included with the Python module
# TODO: Ideally we do not want to do this (likely be removed
# when f90wrap is dropped)

MACRO(GET_GFORTRANLIBS)
    MESSAGE(STATUS "[gfortran libraries]")

    SET(LIBGFORTRAN_NAME "libgfortran")

    # Get the path of libgfortran.so
    EXECUTE_PROCESS(
        COMMAND bash -c "gfortran --print-file-name ${LIBGFORTRAN_NAME}${LIBRARY_OUTPUT_SUFFIX}"

        # Returns just "libgfortran.so" rather than full path if not found
        OUTPUT_VARIABLE LIBGFORTRAN_PATH
    )
    STRING(REGEX REPLACE "\n" "" LIBGFORTRAN_PATH ${LIBGFORTRAN_PATH})

    # Error if gfortran library path not found
    IF(LIBGFORTRAN_PATH STREQUAL ${LIBGFORTRAN_NAME}${LIBRARY_OUTPUT_SUFFIX})
        MESSAGE(FATAL_ERROR "Could not retrieve location of gfortran library")
    ENDIF()

    MESSAGE(STATUS "\tlibgfortran path: ${LIBGFORTRAN_PATH}")

    # Define copy destination for gfortran library
    GET_FILENAME_COMPONENT(LIBGFORTRAN_FILE_NAME ${LIBGFORTRAN_PATH} NAME)
    SET(LIBGFORTRAN_OUTPUT ${PYTHON_LIBS_DIR}/${LIBGFORTRAN_FILE_NAME})

    # Target to copy gfortran library into Python package
    ADD_CUSTOM_TARGET(
        ${LIBGFORTRAN_NAME}
        DEPENDS ${LIBGFORTRAN_OUTPUT}
    )

    ADD_CUSTOM_COMMAND(
        OUTPUT ${LIBGFORTRAN_OUTPUT}
        COMMAND ${CMAKE_COMMAND} -E copy ${LIBGFORTRAN_PATH} ${LIBGFORTRAN_OUTPUT}
    )

    ADD_DEPENDENCIES(${LIBGFORTRAN_NAME} ${F2PY_NAME})
ENDMACRO(GET_GFORTRANLIBS)
