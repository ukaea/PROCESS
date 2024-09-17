# Build Ford Site
# Author  :   K. Zarebski (UKAEA)
# Date    :   last modified 2020-11-05
#
# Build the Python dictionaries and FORD site

MACRO(FORD)
    SET(FORD_NAME "ford")
    SET(FORD_OUTPUT ${CMAKE_BINARY_DIR}/ford_project.pickle)
    SET(FORD_INDEX_INPUT ${CMAKE_SOURCE_DIR}/documentation/ford/index_fast.md)
    MESSAGE(STATUS "[Ford]: ")
    MESSAGE(STATUS "\tInput File: ${FORD_INDEX_INPUT}")
    MESSAGE(STATUS "\tOutput: ${FORD_OUTPUT}")
    ADD_CUSTOM_TARGET(ford DEPENDS ${FORD_OUTPUT})
    ADD_CUSTOM_COMMAND(
        OUTPUT ${FORD_OUTPUT}
        COMMAND ${FORD_NAME} ${FORD_INDEX_INPUT}

        # Custom command needs to re-run whenever the Fortran add_library target
        # detects changes and is recompiled (e.g. Fortran changes are detected)
        # This keeps the dictionaries up-to-date with the Fortran source
        DEPENDS ${PROJECT_NAME}
    )
    ADD_DEPENDENCIES(${FORD_NAME} ${PIP_NAME})
ENDMACRO(FORD)

MACRO(DICTS)
    # Process the pickled Ford project object to create the dicts JSON
    SET(DICTS_NAME "dicts")
    SET(DICTS_OUTPUT_FILE ${CMAKE_BINARY_DIR}/python_fortran_dicts.json)
    SET(DICTS_PYTHON_OUT ${PYTHON_SOURCE_IO_DIR}/python_fortran_dicts.json)
    set(CREATE_DICTS_SCRIPT ${CMAKE_SOURCE_DIR}/scripts/create_dicts.py)
    ADD_CUSTOM_TARGET(${DICTS_NAME} ALL DEPENDS ${DICTS_OUTPUT_FILE} ${DICTS_PYTHON_OUT})

    ADD_CUSTOM_COMMAND(OUTPUT ${DICTS_OUTPUT_FILE} ${DICTS_PYTHON_OUT}

        # The dicts import process, install it (it will be installed again later)
        COMMAND ${Python3_EXECUTABLE} -m pip install ${CMAKE_SOURCE_DIR}

        # The create_dicts script needs to know the Fortran source dir, the pickled
        # Ford project object and the dicts.json file to output to
        COMMAND ${Python3_EXECUTABLE} ${CREATE_DICTS_SCRIPT} ${PROCESS_SRC_DIR} ${FORD_OUTPUT}
        ${DICTS_OUTPUT_FILE}
        COMMAND ${CMAKE_COMMAND} -E copy ${DICTS_OUTPUT_FILE} ${DICTS_PYTHON_OUT}

        # Custom command needs to re-run whenever the Fortran add_library target
        # detects changes and is recompiled (e.g. Fortran changes are detected)
        # This keeps the dictionaries up-to-date with the Fortran source
        # Python sources should also trigger a re-run to keep the dictionaries up-to-date
        # with Python-defined variables
        DEPENDS ${PROJECT_NAME} ${PROCESS_SRC_PATHS} ${PROCESS_PYTHON_SRC_PATHS} ${F2PY_NAME}

        # depends on f2py so that process.fortran exists
    )

    ADD_DEPENDENCIES(${DICTS_NAME} ${FORD_NAME})
ENDMACRO(DICTS)
