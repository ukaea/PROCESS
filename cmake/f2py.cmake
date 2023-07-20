# PROCESS f2py run
# Author    :   K. Zarebski (UKAEA)
# Date      :   last modified 2020-11-09
#
# Run f2py on the given files list

MACRO(F2PY)
    EXECUTE_PROCESS(
        COMMAND bash -c "${Python3_EXECUTABLE} -c \"import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))\""
        OUTPUT_VARIABLE CMAKE_PYTHON_ABI_VERSION
    )
    STRING(STRIP ${CMAKE_PYTHON_ABI_VERSION} CMAKE_PYTHON_ABI_VERSION)

    SET(F2PY_NAME "f2py")
    SET(F2PY_MODULE_NAME "fortran")
    SET(F2PY_SIGNATURE_NAME "f2py_signatures")

    SET(F2PY_SIGNATURE_TARGET ${CMAKE_BINARY_DIR}/${F2PY_MODULE_NAME}.pyf)
    SET(F2PY_TARGET ${CMAKE_BINARY_DIR}/${F2PY_MODULE_NAME}${CMAKE_PYTHON_ABI_VERSION})
    SET(F2PY_OUTPUT ${PYTHON_MODULE_DIR}/${F2PY_MODULE_NAME}${CMAKE_PYTHON_ABI_VERSION})

    IF(CMAKE_HOST_APPLE)
        SET(MD5_COMMAND "md5 -r")
    ELSE()
        SET(MD5_COMMAND "md5sum")
    ENDIF()

    MESSAGE(STATUS "[f2py]: ")
    MESSAGE(STATUS "\tTarget: ${F2PY_TARGET}")
    MESSAGE(STATUS "\tSignature File: ${F2PY_SIGNATURE_TARGET}")

    ADD_CUSTOM_TARGET(
        ${F2PY_SIGNATURE_NAME}
        DEPENDS ${F2PY_SIGNATURE_TARGET}
    )

    ADD_CUSTOM_COMMAND(
        OUTPUT ${F2PY_SIGNATURE_TARGET}
        COMMAND ${F2PY_NAME} ${PREPROCESSED_SOURCE_FILES_PATH} --build-dir ${CMAKE_BINARY_DIR} -m ${F2PY_MODULE_NAME} -h ${F2PY_SIGNATURE_TARGET}.temp --overwrite-signature
        COMMAND /bin/bash -c "if [ -f ${F2PY_SIGNATURE_TARGET} ]; then temp_hash=($(${MD5_COMMAND} ${F2PY_SIGNATURE_TARGET}.temp)); current_hash=($(${MD5_COMMAND} ${F2PY_SIGNATURE_TARGET})); if [ \"$temp_hash\" != \"$current_hash\" ]; then echo \"Hashes are different: $temp_hash vs $current_hash\"; mv ${F2PY_SIGNATURE_TARGET}.temp ${F2PY_SIGNATURE_TARGET}; fi; else echo 'Generating definitions file for the first time'; mv ${F2PY_SIGNATURE_TARGET}.temp ${F2PY_SIGNATURE_TARGET}; fi"

        # the above is run as a bash command as if it were a shell script
        # here, we compare the hash of the existing  definitions file and new file
        # if the hashes are the same, the source code change hasnt affected our interface
        # and a rewrap is not needed
        DEPENDS ${PREPROCESSED_SOURCE_FILES_PATH} # rerun the wrapping when any of the preprocessed source files change
        VERBATIM
    )

    ADD_CUSTOM_TARGET(
        ${F2PY_NAME}
        DEPENDS ${F2PY_TARGET} ${F2PY_OUTPUT}
    )

    IF(CMAKE_HOST_APPLE)
        ADD_CUSTOM_COMMAND(
            OUTPUT ${F2PY_TARGET} ${F2PY_OUTPUT}
            COMMAND echo \"Running f2py:\"\; LDFLAGS=-rpath${CMAKE_SOURCE_DIR}/process/lib ${F2PY_NAME} -L../process/lib/ -l${PROJECT_NAME} -c ${F2PY_SIGNATURE_TARGET} --build-dir ${CMAKE_BINARY_DIR}
            COMMAND ${CMAKE_COMMAND} -E copy ${F2PY_TARGET} ${F2PY_OUTPUT}
            DEPENDS ${F2PY_SIGNATURE_TARGET}
        )
    ELSE()
        IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
            ADD_CUSTOM_COMMAND(
                OUTPUT ${F2PY_TARGET} ${F2PY_OUTPUT}
                COMMAND echo \"Running f2py:\"\; LDFLAGS=-Wl,-rpath=\\$$ORIGIN/lib ${F2PY_NAME} -L../process/lib/ -l${PROJECT_NAME} --opt="-O0" --debug -c ${F2PY_SIGNATURE_TARGET} --build-dir ${CMAKE_BINARY_DIR}
                COMMAND ${CMAKE_COMMAND} -E copy ${F2PY_TARGET} ${F2PY_OUTPUT}
                DEPENDS ${F2PY_SIGNATURE_TARGET} # rerun the wrapping when the signature file changes

                # this means that changes to the source files that do not change the
                # subroutine signature do not force a rewrap of the fortran
            )
        ELSE()
            ADD_CUSTOM_COMMAND(
                OUTPUT ${F2PY_TARGET} ${F2PY_OUTPUT}
                COMMAND echo \"Running f2py:\"\; LDFLAGS=-Wl,-rpath=\\$$ORIGIN/lib ${F2PY_NAME} -L../process/lib/ -l${PROJECT_NAME} -c ${F2PY_SIGNATURE_TARGET} --build-dir ${CMAKE_BINARY_DIR}
                COMMAND ${CMAKE_COMMAND} -E copy ${F2PY_TARGET} ${F2PY_OUTPUT}
                DEPENDS ${F2PY_SIGNATURE_TARGET}
            )
        ENDIF()
    ENDIF()

    ADD_DEPENDENCIES(${F2PY_NAME} ${PIP_NAME} ${PROJECT_NAME})
ENDMACRO()
