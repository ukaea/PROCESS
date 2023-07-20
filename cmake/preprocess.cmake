# Preprocess PROCESS Fortran sources before they are wrapped
# Author: Timothy Nunn
# Date: 2021-11-09

# Note: CMake verbatim commands insert quotes inappropriately meaning
# the preprocessor directives must (annoyingly) be inserted directly inside
# of the command

# Note: when adding new preprocessor directives, you are responsible
# for adding in appropriate ' marks for non-integer/float variables
# e.g. -Dfoo=some string of text
# will not work and must be -Dfoo="'some string of text'"
# where the "" are for CMake to identify this as a string
# and the '' are then directly inserted into the fortran along
# with the text

MACRO(PREPROCESS)
    INCLUDE(${CMAKE_SOURCE_DIR}/cmake/utilities.cmake)

    LIST(TRANSFORM PROCESS_SRCS PREPEND ${PROCESS_SRC_DIR}/ OUTPUT_VARIABLE PROCESS_SOURCE_FILES_PATH)
    LIST(TRANSFORM PROCESS_SRCS PREPEND ${CMAKE_BINARY_DIR}/ OUTPUT_VARIABLE PREPROCESSED_SOURCE_FILES_PATH)
    LIST(TRANSFORM PROCESS_SRCS PREPEND "preprocess_" OUTPUT_VARIABLE PREPROCESS_TARGET_NAMES)

    ensure_string_length(${COMMIT_MSG} 50 COMMIT_MSG) # 1502 line truncation error fix
    # f2py has a smaller line length limit so must be truncated more (this won't affect the way the code works)

    # Uses the CMake 3.17 ZIP_LISTS feature for brevity:
#    FOREACH(target_name source output IN ZIP_LISTS PREPROCESS_TARGET_NAMES PROCESS_SOURCE_FILES_PATH PREPROCESSED_SOURCE_FILES_PATH)
    # Uses features compatible with CMake 3.16 only (next 6 lines):
    LIST(LENGTH PREPROCESS_TARGET_NAMES PREPROCESS_LIST_LENGTH)
    MATH(EXPR PREPROCESS_LIST_LENGTH "${PREPROCESS_LIST_LENGTH}-1")
    FOREACH(index RANGE ${PREPROCESS_LIST_LENGTH})
        LIST(GET PREPROCESS_TARGET_NAMES ${index} target_name)
        LIST(GET PROCESS_SOURCE_FILES_PATH ${index} source)
        LIST(GET PREPROCESSED_SOURCE_FILES_PATH ${index} output)

        ADD_CUSTOM_TARGET (
            ${target_name}
            DEPENDS ${output}
        )
        ADD_CUSTOM_COMMAND (
            OUTPUT ${output}
            COMMAND gfortran -E -cpp -DINSTALLDIR="'${CMAKE_SOURCE_DIR}'" -DCOMMSG="'${COMMIT_MSG}'" -Dbranch_name="'${GIT_BRANCH}'" -Dtagno="'${GIT_TAG}'" -Duntracked=${GIT_DIFF} -Ddp=8 ${source} -o ${output}
            DEPENDS ${source} # rerun preprocessing when the source file changes
        )
    ENDFOREACH()

ENDMACRO(PREPROCESS)
