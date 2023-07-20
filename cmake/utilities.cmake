FUNCTION(ensure_string_length str length var)
    # str: string to ensure the length of
    # length: integer length the string cannot exceed
    # var: the variable to output the truncated (or un-truncated) string into


    STRING(LENGTH ${str} str_length)
    MESSAGE(STATUS "\tensure_string_length() string length : ${str_length}")
    MESSAGE(STATUS "\tensure_string_length() desired length : ${length}")

    IF(str_length GREATER length)
        STRING(SUBSTRING ${str} 0 ${length} str_truncated)
    ELSE()
        SET(str_truncated ${str})
    ENDIF()

    MESSAGE(STATUS "\tensure_string_length() truncated string : ${str_truncated}")

    SET("${var}" ${str_truncated} PARENT_SCOPE)
ENDFUNCTION()
