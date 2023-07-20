#!/usr/bin/bash

SOURCE_FOLDER=$PWD/source/fortran/
BADGE_LABEL="Code Coverage"
ANYBADGE_EXE=$(which anybadge)
COVERAGE_DIR=$PWD/lcov_results

if [ $# -eq 2 ]
then
    ANYBADGE_EXE="python3 $2/anybadge.py"
fi

if [ $# -eq 1 ]
then
    SOURCE_FOLDER=$1
fi

if [ -z "${ANYBADGE_EXE}" ]
then
    echo "ERROR: Could not find 'anybadge' executable, is python module 'anybadge' installed?"
    exit 1
fi

LINE_COV_PERC=$(cat $COVERAGE_DIR/html/index.html | grep "headerCovTableEntryLo" | head -n 1 | cut -d ">" -f 2 | cut -d "<" -f 1 | cut -d " " -f 1)
SCORE_INT=$(echo ${LINE_COV_PERC} | cut -d '.' -f 1)

# Get Badge Colour Based on Score

if [ ${SCORE_INT} -lt 50 ]
then
    BADGE_COLOR="red"
elif [ ${SCORE_INT} -lt 60 ]
then
    BADGE_COLOR="orange"
elif [ ${SCORE_INT} -lt 70 ]
then
    BADGE_COLOR="yellow"
elif [ ${SCORE_INT} -lt 90 ]
then
    BADGE_COLOR="yellowgreen"
else
    BADGE_COLOR="green"
fi

# Generate the badge
${ANYBADGE_EXE} --label="${BADGE_LABEL}" --value="${LINE_COV_PERC}%" --file=coverage.svg --color=${BADGE_COLOR}
