#!/usr/bin/bash

SOURCE_FOLDER=$PWD/source/fortran/
FLINTER_EXE=$(which flint)
BADGE_LABEL="FORTRAN Code Quality"
ANYBADGE_EXE=$(which anybadge)
QUALITY_LOG_DIR=$PWD/code_quality
SUMMARY_OUTPUT_FILE=$QUALITY_LOG_DIR/quality_summary.log

mkdir -p ${QUALITY_LOG_DIR}


# Check options and if flinter and anybadge are available
if [ -z "${FLINTER_EXE}" ]
then
    echo "ERROR: Could not find 'flint' executable, is python module 'flinter' installed?"
    exit 1
fi

if [ $# -eq 2 ]
then
    SOURCE_FOLDER=$1
    ANYBADGE_EXE="python3 $2/anybadge.py"
fi

if [ -z "${ANYBADGE_EXE}" ]
then
    echo "ERROR: Could not find 'anybadge' executable, is python module 'anybadge' installed?"
    exit 1
fi

if [ $# -eq 1 ]
then
    SOURCE_FOLDER=$1
fi


# Run the Linter on the PROCESS Source code
echo "Running flinter on Directory: ${SOURCE_FOLDER}"

${FLINTER_EXE} all-files ${SOURCE_FOLDER} | tee ${SUMMARY_OUTPUT_FILE}

SCORE=$(cat ${SUMMARY_OUTPUT_FILE} | tail -n 2 | cut -d ' ' -f 7 | cut -d '/' -f 1)
SCORE_INT=$(echo ${SCORE} | cut -d '.' -f 1)
SCORE_PERC=$(bc -l <<<"${SCORE}*10")

for file in $(ls ${SOURCE_FOLDER}); do
    echo "Running Flinter on file '${SOURCE_FOLDER}/${file}' in format mode..."
    LABEL=$(echo ${file} | cut -d '.' -f 1)
    ${FLINTER_EXE} fmt ${SOURCE_FOLDER}/${file} >> ${QUALITY_LOG_DIR}/${LABEL}_code_quality_fmt.log
    echo "Running Flinter on file '${SOURCE_FOLDER}/${file}' in pep8 style mode..."
    ${FLINTER_EXE} pep8 ${SOURCE_FOLDER}/${file} >> ${QUALITY_LOG_DIR}/${LABEL}_code_quality_pep8.log
done

# Get Badge Colour Based on Score

if [ ${SCORE_INT} -lt 3 ]
then
    BADGE_COLOR="red"
elif [ ${SCORE_INT} -lt 4 ]
then
    BADGE_COLOR="orange"
elif [ ${SCORE_INT} -lt 6 ]
then
    BADGE_COLOR="yellow"
elif [ ${SCORE_INT} -lt 8 ]
then
    BADGE_COLOR="yellowgreen"
else
    BADGE_COLOR="green"
fi

# Generate the badge
${ANYBADGE_EXE} --label="${BADGE_LABEL}" --value="${SCORE_PERC}%" --file=quality.svg --color=${BADGE_COLOR}
