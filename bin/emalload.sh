#!/bin/sh
#
#  emalload.sh
###########################################################################
#
#  Purpose:
#
#      This script is a wrapper around the EMAL load process
#
Usage="Usage: emalload.sh *load.config"
#
#  Env Vars:
#
#      See the configuration file
#
#  Inputs:  None
#
#  Outputs:
#
#      - Log file (${LOG_DIAG})
#      - Log file (${LOG_PROC})
#      - Log file (${LOG_CUR})
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  Fatal error occurred
#
#  Assumes:  Nothing
#
#  Implementation:
#
#      This script will perform following steps:
#
#      1) Source the configuration file to establish the environment.
#      2) Establish the log file.
#      3) Copy the input file to the Input directory
#
#  Notes:  None
#
###########################################################################

cd `dirname $0`

if [ $# -lt 1 ]
then
    echo ${Usage}
    exit 1
fi

CONFIG=$1

.  ${CONFIG}

#
# Make sure the configuration file exists and source it.
#
if [ -f ${CONFIG} ]
then
    . ${CONFIG}
else
    echo "Missing configuration file: ${CONFIG}"
    exit 1
fi

# Establish the log file.
#
LOG=${LOG_DIAG}
rm -rf ${LOG}
rm -rf ${LOG_CUR}
touch ${LOG}
touch ${LOG_CUR}

#
#  Source the DLA library functions.
#

if [ "${DLAJOBSTREAMFUNC}" != "" ]
then
    if [ -r ${DLAJOBSTREAMFUNC} ]
    then
        . ${DLAJOBSTREAMFUNC}
    else
        echo "Cannot source DLA functions script: ${DLAJOBSTREAMFUNC}" | tee -a ${LOG}
        exit 1
    fi
else
    echo "Environment variable DLAJOBSTREAMFUNC has not been defined." | tee -a ${LOG}
    exit 1
fi

#
# createArchive
#
preload ${OUTPUTDIR}
rm -f ${OUTPUTDIR}/*

#
# There should be a "lastrun" file in the input directory that was created
# the last time the load was run for this input file. If this file exists
# and is more recent than the input file, the load does not need to be run.
#
LASTRUN_FILE=${INPUTDIR}/lastrun
if [ -f ${LASTRUN_FILE} ]
then
    if test ${LASTRUN_FILE} -nt ${SOURCE_INPUT_FILE}; then
       echo "" >> ${LOG_CUR} 2>&1
       echo "LOAD SKIPPED: No new input file: ${SOURCE_INPUT_FILE}" >> ${LOG_CUR} 2>&1
       STAT=0
       checkStatus ${STAT} "LOAD SKIPPED: No new input file ${SOURCE_INPUT_FILE}"
       shutDown
       exit 0
    fi
fi

#
# copy source input file
#
echo "copying source input file..." >> ${LOG}
date >> ${LOG}
rm -rf ${SOURCE_COPY_INPUT_FILE}
cp ${SOURCE_INPUT_FILE} ${SOURCE_COPY_INPUT_FILE}
STAT=$?
checkStatus ${STAT} "Copying input file"

#
# FUNCTION: Check for duplicate lines in an input file and write the lines
#           to the sanity report.
#
checkDupLines ()
{
    FILE=$1    # The input file to check
    REPORT=$2  # The sanity report to write to

    echo "Duplicate Lines" >> ${REPORT}
    echo "---------------" >> ${REPORT}
    sort ${FILE} | uniq -d > ${TMP_FILE1}
    cat ${TMP_FILE1} >> ${REPORT}
    if [ `cat ${TMP_FILE1} | wc -l` -eq 0 ]
    then
        return 0
    else
        return 1
    fi
}

#
# FUNCTION: Check for lines with missing columns in an input file and write
#           the line numbers to the sanity report.
#
checkColumns ()
{
    FILE=$1         # The input file to check
    REPORT=$2       # The sanity report to write to
    NUM_COLUMNS=$3  # The number of columns expected in each input record

    echo "" >> ${REPORT}
    echo "" >> ${REPORT}
    echo "Lines With Missing Columns" >> ${REPORT}
    echo "--------------------------" >> ${REPORT}
    ${PYTHON} ${EMALLOAD}/bin/checkColumns.py ${FILE} ${NUM_COLUMNS} > ${TMP_FILE2}
    cat ${TMP_FILE2} >> ${REPORT}
    if [ `cat ${TMP_FILE2} | wc -l` -eq 0 ]
    then
        return 0
    else
        return 1
    fi

}

#
# FUNCTION: Check an input file to make sure it has a minimum number of lines.
#
checkLineCount ()
{
    FILE=$1        # The input file to check
    REPORT=$2      # The sanity report to write to
    NUM_LINES=$3   # The minimum number of lines expected in the input file

    COUNT=`cat ${FILE} | wc -l | sed 's/ //g'`
    if [ ${COUNT} -lt ${NUM_LINES} ]
    then
        echo "" >> ${REPORT}
        echo "" >> ${REPORT}
        echo "**** WARNING ****" >> ${REPORT}
        echo "${FILE} has ${COUNT} lines." >> ${REPORT}
        echo "Expecting at least ${NUM_LINES} lines." >> ${REPORT}
        return 1
    else
        return 0
    fi
}

#
# Run sanity checks on the input file.
#
#
# Create temporary files and make sure it is removed when this script
# terminates.
#
TMP_FILE1=/tmp/`basename $0`.$$
trap "rm -f ${TMP_FILE}" 0 1 2 15

TMP_FILE2=/tmp/`basename $0`.$$
trap "rm -f ${TMP_FILE}" 0 1 2 15

echo "" >> ${LOG}
date >> ${LOG}
echo "Run sanity checks on the input file" >> ${LOG}
SANITY_ERROR=0

# reset SANITY_RPT
rm -f ${SANITY_RPT}; >${SANITY_RPT}
 
checkDupLines  ${SOURCE_COPY_INPUT_FILE} ${SANITY_RPT}
if [ $? -ne 0 ]
then
    SANITY_ERROR=1
fi

checkColumns  ${SOURCE_COPY_INPUT_FILE} ${SANITY_RPT} ${NUM_COLUMNS}
if [ $? -ne 0 ]
then
    SANITY_ERROR=1
fi

checkLineCount  ${SOURCE_COPY_INPUT_FILE} ${SANITY_RPT} ${FILE_MIN_SIZE}
if [ $? -ne 0 ]
then
    SANITY_ERROR=1
fi

if [ ${SANITY_ERROR} -ne 0 ]
then
    echo "Sanity errors detected. See ${SANITY_RPT}" | tee -a ${LOG}
    shutDown
    exit 1
fi

#
# run pre-processor to do QC and create allele input file
#
${PYTHON} ${PREPROCESSOR} 2>&1 >> ${LOG}
STAT=$?
checkStatus ${STAT} "${PREPROCESSOR}"

#
# run noteload to add colony id notes to existing alleles
# check to make sure the file exists and is not size 0
#
if [ "${LOG_DEBUG}" != "true" ]
then
    echo 'executing noteload' >> ${LOG}
    if [ -s ${CID_NOTE_FILE} ]
    then
	echo 'cid file is not empty' >> ${LOG}
	echo "" >> ${LOG}
	date >> ${LOG}
	${NOTELOAD}/mginoteload.csh ${EMALLOAD}/impc_noteload.config
	STAT=$?
	checkStatus ${STAT} "CID noteload ${CONFIG}"
    fi
fi
#
# Create alleles
#
echo "" >> ${LOG}
date >> ${LOG}
${PYTHON} ./makeAllele.py  2>&1 >> ${LOG}
STAT=$?
checkStatus ${STAT} "makeAllele.py ${CONFIG}"

#
# Archive a copy of the input file, adding a timestamp suffix.
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo "Archive input file" >> ${LOG_DIAG}
TIMESTAMP=`date '+%Y%m%d.%H%M'`
ARC_FILE=`basename ${SOURCE_INPUT_FILE}`.${TIMESTAMP}
cp -p ${SOURCE_INPUT_FILE} ${ARCHIVEDIR}/${ARC_FILE}

#
# Touch the "lastrun" file to note when the load was run.
#
touch ${LASTRUN_FILE}

#
# run postload cleanup and email logs
#
shutDown
exit 0
