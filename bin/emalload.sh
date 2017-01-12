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
echo "archiving..." >> ${LOG}
date >> ${LOG}
preload ${OUTPUTDIR}
echo "archiving complete" >> ${LOG}
date >> ${LOG}

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
# run pre-processor to do QC and create allele input file
#
${PREPROCESSOR} 2>&1 >> ${LOG}
STAT=$?
checkStatus ${STAT} "${PREPROCESSOR}"


# Create alleles
#
echo "" >> ${LOG}
date >> ${LOG}
./makeAllele.py  2>&1 >> ${LOG}
STAT=$?
checkStatus ${STAT} "makeAllele.py ${CONFIG}"

#
# Touch the "lastrun" file to note when the load was run.
#
#touch ${LASTRUN_FILE}

#
# run postload cleanup and email logs
#
shutDown
exit 0
