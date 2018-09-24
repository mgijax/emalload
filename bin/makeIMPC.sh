#!/bin/sh
#
#  makeIMPC.sh
###########################################################################
#
#  Purpose:
#
#      This script is a wrapper around the process that creates the 
#	IMPC Allele file
#
Usage="Usage: makeIMPC.sh config"
#
#  Env Vars:
#
#      See the configuration file
#
#  Inputs:  None
#
#  Outputs:
#
#      - Log file (${LOG})
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
#      2) Verify that the input file exists.
#      3) Establish the log file.
#      4) Call makeIMPC.py to QC and create the allele file
#
#  Notes:  
#	12/15/2016
#	- TR12115
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
# Establish the log file.
#
LOG=${LOG_DIAG}

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

#
# copy source input file
#
#echo "copying source input file..." >> ${LOG}
#date >> ${LOG}
#$rm -rf ${SOURCE_COPY_INPUT_FILE}
#cp ${SOURCE_INPUT_FILE} ${SOURCE_COPY_INPUT_FILE}
#STAT=$?
#checkStatus ${STAT} "Copying input file"

#
# Create the IMPC Allele input files
#
echo "" >> ${LOG}
date >> ${LOG}
echo "Create the IMPC Allele input file (makeIMPC.sh)" | tee -a ${LOG}
./makeIMPC.py 2>&1 >> ${LOG}
STAT=$?
checkStatus ${STAT} "makeIMPC.py ${CONFIG}"
#
# run noteload to add colony id notes to existing alleles
#

if [ "${LOG_DEBUG}" != "true" ]
then
    if [ -s ${CID_NOTE_FILE} ]
    then
	echo "" >> ${LOG}
	date >> ${LOG}
	${NOTELOAD}/mginoteload.csh ${EMALLOAD}/impc_noteload.config
	STAT=$?
	checkStatus ${STAT} "CID noteload ${CONFIG}"
    fi
fi

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
# run postload cleanup and email logs
#
shutDown
exit 0

