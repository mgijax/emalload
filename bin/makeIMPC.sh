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

#
# Establish the log file.
#
LOG=${LOG_DIAG}

#
# Create the IMPC Allele input files
#
echo "" >> ${LOG}
date >> ${LOG}
echo "Create the IMPC Allele input file (makeIMPC.sh)" | tee -a ${LOG}
./makeIMPC.py 2>&1 >> ${LOG}
STAT=$?
if [ ${STAT} -eq 1 ]
then
    exit 1
fi

exit 0

