#
# Program: makeAllele.py
#
# Original Author: sc
#
# Purpose:
#
#	To load new CRISPR Alleles into MGI
#
# Usage:
#	makeAllele.py
#
# Envvars:
#	see config file
# Inputs:
#
#	A tab-delimited file in the format:
#
#       field 1: MGI Marker ID
#	field 2: Marker Symbol (For IMPC allele report)
#       field 3: Mutation Type (IMPC Allele Type)
#       field 4: Allele Description
#       field 5: Colony ID
#       field 6: Strain of Origin
#       field 7: Allele Symbol
#       field 8: Allele Name
#       field 9: Inheritance Mode
#       field 10: Allele Type (IMPC Allele Class)
#       field 11: Allele Subtype
#       field 12: Allele Status
#       field 13: Transmission State
#       field 14: Allele Collection
#       field 15: J Number
#       field 16: Created By

#
# Outputs:
#
#	BCP files:
#	  ALL_Allele
#	  ALL_Allele_Mutation
#	  MGI_Reference_Assoc
#	  ACC_Accession
#	  MGI_Note	(molecular and colony ID)
#	  MGI_NoteChunk (molecular and colony ID)
#	  VOC_Annot (allele subType)
#
#	Log Files
#	  Diagnostic Log
#	  Error Log
#
# Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#
# Implementation:
#
#      This script will perform following steps:
#
#      1) Initialize variables
#      2) Open files.
#      3) Create bcp files
#      4) Close files.
#      5) Execute bcp
#
# History
#
# 09/2018	sc
#	- TR12115 rewrite of load now that we have a 
#		real input file
#
# 01/09/2016	sc
#	- TR12115 CRISPR Allele load
#

import sys
import os
import db
import mgi_utils
import loadlib
import sourceloadlib

#
# from configuration file
#
inputFileName = os.getenv('ALLELE_FILE')
outputDir = os.getenv('OUTPUTDIR')
BCP_COMMAND = os.getenv('PG_DBUTILS') + '/bin/bcpin.csh'
newAlleleRptFileName = os.getenv('NEW_ALLELE_RPT')

DEBUG = os.getenv('LOG_DEBUG')	# if 'true', in debug mode and  bcp files 
                                # will not be bcp-ed into the database. Default is 'false'.

fpDiagFile = ''		# diagnostic file descriptor
fpErrorFile = ''	# error file descriptor
fpInputFile = ''	# input file descriptor
fpAlleleFile = ''       # allele bcp file descriptor
fpMutationFile = ''	# allele mutation bcp file descriptor
fpRefFile = ''          # reference assoc bcp file descriptor
fpAccFile = ''          # accession bcp file descriptor
fpNoteFile = ''		# note bcp file descriptor
fpNoteChunkFile = ''	# note chunk bcp file descriptor
fpAnnotFile = ''	# annotation bcp file descriptor
fpNewAlleleRptFile = '' # new allele report file descriptor

alleleTable = 'ALL_Allele'
mutationTable = 'ALL_Allele_Mutation'
refTable = 'MGI_Reference_Assoc'
accTable = 'ACC_Accession'
noteTable = 'MGI_Note'
noteChunkTable = 'MGI_NoteChunk'
annotTable = 'VOC_Annot'

# reference types
origRefTypeKey = 1011 # Original
molRefTypeKey = 1012  # Molecular

alleleFileName = outputDir + '/' + alleleTable + '.bcp'
mutationFileName = outputDir + '/' + mutationTable + '.bcp'
refFileName = outputDir + '/' + refTable + '.bcp'
accFileName = outputDir + '/' + accTable + '.bcp'
noteFileName = outputDir + '/' + noteTable + '.bcp'
noteChunkFileName = outputDir + '/' + noteChunkTable + '.bcp'
annotFileName = outputDir + '/' + annotTable + '.bcp'

diagFileName = ''	# diagnostic file name
errorFileName = ''	# error file name

alleleKey = 0           # ALL_Allele._Allele_key
refAssocKey = 0		# MGI_Reference_Assoc._Assoc_key
accKey = 0              # ACC_Accession._Accession_key
noteKey = 0		# MGI_Note._Note_key
mgiKey = 0              # ACC_AccessionMax.maxNumericPart
annotKey = 0		# VOC_Annot._Annot_key
mgiNoteSeqNum = 1       # MGI_NoteChunk.sequenceNum
molecularNoteTypeKey = 1021      # MGI_Note._NoteType_key for molecular note
colonyIdNoteTypeKey = 1041   	 # MGI_Note._NoteType_key for colony id note

mgiTypeKey = 11		# Allele
mgiPrefix = 'MGI:'
annotTypeKey = 1014	# Allele SubType
qualifierKey = 1614158  # SubType annotation qualifier key (null as opposed
                        # to NOT
isMixed = 0	 	# allele isMixed value
isExtinct = 0		# allele isExtinct value

loaddate = loadlib.loaddate

def exit(
    # Purpose: prints error 'message' if it is not None
    #     writes to log files and exits with 'status'
    # Returns: nothing
    # Assumes: Nothing
    # Effects: Exits with 'status'

    status,          # numeric exit status (integer)
    message = None   # exit message (str.
    ):

    if message is not None:
        sys.stderr.write('\n' + str(message) + '\n')
 
    try:
        fpDiagFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
        fpErrorFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
        fpDiagFile.close()
        fpErrorFile.close()
        fpInputFile.close()
        
    except:
        pass

    db.useOneConnection(0)
    sys.exit(status)
 
def initialize():
    # Purpose: open file descriptors; write timestamps to log files
    # Returns: 1 if error, else 0
    # Assumes: Nothing
    # Effects: Sets global variables, exits if a file can't be opened,
    #  creates files in the file system

    global fpDiagFile, fpErrorFile, fpInputFile, errorFileName, diagFileName
    global fpAlleleFile, fpMutationFile, fpRefFile
    global fpAccFile, fpNoteFile, fpNoteChunkFile, fpAnnotFile
    global fpNewAlleleRptFile
 
    db.useOneConnection(1)
    #db.set_sqlUser(user)
    #db.set_sqlPasswordFromFile(passwordFileName)
 
    head, tail = os.path.split(inputFileName) 

    diagFileName = outputDir + '/' + tail + '.diagnostics'
    errorFileName = outputDir + '/' + tail + '.error'

    try:
        fpDiagFile = open(diagFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % diagFileName)
                
    try:
        fpErrorFile = open(errorFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % errorFileName)
                
    try:
        fpNewAlleleRptFile = open(newAlleleRptFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % newAlleleRptFileName)

    try:
        fpInputFile = open(inputFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inputFileName)

    try:
        fpAlleleFile = open(alleleFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % alleleFileName)

    try:
        fpMutationFile = open(mutationFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % mutationFileName)

    try:
        fpRefFile = open(refFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % refFileName)

    try:
        fpAccFile = open(accFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % accFileName)

    try:
        fpNoteFile = open(noteFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % noteFileName)

    try:
        fpNoteChunkFile = open(noteChunkFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % noteChunkFileName)

    try:
        fpAnnotFile = open(annotFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % annotFileName)

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    fpDiagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    fpDiagFile.write('Server: %s\n' % (db.get_sqlServer()))
    fpDiagFile.write('Database: %s\n' % (db.get_sqlDatabase()))

    fpErrorFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

    return 0

def closeFiles():
    # Purpose: Close all file descriptors
    # Returns: 1 if error, else 0
    # Assumes: all file descriptors were initialized
    # Effects: Nothing
    # Throws: Nothing

    try:
        fpAlleleFile.close()
        fpMutationFile.close()
        fpRefFile.close()
        fpAccFile.close()
        fpNoteFile.close()
        fpNoteChunkFile.close()
        fpAnnotFile.close()
        fpNewAlleleRptFile.close()
    except:
        return 1 
    return 0

def setPrimaryKeys():
    # Purpose: sets global primary key variables
    # Returns: 1 if error, else 0
    # Assumes: database connection exists
    # Effects: Nothing
    # Throws: Nothing

    global alleleKey, refAssocKey, accKey, noteKey, mgiKey, annotKey, alleleMutationKey

    results = db.sql(''' select nextval('all_allele_mutation_seq') as nextKey ''', 'auto')
    alleleMutationKey = results[0]['nextKey']

    results = db.sql('select max(_Allele_key) + 1 as nextKey from ALL_Allele', 'auto')
    alleleKey = results[0]['nextKey']

    results = db.sql('select max(_Assoc_key) + 1 as nextKey from MGI_Reference_Assoc', 'auto')
    refAssocKey = results[0]['nextKey']

    results = db.sql('select max(_Accession_key) + 1 as nextKey from ACC_Accession', 'auto')
    accKey = results[0]['nextKey']

    results = db.sql('select max(_Note_key) + 1 as nextKey from MGI_Note', 'auto')
    noteKey = results[0]['nextKey']

    results = db.sql('''select max(maxNumericPart) + 1 as nextKey from ACC_AccessionMax 
        where prefixPart = '%s' ''' % (mgiPrefix), 'auto')
    mgiKey = results[0]['nextKey']

    results = db.sql(''' select nextval('voc_annot_seq') as nextKey ''', 'auto')
    annotKey = results[0]['nextKey']

    return 0

def bcpFiles():
    # Purpose: BCPs the data into the database
    # Returns: 1 if error,  else 0
    # Assumes: connection to the database
    # Effects: copies data into the db
    # Throws: Nothing

    if DEBUG == 'true':
        return 0

    closeFiles()

    bcpI = '%s %s %s' % (BCP_COMMAND, db.get_sqlServer(), db.get_sqlDatabase())
    bcpII = '"|" "\\n" mgd'

    bcp1 = '%s %s "/" %s %s' % (bcpI, alleleTable, alleleFileName, bcpII)
    bcp2 = '%s %s "/" %s %s' % (bcpI, mutationTable, mutationFileName, bcpII)
    bcp3 = '%s %s "/" %s %s' % (bcpI, refTable, refFileName, bcpII)
    bcp4 = '%s %s "/" %s %s' % (bcpI, accTable, accFileName, bcpII)
    bcp5 = '%s %s "/" %s %s' % (bcpI, noteTable, noteFileName, bcpII)
    bcp6 = '%s %s "/" %s %s' % (bcpI, noteChunkTable, noteChunkFileName, bcpII)
    bcp7 = '%s %s "/" %s %s' % (bcpI, annotTable, annotFileName, bcpII)

    db.commit()

    for bcpCmd in [bcp1, bcp2, bcp3, bcp4, bcp5, bcp6, bcp7]:
        fpDiagFile.write('%s\n' % bcpCmd)
        os.system(bcpCmd)

    # update serialization on mgi_reference_assoc
    db.sql(''' select setval('mgi_reference_assoc_seq', (select max(_Assoc_key) + 1
            from MGI_Reference_Assoc)) ''', None)
    db.commit()

    # update voc_annot_seq auto-sequence
    db.sql(''' select setval('voc_annot_seq', (select max(_Annot_key) from VOC_Annot)) ''', None)

    # update all_allele_mutation_seq auto-sequence
    db.sql(''' select setval('all_allele_mutation_seq', (select max(_Assoc_key) from ALL_Allele_Mutation)) ''', None)

    return 0

def processFile():
    # Purpose: Read the input file, resolve values to keys. Create bcp files
    # Returns: 1 if error,  else 0
    # Assumes: file descriptors have been initialized
    # Effects: exits if the line does not have 15 columns
    # Throws: Nothing

    global alleleKey, refAssocKey, accKey, noteKey, mgiKey, annotKey
    global alleleLookup

    lineNum = 0
    # For each line in the input file

    for line in fpInputFile.readlines():

        error = 0
        lineNum = lineNum + 1
        print('%s: %s' % (lineNum, line))
        # Split the line into tokens
        tokens = line[:-1].split('\t')
        try:
            markerID = tokens[0]
            markerSymbol = tokens[1]
            mutationType = tokens[2] 	# IMPC allele type
            description = tokens[3]
            colonyID = tokens[4]
            strainOfOrigin =  tokens[5]
            alleleSymbol =  tokens[6]
            alleleName =   tokens[7]
            inheritanceMode =  tokens[8]
            alleleType = tokens[9] 	# IMPC allele class
            alleleSubType  = tokens[10]
            alleleStatus = tokens[11]
            transmission = tokens[12]
            collection = tokens[13]
            jNum = tokens[14]
            createdBy  = tokens[15]

        except:
            print('exiting with invalid line')
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

        print('validating data and getting keys')
        # marker key
        markerKey = loadlib.verifyMarker(markerID, lineNum, fpErrorFile)

        # _vocab_key = 36 (Allele Molecular Mutation)
        mutationList = str.split(mutationType, ';')
        if len(mutationList) > 1:
           print('mutationList: %s' % mutationList)
        mutationKeyList = []
        for m in mutationList:
            mutationKey = loadlib.verifyTerm('', 36, m, lineNum, fpErrorFile)
            if mutationKey != 0:
                mutationKeyList.append(mutationKey)
        if len(mutationKeyList) > 1:
            print('mutationKeyList: %s' % mutationKeyList)
        # strains
        strainOfOriginKey = sourceloadlib.verifyStrain(strainOfOrigin, lineNum, fpErrorFile)


        # _vocab_key = 35 (Allele Inheritance Mode)
        inheritanceModeKey = loadlib.verifyTerm('', 35, inheritanceMode, lineNum, fpErrorFile)

        # _vocab_key = 38 (Allele Type)
        alleleTypeKey = loadlib.verifyTerm('', 38, alleleType, lineNum, fpErrorFile)

        # _vocab_key = 93 (Allele Subtype)
        subTypeList = str.split(alleleSubType, ';')
        if len(subTypeList) > 1:
           print('subTypeList: %s' % subTypeList)
        subTypeKeyList = []
        for s in subTypeList:
            if s != '': # if we have a subtype, get it's key
                subTypeKey = loadlib.verifyTerm('', 93, s, lineNum, fpErrorFile)
                if subTypeKey != 0:
                    subTypeKeyList.append(subTypeKey)
        if len(subTypeKeyList) > 1:
            print('subTypeKeyList: %s' % subTypeKeyList)

        # _vocab_key = 37 (Allele Status)
        alleleStatusKey = loadlib.verifyTerm('', 37, alleleStatus, lineNum, fpErrorFile)

        # _vocab_key = 61 (Allele Transmission)
        transmissionKey = loadlib.verifyTerm('', 61, transmission, lineNum, fpErrorFile)

        # _vocab_key = 92
        collectionKey = loadlib.verifyTerm('', 92, collection, lineNum, fpErrorFile)

        # _vocab_key = 73 (Marker-Allele Association Status)
        # _term_key = 4268545 (Curated)
        markerStatusKey = 4268545

        # reference
        refKey = loadlib.verifyReference(jNum, lineNum, fpErrorFile)

        # creator
        createdByKey = loadlib.verifyUser(createdBy, lineNum, fpErrorFile)
        if createdByKey == 0:
            continue

        print('checking for missing data')
        # if errors, continue to next record
        # errors are stored (via loadlib) in the .error log
        if markerKey == 0 \
                or mutationKeyList == [] \
                or strainOfOriginKey == 0 \
                or inheritanceModeKey == 0 \
                or alleleTypeKey == 0 \
                or alleleStatusKey == 0 \
                or transmissionKey == 0 \
                or collectionKey == 0 \
                or refKey == 0 \
                or createdByKey == 0:
            print('missing data, skipping this line')
            continue

        # if no errors, process the allele
        print('writing to allele file')
        # allele (isWildType = 0)
        fpAlleleFile.write('%d|%s|%s|%s|%s|%s|%s|%s|%s|%s|0|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
            % (alleleKey, markerKey, strainOfOriginKey, inheritanceModeKey, alleleTypeKey, \
            alleleStatusKey, transmissionKey, collectionKey, alleleSymbol, alleleName, \
            isExtinct, isMixed, refKey, markerStatusKey, \
            createdByKey, createdByKey, createdByKey, loaddate, loaddate, loaddate))

        # molecular mutation
        for mutationKey in mutationKeyList:
            fpMutationFile.write('%s|%s|%s|%s|%s\n' \
                % (alleleMutationKey, alleleKey, mutationKey, loaddate, loaddate))
            alleleMutationKey += 1

        # reference associations

        # Original
        fpRefFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
            % (refAssocKey, refKey, alleleKey, mgiTypeKey, origRefTypeKey, \
                        createdByKey, createdByKey, loaddate, loaddate))
        refAssocKey = refAssocKey + 1

        # Molecular
        fpRefFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
            % (refAssocKey, refKey, alleleKey, mgiTypeKey, molRefTypeKey, \
                        createdByKey, createdByKey, loaddate, loaddate))
        refAssocKey = refAssocKey + 1

        # allele subtype
        for subTypeKey in subTypeKeyList:
            fpAnnotFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
                    % (annotKey, annotTypeKey, alleleKey, subTypeKey, \
                            qualifierKey, loaddate, loaddate))
            annotKey = annotKey + 1

        # MGI Accession ID for the allele
        alleleID = '%s%s' % (mgiPrefix, mgiKey)
        fpAccFile.write('%s|%s|%s|%s|1|%d|%d|0|1|%s|%s|%s|%s\n' \
            % (accKey, alleleID, mgiPrefix, mgiKey, alleleKey, mgiTypeKey, \
               createdByKey, createdByKey, loaddate, loaddate))

        # storing data in MGI_Note/MGI_NoteChunk
        # molecular note

        fpNoteFile.write('%s|%s|%s|%s|%s|%s|%s|%s\n' \
            % (noteKey, alleleKey, mgiTypeKey, molecularNoteTypeKey, \
               createdByKey, createdByKey, loaddate, loaddate))

        fpNoteChunkFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
            % (noteKey, mgiNoteSeqNum, description, createdByKey, createdByKey, loaddate, loaddate))

        noteKey = noteKey + 1

        # colony ID note
        fpNoteFile.write('%s|%s|%s|%s|%s|%s|%s|%s\n' \
            % (noteKey, alleleKey, mgiTypeKey, colonyIdNoteTypeKey, \
               createdByKey, createdByKey, loaddate, loaddate))

        fpNoteChunkFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
            % (noteKey, 1, colonyID, createdByKey, createdByKey, loaddate, loaddate))

        noteKey = noteKey + 1

        # Print out a new text file and attach the new MGI Allele IDs 
        # as the last field

        fpNewAlleleRptFile.write('%s\t%s\t%s\t%s\t%s\t%s\n' \
        % (mgi_utils.prvalue(alleleID), \
        mgi_utils.prvalue(alleleSymbol), \
        mgi_utils.prvalue(alleleName), \
        mgi_utils.prvalue(markerID), \
        mgi_utils.prvalue(markerSymbol), \
        mgi_utils.prvalue(colonyID)))

        accKey = accKey + 1
        mgiKey = mgiKey + 1
        alleleKey = alleleKey + 1

    #
    # Update the AccessionMax value
    #
    print('DEBUG: %s' % DEBUG)
    if DEBUG == 'false':
        db.sql('select * from ACC_setMax(%d)' % (lineNum), None)
        db.commit()

    return 0

#
#  MAIN
#

if initialize() != 0:
    sys.exit(1)

if setPrimaryKeys() != 0:
    sys.exit(1)

if processFile() != 0:
    sys.exit(1)

if bcpFiles() != 0:
    sys.exit(1)

sys.exit(0)
