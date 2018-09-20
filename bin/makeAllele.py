#!/usr/local/bin/python

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
#       field 2: Mutation Type
#       field 3: Allele Description
#       field 4: Colony ID
#       field 5: Strain of Origin
#       field 6: Allele Symbol
#       field 7: Allele Name
#       field 8: Inheritance Mode
#       field 9: Allele Type
#       field 10: Allele Subtype
#       field 11: Allele Status
#       field 12: Transmission State
#       field 13: Allele Collection
#       field 14: J Number
#       field 15: Created By

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
# Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#
# Implementation:
#
#      This script will perform following steps:
#
#      1) Initialize variables.
#      2) Open files.
#      3) Create bcp files
#      4) Close files.
#
# History
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

#globals

#
# from configuration file
#
user = os.environ['MGD_DBUSER']
passwordFileName = os.environ['MGD_DBPASSWORDFILE']
inputFileName = os.environ['ALLELE_FILE']
outputDir = os.environ['OUTPUTDIR']
BCP_COMMAND = os.environ['PG_DBUTILS'] + '/bin/bcpin.csh'

DEBUG = 0		# if 0, not in debug mode

bcpon = 1		# can the bcp files be bcp-ed into the database?  default is yes.

diagFile = ''		# diagnostic file descriptor
errorFile = ''		# error file descriptor
inputFile = ''		# file descriptor
alleleFile = ''         # file descriptor
mutationFile = ''	# file descriptor
refFile = ''            # file descriptor
accFile = ''            # file descriptor
noteFile = ''		# file descriptor
noteChunkFile = ''	# file descriptor
annotFile = ''		# file descriptor
newAlleleFile = ''      # file descriptor

alleleTable = 'ALL_Allele'
mutationTable = 'ALL_Allele_Mutation'
refTable = 'MGI_Reference_Assoc'
accTable = 'ACC_Accession'
noteTable = 'MGI_Note'
noteChunkTable = 'MGI_NoteChunk'
annotTable = 'VOC_Annot'

alleleFileName = outputDir + '/' + alleleTable + '.bcp'
mutationFileName = outputDir + '/' + mutationTable + '.bcp'
refFileName = outputDir + '/' + refTable + '.bcp'
accFileName = outputDir + '/' + accTable + '.bcp'
noteFileName = outputDir + '/' + noteTable + '.bcp'
noteUpdateFileName = 
noteChunkFileName = outputDir + '/' + noteChunkTable + '.bcp'
annotFileName = outputDir + '/' + annotTable + '.bcp'

diagFileName = ''	# diagnostic file name
errorFileName = ''	# error file name
newAlleleFileName = ''	# output file with new accession ids

alleleKey = 0           # ALL_Allele._Allele_key
refAssocKey = 0		# MGI_Reference_Assoc._Assoc_key
accKey = 0              # ACC_Accession._Accession_key
noteKey = 0		# MGI_Note._Note_key
mgiKey = 0              # ACC_AccessionMax.maxNumericPart
annotKey = 0		# VOC_Annot._Annot_key
mgiNoteSeqNum = 1       # MGI_NoteChunk.sequenceNum
molecularNoteTypeKey = 1021   # MGI_Note._NoteType_key
colonyIdNoteTypeKey = 1041   	 # MGI_Note._NoteType_key

mgiTypeKey = 11		# Allele
mgiPrefix = 'MGI:'
annotTypeKey = 1014	# Allele SubType
qualifierKey = 1614158  # SubType annotation qualifier key (null as opposed
			# to NOT
isMixed = 0
isExtinct = 0

loaddate = loadlib.loaddate

#
# Purpose: prints error message and exits
#
def exit(
    status,          # numeric exit status (integer)
    message = None   # exit message (string)
    ):

    if message is not None:
        sys.stderr.write('\n' + str(message) + '\n')
 
    try:
        diagFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
        errorFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
        diagFile.close()
        errorFile.close()
	inputFile.close()
	
    except:
        pass

    db.useOneConnection(0)
    sys.exit(status)
 
#
# Purpose: open file descriptors
#
def initialize():
    global diagFile, errorFile, inputFile, errorFileName, diagFileName
    global alleleFile, mutationFile, refFile
    global accFile, noteFile, noteChunkFile, annotFile
    global newAlleleFile
 
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)
 
    head, tail = os.path.split(inputFileName) 

    diagFileName = outputDir + '/' + tail + '.diagnostics'
    errorFileName = outputDir + '/' + tail + '.error'
    newAlleleFileName = outputDir + '/' + tail + '.new'

    try:
        diagFile = open(diagFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % diagFileName)
		
    try:
        errorFile = open(errorFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % errorFileName)
		
    try:
        newAlleleFile = open(newAlleleFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % newAlleleFileName)

    try:
        inputFile = open(inputFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inputFileName)

    try:
        alleleFile = open(alleleFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % alleleFileName)

    try:
        mutationFile = open(mutationFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % mutationFileName)

    try:
        refFile = open(refFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % refFileName)

    try:
        accFile = open(accFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % accFileName)

    try:
        noteFile = open(noteFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % noteFileName)

    try:
        noteChunkFile = open(noteChunkFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % noteChunkFileName)

    try:
        annotFile = open(annotFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % annotFileName)

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    diagFile.write('Server: %s\n' % (db.get_sqlServer()))
    diagFile.write('Database: %s\n' % (db.get_sqlDatabase()))

    errorFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

    return 0
#
# Purpose: Close files.
#
def closeFiles():

    alleleFile.close()
    mutationFile.close()
    refFile.close()
    accFile.close()
    noteFile.close()
    noteChunkFile.close()
    annotFile.close()
    newAlleleFile.close()
 
    return 0
#
# Purpose:  sets global primary key variables
#
def setPrimaryKeys():

    global alleleKey, refAssocKey, accKey, noteKey, mgiKey, annotKey

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

    results = db.sql('select max(_Annot_key) + 1 as nextKey from VOC_Annot', 'auto')
    annotKey = results[0]['nextKey']

    return 0
#
# Purpose:  BCPs the data into the database
#
def bcpFiles():

    bcpdelim = "|"

    if DEBUG or not bcpon:
	return

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
	diagFile.write('%s\n' % bcpCmd)
	os.system(bcpCmd)

    # update serialization on mgi_reference_assoc
    db.sql(''' select setval('mgi_reference_assoc_seq', (select max(_Assoc_key) + 1
            from MGI_Reference_Assoc)) ''', None)
    db.commit()

    return 0

#
# Purpose:  processes data
#
def processFile():

    global alleleKey, refAssocKey, accKey, noteKey, mgiKey, annotKey
    global alleleLookup

    lineNum = 0
    # For each line in the input file

    for line in inputFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = line[:-1].split('\t')
        try:
	    markerID = tokens[0]
	    mutationType = tokens[1]
	    description = tokens[2]
	    colonyID = tokens[3]
	    strainOfOrigin =  tokens[4]
	    alleleSymbol =  tokens[5]
	    alleleName =   tokens[6]
	    inheritanceMode =  tokens[7]
	    alleleType = tokens[8]
	    alleleSubType  = tokens[9]
	    alleleStatus = tokens[10]
	    transmission = tokens[11]
	    collection = tokens[12]
	    jNum = tokens[13]
	    createdBy  = tokens[14]

        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

        # marker key
        markerKey = loadlib.verifyMarker(markerID, lineNum, errorFile)

        # _vocab_key = 36 (Allele Molecular Mutation)
        mutationKey = loadlib.verifyTerm('', 36, mutationType, lineNum, errorFile)

        # strains
        strainOfOriginKey = sourceloadlib.verifyStrain(strainOfOrigin, lineNum, errorFile)


        # _vocab_key = 35 (Allele Status)
        inheritanceModeKey = loadlib.verifyTerm('', 35, inheritanceMode, lineNum, errorFile)

        # _vocab_key = 38 (Allele Type)
        alleleTypeKey = loadlib.verifyTerm('', 38, alleleType, lineNum, errorFile)

        # _vocab_key = 37 (Allele Status)
        alleleStatusKey = loadlib.verifyTerm('', 37, alleleStatus, lineNum, errorFile)

	# _vocab_key = 61 (Allele Transmission)
        transmissionKey = loadlib.verifyTerm('', 61, transmission, lineNum, errorFile)

	# _vocab_key = 92
	collectionKey = loadlib.verifyTerm('', 92, collection, lineNum, errorFile)

	# _vocab_key = 73 (Marker-Allele Association Status)
	# _term_key = 4268545 (Curated)
	markerStatusKey = 4268545

	# reference
	refKey = loadlib.verifyReference(jNum, lineNum, errorFile)

        # creator
        createdByKey = loadlib.verifyUser(createdBy, lineNum, errorFile)
        if createdByKey == 0:
            continue


        # if errors, continue to next record
	# errors are stored (via loadlib) in the .error log
        if markerKey == 0 \
	        or mutationKey == 0 \
	  	or strainOfOriginKey == 0 \
                or inheritanceModeKey == 0 \
                or alleleTypeKey == 0 \
                or alleleStatusKey == 0 \
		or transmissionKey == 0 \
	 	or collectionKey == 0 \
		or refKey == 0 \
		or createdByKey == 0:
            continue

        # if no errors, process the allele

	# allele (isWildType = 0)
        alleleFile.write('%d|%s|%s|%s|%s|%s|%s|%s|%s|%s|0|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
            % (alleleKey, markerKey, strainOfOriginKey, inheritanceModeKey, alleleTypeKey, \
	    alleleStatusKey, transmissionKey, collectionKey, alleleSymbol, alleleName, \
	    isExtinct, isMixed, refKey, markerStatusKey, \
	    createdByKey, createdByKey, createdByKey, loaddate, loaddate, loaddate))

	# molecular mutation
	mutationFile.write('%s|%s|%s|%s\n' \
	    	% (alleleKey, mutationKey, loaddate, loaddate))

	# reference association
	refAssocTypeKey = 1011 # Original

	refFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
	    % (refAssocKey, refKey, alleleKey, mgiTypeKey, refAssocTypeKey, \
	       		createdByKey, createdByKey, loaddate, loaddate))
	refAssocKey = refAssocKey + 1

	#
	# allele subtype
	#
		# _vocab_key = 93 (Allele Subtype)
	alleleSubtypeKey = loadlib.verifyTerm('', 93, alleleSubType, lineNum, errorFile)

	annotFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
		% (annotKey, annotTypeKey, alleleKey, alleleSubtypeKey, \
			qualifierKey, loaddate, loaddate))
	annotKey = annotKey + 1

        # MGI Accession ID for the allelearker

        accFile.write('%s|%s%d|%s|%s|1|%d|%d|0|1|%s|%s|%s|%s\n' \
            % (accKey, mgiPrefix, mgiKey, mgiPrefix, mgiKey, alleleKey, mgiTypeKey, \
	       createdByKey, createdByKey, loaddate, loaddate))

	# storing data in MGI_Note/MGI_NoteChunk
	# molecular note

	noteFile.write('%s|%s|%s|%s|%s|%s|%s|%s\n' \
	    % (noteKey, alleleKey, mgiTypeKey, molecularNoteTypeKey, \
	       createdByKey, createdByKey, loaddate, loaddate))

	noteChunkFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
	    % (noteKey, mgiNoteSeqNum, description, createdByKey, createdByKey, loaddate, loaddate))

	noteKey = noteKey + 1

	# colony ID note
	noteFile.write('%s|%s|%s|%s|%s|%s|%s|%s\n' \
	    % (noteKey, alleleKey, mgiTypeKey, colonyIdNoteTypeKey, \
	       createdByKey, createdByKey, loaddate, loaddate))

	noteChunkFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
	    % (noteKey, 1, colonyID, createdByKey, createdByKey, loaddate, loaddate))

	noteKey = noteKey + 1

	# Print out a new text file and attach the new MGI Allele IDs 
	# as the last field

	newAlleleFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s\n' \
	% (mgi_utils.prvalue(markerID), \
	mgi_utils.prvalue(alleleSymbol), \
	mgi_utils.prvalue(alleleName), \
	mgi_utils.prvalue(alleleStatus), \
	mgi_utils.prvalue(alleleType), \
	mgi_utils.prvalue(alleleSubType), \
	mgi_utils.prvalue(collection), \
	mgi_utils.prvalue(transmission), \
	mgi_utils.prvalue(jNum), \
	mgi_utils.prvalue(strainOfOrigin), \
	mgi_utils.prvalue(mutationType), \
	mgi_utils.prvalue(inheritanceMode), \
	mgi_utils.prvalue(isMixed), \
	mgi_utils.prvalue(isExtinct), \
	mgi_utils.prvalue(createdBy), \
	mgi_utils.prvalue(mgiPrefix), mgi_utils.prvalue(mgiKey)))

        accKey = accKey + 1
        mgiKey = mgiKey + 1
        alleleKey = alleleKey + 1

    #
    # Update the AccessionMax value
    #

    if not DEBUG:
        db.sql('select * from ACC_setMax(%d)' % (lineNum), None)
	db.commit()

    return 0
#
# Main
#

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
