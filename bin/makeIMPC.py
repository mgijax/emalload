#!/usr/local/bin/python
#
    #  makeIMPC.py
###########################################################################
#
#  Purpose:
#
#      This script will use the records in the IMPC input file to:
#
#      1) create an Allele input file
#
#  Usage:
#
#      makeIMPC.py
#
#  Env Vars:
#
#      see impc.config
#
#  Inputs:
#
#      IMPC file ($SOURCE_COPY_INPUT_FILE)
#
#       field 1: Marker Symbol 
#       field 2: MGI Marker ID
#       field 3: ES CellLine (blank)
#	field 4: Mutation Type
#	field 5: Allele Description
#	field 6: Colony Name     
#	field 7: Colony Background Strain
#	field 8: Production Center
#	field 9: Allele MGI ID (optional)
#	field 10: Allele Symbol
#
#  Outputs: 
#
#	Allele file ($ALLELE_FILE):
#
#	field 1: MGI Marker ID
#       field 2: Mutation Type
#       field 3: Allele Description
#       field 4: Colony Name
#       field 5: Allele Description
#       field 6: Colony Name
#       field 7: Colony Background Strain
#       field 8: Allele Symbol
#       field 9: Lab Name
#       field 10: Inheritance Mode
#       field 11: Allele Type
#       field 12: Allele Subtype
#       field 13: Allele Status
#       field 14: Transmission State
#       field 15: Allele Collection

#  Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#
#  Implementation:
#
#      This script will perform following steps:
#
#      1) Initialize variables.
#      2) Open files.
#      3) perform QC
#      3) Create Allele input file
#      4) Close files.
#
#  Notes:
#       01/09/2017	sc
#	 - TR12115
#
###########################################################################

import sys 
import os
import db
import re
import string


CRT = '\n'
TAB = '\t'

# Values common to all alleles
inHeritMode = 'Not Specified'
alleleType = 'Endonuclease-mediated'
alleleSubType = 'Null/knockout'
alleleStatus = 'Approved'
transmissionState = 'Not Specified' # (no ES Cells are used)
alleleCollection = 'IMPC'

# from config
jNumber = ''
createdBy = ''

logDiagFile = None
logCurFile = None

qcFile = None
impcFile = None
alleleFile = None

# file pointers
fpLogDiag = None
fpLogCur = None

# IMPC input file
fpIMPC = None

# allele file created from IMPC Input file
fpAllele = None

#
# lookups
#

# colonyId to allele symbol in database
# {colonyID: [a1, ...an], ...}
colonyToAlleleDict = {}

# em allele lookup
# {symbol:Allele, ...}
emAlleleDict = {}

# list of lab codes in the database
# {labCode:labName, ...}
# labCode is vocab abbreviation, name is term
labCodeDict = {}

#
# QC structures
#
missingRequiredValueList = []	# 7.2.1
multiAlleleForCidList = []	# 7.2.2.1
markerIdMismatchList = []	# 7.2.2.2
alleleIdMismatchList = []	# 7.2.2.3a
alleleSymbolMismatchList = []	# 7.2.2.3b
noColIdButAlleleMgiIDList = []	# 7.2.3.1
noColIdButAlleleMatchList = []	# 7.2.3.2
labCodeNotInMgiList = []	# 7.2.3.3

# convenience object for allele information
#
class Allele:
    def __init__(self, alleleID,    # string - allele  MGI ID
            alleleSymbol,           # string - allele symbol
            markerID,               # string - marker MGI ID
	    markerSymbol,	    # string - marker symbol
	    colonyID):		    # string - pipe delim colony ID string
        self.a = alleleID
        self.s = alleleSymbol
        self.m = markerID
	self.mi = markerSymbol
	self.n = colonyID
    def toString(this):
	return '%s, %s, %s, %s, %s' % (this.a, this.s, this.m, this.mi, this.n)
#
# Purpose: Initialization
#
def initialize():
    global logDiagFile, logCurFile, qcFile
    global impcFile, alleleFile, qcFile
    global jNumber, createdBy
    global colonyToAlleleDict, emAlleleDict, labCodeDict

    logDiagFile = os.getenv('LOG_DIAG')
    logCurFile = os.getenv('LOG_CUR')
    qcFile = os.getenv('QC_FILE')
    impcFile = os.getenv('SOURCE_COPY_INPUT_FILE')
    alleleFile = os.getenv('ALLELE_FILE')
    jNumber = os.getenv('JNUMBER')
    createdBy = os.getenv('CREATEDBY')

    if openFiles() != 0:
	sys.exit(1)
    print 'querying for alleles with colony id'
    # Query for IKMC Allele Colony Name - there are multi per allele
    results = db.sql('''select distinct nc.note, a.symbol as alleleSymbol, m.symbol as markerSymbol, a1.accid as alleleID, a2.accid as markerID, a1.preferred as allelePref, a2.preferred as markerPref
	from MGI_Note n, MGI_NoteChunk nc, ALL_Allele a, MRK_Marker m, ACC_Accession a1, ACC_Accession a2
	where n._NoteType_key = 1041
	and n._Note_key = nc._Note_key
	and n._Object_key = a._Allele_key
	and a._Marker_key = m._Marker_key
	and a._Allele_key = a1._Object_key
	and a1._MGIType_key = 11
	and a1._LogicalDB_key = 1
	and a1.prefixPart = 'MGI:' 
	and a1.preferred = 1
	and a._Marker_key = a2._Object_key
	and a2._MGIType_key = 2
	and a2._LogicalDB_key = 1
	and a2.prefixPart = 'MGI:' 
	and a2.preferred = 1''', 'auto')
    for r in results:
	colonyIDString = string.strip(r['note'])
        alleleSymbol = r['alleleSymbol']
        alleleID  = r['alleleID']
	markerSymbol = r['markerSymbol']
        markerID = r['markerID']
	colonyID = r['note']
        # create allele object
        allele = Allele(alleleID, alleleSymbol, markerID, markerSymbol, colonyID)
	colonyIDList = string.split(colonyIDString, '|')
	
	# map the allele to each colony ID and create lookup
	for c in colonyIDList:
	    if c not in colonyToAlleleDict:
		colonyToAlleleDict[c] = []
		colonyToAlleleDict[c].append(allele)
	    else: # there are dup colony ID associations to alleles
		symbolList = []
		alleles =  colonyToAlleleDict[c]
		for a in alleles:
		    symbolList.append(a.s)
		if alleleSymbol not in symbolList:
		    colonyToAlleleDict[c].append(allele)
 
    # Query for em alleles and create lookup
    results = db.sql('''select a._Allele_key, a.symbol, a1.accid as alleleID, 
	    a2.accid as markerID
	from ALL_Allele a,  ACC_Accession a1, ACC_Accession a2
	where a._Allele_Type_key =  11927650
	and a.symbol like '%<em%'
	and a._Allele_key = a1._Object_key
	and a1._MGIType_key = 11
	and a1.preferred = 1
	and a1._LogicalDB_key = 1 
	and a._Marker_key = a2._Object_key
	and a2._MGIType_key = 2
	and a2.preferred = 1
	and a2._LogicalDB_key = 1''', 'auto')
    for r in results:
	alleleSymbol = r['symbol']
	alleleID  = r['alleleID']
	markerID = r['markerID']
	# create allele object
	allele = Allele(alleleID, alleleSymbol, markerID, None, None)
	emAlleleDict[alleleSymbol] = allele

    # Query for lab codes and create lookup
    results = db.sql('''select term, abbreviation from VOC_Term
	where _Vocab_key = 71''', 'auto')
    for r in results:
	labCodeDict[r['abbreviation']] = r['term']

    return 0

#
# Purpose: Open files.
#
def openFiles():
    global fpLogDiag, fpLogCur, fpQC
    global fpIMPC, fpAllele

    #
    # Open the Log Diag file; append to existing file
    #
    try:
        fpLogDiag = open(logDiagFile, 'a+')
    except:
        print 'Cannot open file: ' + logDiagFile
        return 1

    #
    # Open the Log Cur file; append to existing file
    #
    try:
        fpLogCur = open(logCurFile, 'a+')
    except:
        print 'Cannot open file: ' + logCurFile
        return 1

    #
    # Open the QC file
    #
    try:
        fpQC = open(qcFile, 'w')
    except:
        print 'Cannot open file: ' + qcFile
        return 1

    #
    # Open the IMPC file
    #
    try:
        fpIMPC = open(impcFile, 'r')
    except:
        print 'Cannot open file: ' + impcFile
        return 1

    #
    # Open the IMPC file with genotype sequence #
    #
    try:
        fpAllele = open(alleleFile, 'w')
    except:
        print 'Cannot open file: ' + alleleFile
        return 1


    return 0


#
# Purpose: Close files.
#
def closeFiles():

    if fpLogDiag:
        fpLogDiag.close()

    if fpLogCur:
        fpLogCur.close()
    
    if fpQC:
	fpQC.close()

    if fpIMPC:
        fpIMPC.close()

    if fpAllele:
        fpAllele.close()

    return 0


#
# Purpose: Read the IMPC file and QC. Create a Allele input file
#
def createAlleleFile():
    global missingRequiredValueList, multiAlleleForCidList, noColIdButAlleleMatchList
    global noColIdButAlleleMgiIDList, labCodeNotInMgiList

    lineNum = 0
    for line in fpIMPC.readlines():
	lineNum += 1
	hasError = 0
        tokens = line[:-1].split('\t')

	markerSymbol = tokens[0]
        markerID = tokens[1]
	# don't care about column 3 ES Cell Line; should be blank
	mutationType = tokens[3]
	alleleDescription = tokens[4]
	colonyID = tokens[5]
	strain = tokens[6] # background strain
	prodCtr = tokens[7]
	alleleID = tokens[8] # can be blank
	alleleSymbol = tokens[9]

	# Requirement 7.2.1
        missingDataList = []
	# report missing required values
	if markerID == '':
	    missingDataList.append('Marker ID')
	if mutationType == '':
	    missingDataList.append('Mutation Type')
	if alleleDescription == '':
	    missingDataList.append('Allele Description')
	if colonyID == '':
	    missingDataList.append('Colony ID')
	if strain == '':
	    missingDataList.append('Strain')
	if alleleSymbol == '':
	    missingDataList.append('Allele Symbol')

	if len(missingDataList):
	    missingRequiredValueList.append('%s%s%s%s%s' % (lineNum, TAB, string.join(missingDataList, ', '), TAB, line))
	    continue	# If missing fields skip remainder of QC

	# Is the colony ID in the database?
	alleleList = []
	print 'colonyID: "%s"' % colonyID
	if colonyID in colonyToAlleleDict:
	    alleleList = colonyToAlleleDict[colonyID]
	    print 'alleleList: %s' % alleleList
	    for a in alleleList:
		print a.toString()	
	
	# if colony ID associated with multiple alleles in db, report
	# Requirement 7.2.2.1
	if len(alleleList) > 1:
	    for a in alleleList:
		# report multiple alleles for a colony ID
		multiAlleleForCidList.append('%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (lineNum, TAB, a.a, TAB, a.s, TAB, a.n, TAB, a.m, TAB, a.mi, TAB, line, CRT))
	    hasError = 1
	    continue # multiple alleles for colony ID, don't do furthur checks

	# if the colony ID is not associated with an allele in the DB
	if alleleList == []:
	    # Requirement 7.2.3.1
	    if alleleID != '':
		noColIdButAlleleMgiIDList.append('%s%s%s%s' % (lineNum, TAB, line, CRT))
		hasError = 1
	    # not in database remove "(.*)" from incoming markerSymbol; compare
	    # to database
	    #print 'alleleSymbol: %s' % alleleSymbol
	    compareSymbol = re.sub('\(.*\)', '', alleleSymbol)
	    #print 'compareSymbol: %s' % compareSymbol
	    # Requirement 7.2.3.2
	    if compareSymbol in emAlleleDict:
		a = emAlleleDict[compareSymbol]
		noColIdButAlleleMatchList.append('%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (lineNum, TAB, a.a, TAB, a.s, TAB, a.n, TAB, a.m, TAB, a.mi, TAB, line, CRT))
		#print 'allele in database %s' % (compareSymbol)
		hasError = 1
	else:
	    #print 'Associated'
	    # if associated, check for alleleID, alleleSymbol, markerID consistency
	    # At this point we know there is only one allele associated with the colony ID

	    # this is the allele associated with the colony  ID in the database
	    dbAllele = alleleList[0]
	    dbAlleleID = dbAllele.a
	    dbAlleleSymbol = dbAllele.s
	    dbMarkerID = dbAllele.m
	    dbMarkerSymbol = dbAllele.mi
	    dbColonyID = dbAllele.n
	
	    # Colony ID matches IKMC Allele ColonyID Note, Marker ID Mismatch
	    # Requirement 7.2.2.2
	    #print 'markerID: %s dbMarkerID: %s' % (markerID, dbMarkerID)
	    if markerID != dbMarkerID:
                markerIdMismatchList.append('%s%s%s%s%s%s%s%s%s%s%s%s%s' % (lineNum, TAB, dbAlleleID, TAB, dbAlleleSymbol, TAB, dbColonyID, TAB, dbMarkerID, TAB, dbMarkerSymbol, TAB, line))
		hasError = 1
	    # Colony ID matches IKMC Allele ColonyID Note, Allele ID Mismatch
	    # Requirement 7.2.2.3a
	    #print 'alleleID: %s dballeleID: %s' % (alleleID, dbAlleleID)
            if alleleID != dbAlleleID:
                alleleIdMismatchList.append('%s%s%s%s%s%s%s%s%s%s%s%s%s' % (lineNum, TAB, dbAlleleID, TAB, dbAlleleSymbol, TAB, dbColonyID, TAB, dbMarkerID, TAB, dbMarkerSymbol, TAB, line))
		hasError = 1
	    # Colony ID matches IKMC Allelel ColonyID Note, Allele Symbol Mismatch
	    # Requirement 7.2.2.3b
	    #print 'alleleSymbol: %s dbAlleleSymbol: %s' % (alleleSymbol, dbAlleleSymbol)
	    if alleleSymbol != dbAlleleSymbol:
		alleleSymbolMismatchList.append('%s%s%s%s%s%s%s%s%s%s%s%s%s' % (lineNum, TAB, dbAlleleID, TAB, dbAlleleSymbol, TAB, dbColonyID, TAB, dbMarkerID, TAB, dbMarkerSymbol, TAB, line))
		hasError = 1

	# Requirement 7.2.3.3
	labCode = ''
	labName = ''
	labCodeFinder = re.compile ('\)(\w*)>')
 	match = labCodeFinder.search(alleleSymbol)
	if match:
	    labCode = match.group(1)

	if labCode not in labCodeDict:
	    labCodeNotInMgiList.append('%s%s%s' % (lineNum, TAB, line))
	    hasError = 1
	else:
	    labName = labCodeDict[labCode]
	# If no errors write out to allele file
	if hasError == 0:
	    #print 'no errors'
	    fpAllele.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (markerID, TAB, mutationType, TAB, alleleDescription, TAB, colonyID, TAB, strain, TAB, alleleSymbol, TAB, labName, TAB, inHeritMode, TAB, alleleType, TAB, alleleSubType, TAB, alleleStatus, TAB, transmissionState, TAB, alleleCollection, TAB, jNumber, TAB, createdBy, CRT))

    return 0

def writeQCReport():
    # Requirement 7.2.1
    fpQC.write('7.2.1 Missing Required Value%s%s' % (CRT, CRT))
    fpQC.write('Line#%s Missing Value(s)%sInput Line%s' % (TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(missingRequiredValueList):
	fpQC.write(string.join(missingRequiredValueList, CRT))
    fpQC.write('Total: %s' % len(missingRequiredValueList))

    # Requirement 7.2.2.1
    fpQC.write('%s%s7.2.2.1 Colony ID matches IKMC Allele Colony ID Note, Multiple Alleles Matched%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sAlleleID%sAlleleSymbol%sColonyID%sMarkerID%sMarkerSymbol%sInput Line%s' % (TAB, TAB, TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(multiAlleleForCidList):
	fpQC.write(string.join(multiAlleleForCidList, CRT))
    fpQC.write('Total: %s' % len(multiAlleleForCidList))

    # Requirement 7.2.2.2
    fpQC.write('%s%s7.2.2.2 Colony ID matches IKMC Allele Colony ID Note, Marker ID Mismatch %s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sAlleleID%sAlleleSymbol%sColonyID%sMarkerID%sMarkerSymbol%sInput Line%s' % (TAB, TAB, TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(markerIdMismatchList):
	fpQC.write(string.join(markerIdMismatchList, CRT))
    fpQC.write('Total: %s' % len(markerIdMismatchList))

    # Requirement 7.2.2.3a
    fpQC.write('%s%s7.2.2.3a Colony ID matches IKMC Allele Colony ID Note, Allele ID Mismatch %s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sAlleleID%sAlleleSymbol%sColonyID%sMarkerID%sMarkerSymbol%sInput Line%s' % (TAB, TAB, TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdMismatchList):
	fpQC.write(string.join(alleleIdMismatchList, CRT))
    fpQC.write('Total: %s' % len(alleleIdMismatchList))

    # Requirement 7.2.2.3b
    fpQC.write('%s%s7.2.2.3b Colony ID matches IKMC Allele Colony ID Note, Allele Symbol Mismatch %s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sAlleleID%sAlleleSymbol%sColonyID%sMarkerID%sMarkerSymbol%sInput Line%s' % (TAB, TAB, TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleSymbolMismatchList):
	fpQC.write(string.join(alleleSymbolMismatchList, CRT))
    fpQC.write('Total: %s' % len(alleleSymbolMismatchList))

    # Requirement 7.2.3.1
    fpQC.write('%s%s7.2.3.1 Colony ID does not match IKMC Allele Colony ID Note, Allele MGI ID Present%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Line%s' % (TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(noColIdButAlleleMgiIDList):
        fpQC.write(string.join(noColIdButAlleleMgiIDList, CRT))
    fpQC.write('Total: %s' % len(noColIdButAlleleMgiIDList))

    # Requirement 7.2.3.2
    fpQC.write('%s%s7.2.3.2 Colony ID does not match IKMC Allele Colony ID Note, Marker/LabCode/Seq# Discrepancy%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sAlleleID%sAlleleSymbol%sColonyID%sMarkerID%sMarkerSymbol%sInput Line%s' % (TAB, TAB, TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(noColIdButAlleleMatchList):
	fpQC.write(string.join(noColIdButAlleleMatchList, CRT))
    fpQC.write('Total: %s' % len(noColIdButAlleleMatchList))

    # Requirement 7.2.3.3
    fpQC.write('%s%s7.2.3.3 Colony ID does not match IKMC Allele Colony ID Note, Lab Code Not Present%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Line%s' % (TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(labCodeNotInMgiList):
	fpQC.write(string.join(labCodeNotInMgiList, CRT))
    fpQC.write('Total: %s' % len(labCodeNotInMgiList))

    return 0

#
#  MAIN
#

if initialize() != 0:
    sys.exit(1)

if createAlleleFile() != 0:
    closeFiles()
    sys.exit(1)

if writeQCReport() != 0:
    closeFiles()
    sys.exit(1)

closeFiles()
sys.exit(0)

