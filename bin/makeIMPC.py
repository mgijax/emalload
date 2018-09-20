#!/usr/local/bin/python
#
###########################################################################
#
#
# Program: makeIMPC.py
#
# Original Author: sc
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
#       field 1: Marker Symbol (not used by load)
#       field 2: MGI Marker ID
#       field 3: ES Cell Line (not used by load)
#       field 4: Colony ID
#       field 5: Colony Background Strain
#       field 5: Project Name (not used by load)
#       field 6: Production Center (not used by load)
#       field 7: Allele Class (aka allele type)
#       field 8: Allele Type (aka mutation type)
#       field 9: Allele Subtype
#       field 10: Allele Description
#       field 11: Allele Name (aka allele superscript)
#       field 12: MGI Allele ID (can be blank)

#  Outputs: 
#
#	Allele file ($ALLELE_FILE):
#
#      field 1: MGI Marker ID
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
#      field 14: J Number
#      field 15: Created By

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
#	09/2018		sc
#	 - reimplementation based on new file
#
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
# from config
inHeritMode = ''
transmissionState = ''
alleleStatus =   '' 
alleleCollection = ''
jNumber = ''
createdBy = ''

# temp - will be getting this from the file
alleleSubType = ''

logDiagFile = None
logCurFile = None

qcFile = None
impcFile = None
alleleFile = None
notloadFile = None

# file pointers
fpLogDiag = None
fpLogCur = None

# IMPC input file
fpIMPC = None

# allele file created from IMPC Input file
fpAllele = None

# allele colony ID note file
# allele MGI ID\tColony ID
fpNoteload = None

#
# lookups
#

# colonyId to allele symbol in database
# {colonyID: [a1, ...an], ...}
colonyToAlleleDict = {}

# allele lookup by symbol
# {symbol:Allele, ...}
alleleBySymbolDict = {}

# allele lookkup by mgiID
# {mgiID:Allele, ...}
alleleByIDDict = {}

# list of lab codes in the database
# {labCode:labName, ...}
# labCode is vocab abbreviation, laabName is term
labCodeDict = {}

# marker ID to marker name lookup (for Allele Name construction)
# {markerID: markerName|symbol, ...}
markerDict = {}

# {alleleKey:colonyIDNote, ...}
colonyDict = {}

# list of standard and non-private strains in the database
strainList = []

# template for creating allelel name for new alleles
# marker name, sequenceNum, lab code name
alleleNameTemplate = '%s; endonuclease-mediated mutation %s, %s'

#
# QC lists for reporting errors
#
missingRequiredValueList = []	
labCodeNotInMgiList = []	
markerIdNotInMgiList = []	
strainNotInMgiList = []		
unknownAlleleClassList = []     
unknownAlleleTypeList = []      
alleleIdNotInMGIList = []	
alleleIdMatchAlleleStatusDiscrepList = []	
alleleIdMatchMarkerIdMismatchList = []		
alleleIdMatchAlleleSSMismatchList = []  	
alleleIdMatchColonyIDMismatchList = []  	
alleleIdMatchColonyIdMatchToMultiList = []	
alleleIdMatchColonyIdMatchToDiffAlleleList = []	
cidMatchToMultiList = []		
cidMatchMarkerIdMismatchList = []		
cidMatchAlleleSSMismatchList = []               
cidMatchAlleleStatusDiscrepList = []
symbolMatchAlleleStatusDiscrepList = []		
symbolMatchColonyIdMismatchList = []		
symbolMatchMultiAlleleList  = []	

# not an error, just info for testing, maybe curators would like to see too
addCidSymbolMatchList = []
addCidAlleleIDMatchList = []

#
# convenience object for allele information
#
class Allele:
    def __init__(self, alleleID,    # string - allele  MGI ID
            alleleSymbol,           # string - allele symbol
	    alleleStatus,	    # string - allele status
	    alleleType,		    # string - allele type
            markerID,               # string - marker MGI ID
	    markerSymbol,	    # string - marker symbol
	    colonyID):		    # string - pipe delim colony ID string
        self.aid = alleleID
        self.asym = alleleSymbol
	self.ast = alleleStatus
	self.at = alleleType
        self.mid = markerID
	self.ms = markerSymbol
	self.cid = colonyID
    def toString(this):
	return '%s, %s, %s, %s, %s, %s' % (this.aid, this.asym, this.ast, this.mid, this.ms, this.cid)
#
# Purpose: Initialization
#
def initialize():
    global logDiagFile, logCurFile, qcFile, impcFile, alleleFile, noteloadFile
    global jNumber, createdBy, alleleSubType, inHeritMode, alleleStatus
    global transmissionState, alleleCollection, strainList
    global colonyToAlleleDict, alleleBySymbolDict, labCodeDict, markerDict, colonyDict

    logDiagFile = os.getenv('LOG_DIAG')
    logCurFile = os.getenv('LOG_CUR')
    qcFile = os.getenv('QC_FILE')
    impcFile = os.getenv('SOURCE_COPY_INPUT_FILE')
    alleleFile = os.getenv('ALLELE_FILE')
    noteloadFile = os.getenv('CID_NOTE_FILE')
    jNumber = os.getenv('JNUMBER')
    createdBy = os.getenv('CREATEDBY')
    alleleSubType = os.getenv('ALLELE_SUBTYPE') # temp, will be getting this from the file
    inHeritMode = os.getenv('INHERIT_MODE')
    alleleStatus = os.getenv('ALLELE_STATUS')
    transmissionState = os.getenv('TRANSMISSION_STATE')
    alleleCollection = os.getenv('ALLELE_COLLECTION')

    if openFiles() != 0:
	sys.exit(1)
    
    # Query for IKMC Allele Colony Name - there are multi per allele
    results = db.sql('''select distinct nc.note as cidNote, a.symbol as alleleSymbol, t.term as alleleStatus, 
	    t2.term as alleleType, m.symbol as markerSymbol, a1.accid as alleleID, a2.accid as markerID, 
	    a1.preferred as allelePref, a2.preferred as markerPref
	from MGI_Note n, MGI_NoteChunk nc, ALL_Allele a, MRK_Marker m, ACC_Accession a1, ACC_Accession a2, 
	    VOC_Term t, VOC_Term t2
	where n._NoteType_key = 1041
	and n._Note_key = nc._Note_key
	and n._Object_key = a._Allele_key
	and a._Marker_key = m._Marker_key
	and a._Allele_Status_key = t._Term_key
	and a._Allele_Type_key = t2._Term_key
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
	colonyIDString = string.strip(r['cidNote'])
        alleleSymbol = r['alleleSymbol']
	alleleStatus = r['alleleStatus']
	alleleType = r['alleleType']
        alleleID  = r['alleleID']
	markerSymbol = r['markerSymbol']
        markerID = r['markerID']
        # create allele object
        allele = Allele(alleleID, alleleSymbol, alleleStatus, alleleType, markerID, markerSymbol, colonyIDString)
	colonyIDList = string.split(colonyIDString, '|')
	
	# map the allele to each colony ID and create lookup
	for c in colonyIDList:
	    if c not in colonyToAlleleDict:
		colonyToAlleleDict[c] = []
		colonyToAlleleDict[c].append(allele)
	    else: # this colony ID is assoc w/>1 allele (or there is a dupe)
		symbolList = []
		alleles =  colonyToAlleleDict[c]
		for a in alleles:
		    symbolList.append(a.asym)
		if alleleSymbol not in symbolList: # don't add duplicate allele to list
		    colonyToAlleleDict[c].append(allele)

    # Query for alleles with colony IDs
    results = db.sql('''select n._Object_key as alleleKey, nc.note
        from MGI_Note n, MGI_NoteChunk nc
        where n._NoteType_key = 1041
        and n._Note_key = nc._Note_key''', 'auto')
    for r in results:
        colonyDict[r['alleleKey']] =  r['note']
 
    # Query for alleles and create lookup
    results = db.sql('''select a._Allele_key, a.symbol as alleleSymbol, 
	    t.term as alleleStatus, t2.term as alleleType, a1.accid as alleleID, 
	    a2.accid as markerID, m.symbol as markerSymbol
	from ALL_Allele a,  ACC_Accession a1, ACC_Accession a2, MRK_Marker m,
	    VOC_Term t, VOC_Term t2
	where a._Allele_Status_key = t._Term_key
	and a._Allele_Type_key = t2._Term_key
	and a._Marker_key = m._Marker_key
	and a._Allele_key = a1._Object_key
	and a1._MGIType_key = 11
	and a1.preferred = 1
	and a1._LogicalDB_key = 1 
	and a._Marker_key = a2._Object_key
	and a2._MGIType_key = 2
	and a2.preferred = 1
	and a2._LogicalDB_key = 1''', 'auto')
    for r in results:
	alleleID  = r['alleleID']
        alleleKey = r['_Allele_key']
	alleleSymbol = r['alleleSymbol']
	alleleStatus = r['alleleStatus']
	alleleType = r['alleleType']
	markerID = r['markerID']
	markerSymbol = r['markerSymbol']
	colonyID = ''
	if alleleKey in colonyDict:
	    colonyID = colonyDict[alleleKey]
	    print 'colony ID from colonyDict: %s' % colonyID
	# create allele object
	allele = Allele(alleleID, alleleSymbol, alleleStatus, alleleType, markerID, markerSymbol, colonyID)
	alleleBySymbolDict[alleleSymbol] = allele
	alleleByIDDict[alleleID] = allele

    # Query for lab codes and create lookup
    results = db.sql('''select term, abbreviation from VOC_Term
	where _Vocab_key = 71''', 'auto')
    for r in results:
	labCodeDict[r['abbreviation']] = r['term']
    
    # Query for markers and create lookup
    results = db.sql('''select a.accid, m.symbol, m.name
	from MRK_Marker m, ACC_Accession a
	where m._Marker_Status_key = 1
	and m._Marker_Type_key in (1, 7)
	and m._Marker_key = a._Object_key
	and a._MGIType_key = 2
	and a._LogicalDB_key = 1
	and a.prefixPart = 'MGI:'
	and a.preferred = 1''', 'auto')
    for r in results:
	markerDict[r['accid']] = '%s|%s' % (r['name'], r['symbol'])

    # Query for strains
    results = db.sql('''select strain from PRB_Strain
	where standard = 1 and private = 0''', 'auto')
    for r in results:
	strainList.append(r['strain'])

    return 0

#
# Purpose: Open files.
#
def openFiles():
    global fpLogDiag, fpLogCur, fpQC
    global fpIMPC, fpAllele, fpNoteload

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

    #
    # Open Colony ID noteload file
    #
    try:
        fpNoteload = open(noteloadFile, 'w')
    except:
        print 'Cannot open file: ' + noteloadFile
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

    if fpNoteload:
	fpNoteload.close()

    return 0

#
# Purpose: query for the MGI Type of an accession id
#
def queryMGIType(id):
    # exclude VOC_Evidence (25)
    results = db.sql('''select am.tableName
	from ACC_Accession a, ACC_MGIType am
	where a.accid = '%s'
	and a._LogicalDB_key = 1
	and a.prefixPart = 'MGI:'
	and a._MGIType_key not in (25) 
	and a._MGIType_key = am._MGIType_key''' % id, 'auto')
    if results == []:
	return ''
    else:
	typeList = []
	for r in results:
	    typeList.append(r['tableName'])
	return string.join(typeList, ', ')

#
# Purpose: find the lab code in an allele subscript
#
def findLabCode(alleleSS):
    labCode = ''
    labCodeFinder = re.compile ('\)(\w*)')
    match = labCodeFinder.search(alleleSS)
    if match:
	labCode = match.group(1)
    return labCode

#
# Purpose: query for an allele by its symbol
#
def findAlleleBySymbol(symbol):
    results = db.sql('''select t.term as status, a.symbol, aa.accid
	from ALL_Allele a, VOC_Term t, ACC_Accession aa
	where a.symbol  = '%s'
	and a._Allele_Status_key = t._Term_key
	and aa._Object_key = a._Allele_key
	and aa._MGIType_key = 11
	and aa._LogicalDB_key = 1
	and aa.preferred = 1
	and aa.prefixPart = 'MGI:' ''' % symbol, 'auto')
    return results
	
#
# Purpose: Read the IMPC file and QC. Create a Allele input file
#
def createAlleleFile():
    global missingRequiredValueList, cidMatchToMultiList, cidMatchMarkerIdMismatchList
    global cidMatchAlleleStatusDiscrepList
    global labCodeNotInMgiList, markerIdNotInMgiList, alleleIdNotInMGIList, strainNotInMgiList
    global unknownAlleleClassList, unknownAlleleTypeList, alleleIdMatchAlleleStatusDiscrepList
    global alleleIdMatchMarkerIdMismatchList, alleleIdMatchAlleleSSMismatchList
    global alleleIdMatchColonyIDMismatchList, alleleIdMatchColonyIdMatchToMultiList
    global alleleIdMatchColonyIdMatchToDiffAlleleList, cidMatchAlleleSSMismatchList
    global SymbolMatchAlleleNotApprovedList, symbolMatchColonyIdMismatchList
    global symbolMatchMultiAlleleList, addCidSymbolMatchList, addCidAlleleIDMatchList

    lineNum = 0
    header = fpIMPC.readline()
    for line in fpIMPC.readlines(): 
	lineNum += 1
	hasError = 0
        tokens = map(string.strip, line[:-1].split('\t'))
	print '#### Split input line: %s' % tokens
	
	# tokens[0] -  marker symbol, not used by the load
        markerID = tokens[1]
	# tokens[2] - es cell line, not used by load
	colonyID = tokens[3]
	strain = string.strip(tokens[4]) # colony background strain
	# tokens[5] - project name, not used by load
	# tokens[6] - production center, not used by load
	alleleClass = tokens[7] # formerly allele type
	alleleType = tokens[8] # formerly mutation type
	#alleleSubType = tokens[9] # temp, will be using this when req done
	alleleDescription = tokens[10]
	alleleSuperScript = tokens[11] # was symbol, now just superscript
	alleleID = tokens[12] # can be blank

	# Translate colony background strain; 3 cases
	if strain == 'C57BL/6NTac/Den':
	    strain = 'C57BL/6NTac'
	elif strain == 'C57BL/6NTac/USA':
	    strain = 'C57BL/6NTac'
	elif strain == 'C57BL6/NCrl':
	    strain = 'C57BL/6NCrl'

	# Requirement 7.2A1 Missing or Rejected Values for Required Fields
	missingDataList = []
	# report missing required values
	if markerID == '':
	    missingDataList.append('Marker ID')
        if colonyID == '':
            missingDataList.append('Colony ID')
        if strain == '':
            missingDataList.append('Strain') 
	if alleleClass == '':
	    missingDataList.append('Allele Class (type)')
	if alleleType == '':
	    missingDataList.append('Allele (mutation) Type')
	if alleleSubType == '':
            missingDataList.append('Allele SubType')
	if alleleDescription == '':
	    missingDataList.append('Allele Description')
	if alleleSuperScript == '':
	    missingDataList.append('Allele Superscript')

	if len(missingDataList):
	    missingRequiredValueList.append('%s%s%s%s%s' % (lineNum, TAB, string.join(missingDataList, ', '), TAB, line))
	    print '  ### missing fields in input file, skip remaining QC'
	    continue	# If missing fields skip remainder of QC

	# Requirement 7.2A1 col2
        if markerID not in markerDict:  
            markerIdNotInMgiList.append('%s%s%s' % (lineNum, TAB, line))
            hasError = 1

	# Requirement 7.2A1 col5
	if strain not in strainList:
	    strainNotInMgiList.append('%s%s%s' % (lineNum, TAB, line))
            hasError = 1

	 # Requirement 7.2A1 col8
	if string.lower(alleleClass) != 'endonuclease-mediated':
	    unknownAlleleClassList.append('%s%s%s' % (lineNum, TAB, line))
            hasError = 1 

	# Requirement 7.2A1 col9
	if string.lower(alleleType) not in ('deletion', 'hdr', 'hr', 'indel'):
	    unknownAlleleTypeList.append('%s%s%s' % (lineNum, TAB, line))
            hasError = 1

	if hasError: # skip to next line if any of the above checks fails
	    print '  ### unexpected data in input file, skip remaining QC'
	    continue

	#
	# get the marker name and symbol and calculate the allele symbol
	#
	marker = markerDict[markerID] # we've checked that markerID is in DB above
	markerName, markerSymbol = string.split(marker, '|')
	#print 'markerName: %s markerSymbol: %s' % (markerName, markerSymbol)

	calcAlleleSymbol = '%s<%s>' % (markerSymbol, alleleSuperScript)

	# BEGIN ALLELE ID PRESENT IN INPUT
	# Requirement  7.2.C
	# if the allele MGI ID is present in the input, check if it is
	# 1. in the database
	# 2. if so, is it an approved allele
	# 3. if not, is it another object type
	if alleleID != '': # Allele ID present Flow Diagram Box B
	    print '  #### Allele ID Present: %s' % alleleID
	    if alleleID in alleleByIDDict:    # Requirement 7.2.D Allele ID in MGI check
		dbA = alleleByIDDict[alleleID]
		print 'dbA.asym: %s' % dbA.asym
		# if not 'Approved', don't do any other checks.
		if dbA.ast != 'Approved':    # Requirement 7.2.D3 Allele ID status check
		    alleleIdMatchAlleleStatusDiscrepList.append('%s%s%s%s%s' % \
			(lineNum, TAB, dbA.ast, TAB, line))
                    hasError = 1
		else:   
		    # Requirement 7.2.D1 Marker ID check
		    if markerID != dbA.mid:
		        alleleIdMatchMarkerIdMismatchList.append( '%s%s%s%s%s%s%s%s%s%s%s' % \
			    (lineNum, TAB, alleleID, TAB, calcAlleleSymbol, TAB, dbA.mid, TAB, dbA.ms, TAB, line))
			hasError = 1
		    # Requirement 7.2.D2 Allele symbol check
		    print 'alleleSuperScript: %s' % alleleSuperScript
		    print 'dbAlleleSymbol: %s' % dbA.asym
		    if string.find(dbA.asym, alleleSuperScript) == -1:
			alleleIdMatchAlleleSSMismatchList.append('%s%s%s%s%s%s%s' % \
			    (lineNum, TAB, dbA.aid, TAB, dbA.asym, TAB, line))
			hasError = 1
		    # Requirement 7.2.D4 Colony Name/ID check
		    # From the set of cid(s) (0..n) associated with allele ID in the 
		    # db, incoming cid must be in this set (if any exist)
		    # if incoming cid not in the set - add it to the db
		    # Additional cids in db is OK as long as input cid in db
		    
		    print 'dbA.cid: %s' % dbA.cid
		    dbColonyIDList = []
		    cidError = 0 # assume there's no error
		    if dbA.cid != '':
			dbColonyIDList.append(dbA.cid)
		    print 'dbColonyIDList: %s' % dbColonyIDList
		    print 'IncColonyID: %s' % colonyID

		    # Requirement 7.2.D4a Allele ID match, Colony ID Mismatch
		    if dbColonyIDList != [] and colonyID not in dbColonyIDList:
			alleleIdMatchColonyIDMismatchList.append('%s%s%s%s%s%s%s%s%s' % (lineNum, TAB, alleleID, TAB, dbA.asym, TAB, dbA.cid, TAB, line))
			hasError = 1
			cidError = 1

		    # Requirement 7.2.D4b Colony ID matches MULTIPLE  alleles in the database
		    # get the set of allele(s) (0..n) associated with incoming cid
		    # if there are multiple alleles in the set report
		    if colonyID in colonyToAlleleDict:
			print '%s in colonyToAlleleDict' % colonyID
			allelesByCidList = colonyToAlleleDict[colonyID]
			if len(allelesByCidList) > 1:
			    for aByCid in allelesByCidList:
				alleleIdMatchColonyIdMatchToMultiList.append('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (lineNum, TAB, alleleID, TAB, dbA.asym, TAB, dbA.cid, TAB, aByCid.aid, TAB, aByCid.asym, TAB, aByCid.at, TAB, aByCid.cid, TAB, line))
				hasError = 1
				cidError = 1
			else: # 7.2.D4b  Colony ID matches SINGLE allele in the database
			    aByCid = allelesByCidList[0]
			    if alleleID != aByCid.aid:
				alleleIdMatchColonyIdMatchToDiffAlleleList.append('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (lineNum, TAB, alleleID, TAB, dbA.asym, TAB, dbA.cid, TAB, aByCid.aid, TAB, aByCid.asym, TAB, aByCid.at, TAB, aByCid.cid, TAB, line))
				hasError = 1
				cidError = 1

		    if hasError == 0 and cidError == 0 and dbColonyIDList == []:
			# Requirement 7.2.D4 if no error and no cid in the database, add a 
			# new note to the allele
		    	fpNoteload.write('%s%s%s%s' % (alleleID, TAB, colonyID, CRT))
		    	addCidAlleleIDMatchList.append('%s%s%s%s%s%s%s%s%s%s%s' % \
	                       (lineNum, TAB, alleleID, TAB, calcAlleleSymbol, TAB, colonyID, TAB, colonyID, TAB, line))


	    else: # Requirement 7.2.C1 Allele ID not in MGI OR matches different object type
		print 'Allele ID not in MGI OR matches a different object type'
		objectType = queryMGIType(alleleID)
		# report: 
		# error type: 'MGI Allele Accession present, No MGI Allele Match
		# if different MGI Type, report this line from input
		alleleIdNotInMGIList.append('%s%s%s%s%s' % \
		    (lineNum, TAB, objectType, TAB, line))
		hasError = 1
	    # END ALLELE ID PRESENT IN INPUT

	# BEGIN ALLELE ID NOT PRESENT IN INPUT
	else: 

	    # BEGIN COLONY ID MATCH 7.2.E
	    print '  #### Allele ID not in input, colonyID is: %s' % colonyID
	    if colonyID in colonyToAlleleDict: 
	        alleleList = colonyToAlleleDict[colonyID]
		# Requirement 7.2.F1 Colony ID Matches Multiple Alleles in MGI
		if len(alleleList) > 1:
		    for dbA in alleleList:
			# report multiple alleles for a colony ID
			cidMatchToMultiList.append('%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (lineNum, TAB, dbA.aid, TAB, dbA.asym, TAB, line, CRT))
			hasError = 1
			print '  ###  multiple alleles for colony ID, skip remaining checks'
			continue # multiple alleles for colony ID, don't do further checks

		dbA = alleleList[0] # there is only one
		# Requirement 7.2.F3 allele Status Check
		if dbA.ast != 'Approved':  
                    cidMatchAlleleStatusDiscrepList.append('%s%s%s%s%s%s%s%s%s' % \
                        (lineNum, TAB,  dbA.aid, TAB, dbA.asym, TAB, dbA.ast, TAB, line))
                    hasError = 1
		else:
		    # The following two checks could be replaced with a calculated allele symbol match
		    # Requirement 7.2.F2 Marker ID check
		    if markerID != dbA.mid:
			cidMatchMarkerIdMismatchList.append( '%s%s%s%s%s%s%s%s%s%s%s' % \
			    (lineNum, TAB, dbA.aid, TAB, dbA.asym, TAB, line))
			hasError = 1
		    # Requirement 7.2.F2 Allele superscript check
		    print ' cid match, alleleSuperScript: %s' % alleleSuperScript
		    print 'cid match dbAlleleSymbol: %s' % dbA.asym
		    if string.find(dbA.asym, alleleSuperScript) == -1:
			cidMatchAlleleSSMismatchList.append('%s%s%s%s%s%s%s' % \
			    (lineNum, TAB, dbA.aid, TAB, dbA.asym, TAB, line))
			hasError = 1
		# END COLONY ID MATCH

	    # BEGIN NO COLONY ID MATCH Requirement 7.2.G
	    else:
		print '  #### No allele ID, no cid match, calculate allele symbol and check in DB'
		print 'calcAlleleSymbol: %s' % calcAlleleSymbol

		symbolError = 0 # default to no error

		# if match found in database, these will be assigned.
		aID = '' 
		symbol = ''

		# Requirement 7.2.H Allele symbol check
		results = findAlleleBySymbol(calcAlleleSymbol)

		# Requirement 7.2.H No CID Match, Allele Symbol Match
		#elif len(results) == 1:
		if len(results) == 1:
		    print 'no cid match, but symbol match: %s' % results
		    status = results[0]['status']
		    aID = results[0]['accid']
		    symbol = results[0]['symbol']

		    # Requirement 7.2.H1  Allele Status Check
		    if status != 'Approved':
			symbolMatchAlleleStatusDiscrepList.append('%s%s%s%s%s%s%s%s%s' % \
			    (lineNum, TAB, aID, TAB, symbol, TAB, status, TAB, line))
			symbolError = 1
			hasError = 1

		    # Requirement 7.2.H2 Colony Name/ID Check
		    if aID in alleleByIDDict: # has to be, but good to check
			allele = alleleByIDDict[aID]
			if allele.cid != '': # if there is a cid for the symbol matched allele
					     # it has to be a mismatch with the inc cid
			    symbolMatchColonyIdMismatchList.append('%s%s%s%s%s%s%s%s%s' % \
				(lineNum, TAB, aID, TAB, symbol, TAB, allele.cid, TAB, line))
			    symbolError = 1
			    hasError = 1
		    # Requirement 7.2.H4 No CID Match, Symbol match, and no errors
		    # add new note to the allele
		    if symbolError == 0 and hasError == 0:
			fpNoteload.write('%s%s%s%s' % (aID, TAB, colonyID, CRT))
			addCidSymbolMatchList.append('%s%s%s%s%s%s%s%s%s%s%s' % \
			    (lineNum, TAB, aID, TAB, symbol, TAB, colonyID, TAB, colonyID, TAB, line))

		# Requirement 7.2.H3 check for multiple (duplicate) alleles in the database
		#else: # len(results) > 1:
		elif len(results) > 1:
		    print 'no cid match, symbol match to dupe alleles in database: %s' % results
 		    for r in results:
			print r
		    	aID = r['accid']
			symbol = r['symbol']
			symbolMatchMultiAlleleList.append('%s%s%s%s%s%s%s' % \
                            (lineNum, TAB, aID, TAB, symbol, TAB, line))
			symbolError = 1
                        hasError = 1

		# else no results, so create allele if labCode in DB

	# Requirement 7.2.I
	labName = ''
	labCode = findLabCode(alleleSuperScript)
	print '  #### checking lab code: %s' % labCode
	if labCode not in labCodeDict:
	    labCodeNotInMgiList.append('%s%s%s%s%s' % (lineNum, TAB, labCode, TAB, line))
	    print '  #### bad lab code, not createing allele'
	    hasError = 1
	
	#
	# If no errors write out to allele file
	#
	if hasError == 0:
	    print '  #### No allele identified in DB and no errors; translate stuff and create allele'
	    # translate allele type
	    if string.lower(alleleType) in ('deletion', 'hdr', 'hr'):
		alleleType = 'intragenic deletion' 
	    else:
		alleleType = 'intragenic deletion|insertion'

	    # get the sequencNum from the allele
	    seqNumFinder = re.compile ( 'em(.*)\(' )
	    match = seqNumFinder.search(alleleSuperScript)
	    sequenceNum = match.group(1)

	    # get the lab name from the lab code
	    labName = labCodeDict[labCode]

	    # calculate allele name
	    alleleName = alleleNameTemplate % (markerName, sequenceNum, labName)

	    fpAllele.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (markerID, TAB, alleleType, TAB, alleleDescription, TAB, colonyID, TAB, strain, TAB, calcAlleleSymbol, TAB, alleleName, TAB, inHeritMode, TAB, alleleClass, TAB, alleleSubType, TAB, alleleStatus, TAB, transmissionState, TAB, alleleCollection, TAB, jNumber, TAB, createdBy, CRT))

    return 0

def writeQCReport():
    fpQC.write('7.2.A.1 Required Value Missing or Invalid%s%s' % (CRT, CRT))
    fpQC.write('Line#%s Missing Value(s)%sInput Line%s' % (TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(missingRequiredValueList):
	fpQC.write(string.join(missingRequiredValueList, CRT))
    fpQC.write('Total: %s' % len(missingRequiredValueList))

    fpQC.write('%s%s7.2.C1 MGI Allele ID present, No MGI Allele Match%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Line%s' % (TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(markerIdNotInMgiList):
	 fpQC.write(string.join(markerIdNotInMgiList))
    fpQC.write('Total: %s' % len(markerIdNotInMgiList))

    fpQC.write('%s%s7.2.A1 Colony Background Strain not in MGI%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Line%s' % (TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(strainNotInMgiList):
	 fpQC.write(string.join(strainNotInMgiList))
    fpQC.write('Total: %s' % len(strainNotInMgiList))

    fpQC.write('%s%s7.2.A1 Allele Class not endonuclease-mediated%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Line%s' % (TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(unknownAlleleClassList):
         fpQC.write(string.join(unknownAlleleClassList))
    fpQC.write('Total: %s' % len(unknownAlleleClassList))
   
    fpQC.write('%s%s7.2.A1 Allele (mutation) Type not in Translated Set%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Line%s' % (TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(unknownAlleleTypeList):
         fpQC.write(string.join(unknownAlleleTypeList))
    fpQC.write('Total: %s' % len(unknownAlleleTypeList)) 

    fpQC.write('%s%s7.2.C1 MGI Allele ID present, No MGI Allele Match%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sObjectType%sInput Line%s' % (TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdNotInMGIList):
         fpQC.write(string.join(alleleIdNotInMGIList))
    fpQC.write('Total: %s' % len(alleleIdNotInMGIList))

    fpQC.write('%s%s7.2.D3 Allele ID Match, Allele Status Discrepancy"%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sAllele Status%sInput Line%s' % (TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdMatchAlleleStatusDiscrepList):
         fpQC.write(string.join(alleleIdMatchAlleleStatusDiscrepList))
    fpQC.write('Total: %s' % len(alleleIdMatchAlleleStatusDiscrepList))

    fpQC.write('%s%s7.2.D1 Allele ID Match, Marker ID Mismatch%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sAllele ID%sAllele Symbol%sDB Marker ID%sDB Marker Symbol%sInput Line%s' % (TAB, TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdMatchMarkerIdMismatchList):
         fpQC.write(string.join(alleleIdMatchMarkerIdMismatchList))
    fpQC.write('Total: %s' % len(alleleIdMatchMarkerIdMismatchList))

    fpQC.write('%s%s7.2.D2 Allele ID Match, Allele Symbol Superscript Mismatch%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sInput Line%s' % (TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdMatchAlleleSSMismatchList):
         fpQC.write(string.join(alleleIdMatchAlleleSSMismatchList))
    fpQC.write('Total: %s' % len(alleleIdMatchAlleleSSMismatchList))

    fpQC.write('%s%s7.2.D4a Allele ID Match, Colony ID Mismatch%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sAllele ID%sDB Allele Symbol%sDB CID%sInput Line%s' % (TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdMatchColonyIDMismatchList):
         fpQC.write(string.join(alleleIdMatchColonyIDMismatchList))
    fpQC.write('Total: %s' % len(alleleIdMatchColonyIDMismatchList))

    fpQC.write('%s%s7.2.D4b Allele ID match, Colony ID Match to Multi MGI Alleles%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Allele ID%sInput Allele Symbol%sInput Colony ID%sDB Allele ID%sDB Allele Symbol%sDB Allele Type%sDB Colony ID%sInput Line%s' % (TAB, TAB, TAB, TAB, TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdMatchColonyIdMatchToMultiList):
         fpQC.write(string.join(alleleIdMatchColonyIdMatchToMultiList))
    fpQC.write('Total: %s' % len(alleleIdMatchColonyIdMatchToMultiList))

    fpQC.write('%s%s7.2.D4b Allele ID match, Colony ID Match to Different Allele%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Allele ID%sInput Allele Symbol%sInput Colony ID%sDB Allele ID%sDB Allele Symbol%sDB Allele Type%sDB Colony ID%sInput Line%s' % (TAB, TAB, TAB, TAB, TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdMatchColonyIdMatchToDiffAlleleList):
         fpQC.write(string.join(alleleIdMatchColonyIdMatchToDiffAlleleList))
    fpQC.write('Total: %s' % len(alleleIdMatchColonyIdMatchToDiffAlleleList))

    fpQC.write('%s%s7.2.F1 Colony ID Matches Multiple Alleles%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sInput Line%s' % (TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(cidMatchToMultiList):
         fpQC.write(string.join(cidMatchToMultiList))
    fpQC.write('Total: %s' % len(cidMatchToMultiList))

    fpQC.write('%s%s7.2.F2a Colony ID Match, Marker ID Mismatch%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sInput Line%s' % (TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(cidMatchMarkerIdMismatchList):
         fpQC.write(string.join(cidMatchMarkerIdMismatchList))
    fpQC.write('Total: %s' % len(cidMatchMarkerIdMismatchList))

    fpQC.write('%s%s7.2.F2b Colony ID Match, Allele Symbol Superscript Mismatch%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sInput Line%s' % (TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(cidMatchAlleleSSMismatchList):
         fpQC.write(string.join(cidMatchAlleleSSMismatchList))
    fpQC.write('Total: %s' % len(cidMatchAlleleSSMismatchList))

    fpQC.write('%s%s7.2.F3 Colony ID Match, Allele Status Discrepancy%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sDB Allele Status%sInput Line%s' % (TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(cidMatchAlleleStatusDiscrepList):
         fpQC.write(string.join(cidMatchAlleleStatusDiscrepList))
    fpQC.write('Total: %s' % len(cidMatchAlleleStatusDiscrepList))

    fpQC.write('%s%s7.2.H1 Allele Symbol Match, Allele Status Discrepancy%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sDB Allele Status%sInput Line%s' % (TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(symbolMatchAlleleStatusDiscrepList):
         fpQC.write(string.join(symbolMatchAlleleStatusDiscrepList))
    fpQC.write('Total: %s' % len(symbolMatchAlleleStatusDiscrepList))

    fpQC.write('%s%s7.2.H2 Allele Symbol Match, Colony ID Mismatch%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sDB Allele CID%sInput Line%s' % (TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(symbolMatchColonyIdMismatchList):
         fpQC.write(string.join(symbolMatchColonyIdMismatchList))
    fpQC.write('Total: %s' % len(symbolMatchColonyIdMismatchList))

    fpQC.write('%s%s7.2.H3 Allele Symbol Match to Multiple Alleles%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sInput Line%s' % (TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(symbolMatchMultiAlleleList):
         fpQC.write(string.join(symbolMatchMultiAlleleList))
    fpQC.write('Total: %s' % len(symbolMatchMultiAlleleList))

    fpQC.write('%s%s7.2.I No Allele Match, Lab Code not Present%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sLab Code%sInput Line%s' % (TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(labCodeNotInMgiList):
        fpQC.write(string.join(labCodeNotInMgiList, CRT))
    fpQC.write('Total: %s' % len(labCodeNotInMgiList))

    fpQC.write('%s%s7.2.H4 CID added to Allele where Allele ID was Matched%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sAllele ID%sAllele Symbol%sColony ID%sInput Line%s' % (TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(addCidAlleleIDMatchList):
         fpQC.write(string.join(addCidAlleleIDMatchList))
    fpQC.write('Total: %s' % len(addCidAlleleIDMatchList))

    fpQC.write('%s%s7.2.H4 CID added to Allele where Symbol was Matched%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sAllele ID%sAllele Symbol%sColony ID%sInput Line%s' % (TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(addCidSymbolMatchList):
         fpQC.write(string.join(addCidSymbolMatchList))
    fpQC.write('Total: %s' % len(addCidSymbolMatchList))

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

