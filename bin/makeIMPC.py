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
#      IMPC file ($SOURCE_COPY_INPUT_FILE) - the latest 9 col version from GenTar
#           text is from the column header       
#       field 1: Gene Symbol
#       field 2: Gene MGI Accession ID
#       field 3: Colony Name
#       field 4: Colony Background Strain
#       field 5: Mutation Class
#       field 6: Mutation Type
#       field 7: Mutation Subtype
#       field 8: Mutation Symbol
#       field 9: Mutation MGI Accession ID

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
#       01/25/2022    sc
#        - wts2-767 remove mgi_notechunk, note now in mgi_note
#
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
host = ''
alleleDescription = '' # gentar project - do not create molecular note

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
#

# expected values for impc allele type (mutation type)
impcAlleleTypeList = []

# expected values for impc allele subtype
impcSubTypeList = []

# allele type/subtype translation
alleleTypeTransDict = {}

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

# regex to find labcode in symbol
labCodeFinder = re.compile ('\)(\w+)')

# regex to check for multiple < or > in symbol
pattern1 = r'<'
pattern2 = r'>'

# marker ID to marker name lookup (for Allele Name construction)
# {markerID: markerName|symbol, ...}
markerDict = {}

# {alleleKey:colonyIDNote, ...}
colonyDict = {}

# list of standard and non-private strains in the database
strainList = []

# template for creating allelel name for new alleles
# marker name, sequenceNum, lab code name
alleleNameTemplate = 'endonuclease-mediated mutation %s, %s'

# Total number of input lines skipped because of discrepancies
linesLoadedCt = 0

# Total number of input lines loaded
linesSkippedCt = 0

# Total number of alleles found to be in the DB
allelesFoundCt = 0

# current line number - used for Total count in reporting
lineNum = 0

# current set of alleles seen in the input, this is the calculated symbol
# using marker symbol and allele subscript from the input
# key = calculated symbol, value = list of lists
# each member of list represents one line from input
# inner list as three parts
# 1. alleleFile line, 2. input line number, 3. input line
calcAlleleDict = {}

#
# QC lists for reporting errors
#
missingRequiredValueList = []	
labCodeNotInMgiList = []	
badNomenList = []
markerIdNotInMgiList = []	
strainNotInMgiList = []		
unknownAlleleClassList = []     
unknownAlleleTypeList = []      
unknownSubTypeList = []
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
dupeAlleleInInputList = []
atTransKeyNotInMgiList = []

class Allele:
    #
    # Is: data object for a Allele
    # Has: a set of allele attributes
    # Does: provides direct access to its attributes
    #
    def __init__(self, alleleID,    # str.- allele  MGI ID
            alleleSymbol,           # str.- allele symbol
            alleleStatus,	    # str.- allele status
            alleleType,		    # str.- allele type
            markerID,               # str.- marker MGI ID
            markerSymbol,	    # str.- marker symbol
            markerKey,              # integer - marker primary key
            colonyID):		    # str.- pipe delim colony ID string
        self.aid = alleleID
        self.asym = alleleSymbol
        self.ast = alleleStatus
        self.at = alleleType
        self.mid = markerID
        self.ms = markerSymbol
        self.mk = markerKey
        self.cid = colonyID
    def toString(this):
        return '%s, %s, %s, %s, %s, %s, %s' % (this.aid, this.asym, this.ast, this.mid, this.ms, this.mk, this.cid)

def initialize():
    # Purpose: create lookups, open files
    #   get max keys from the db
    # Returns: 1 if error, else 0
    # Assumes: Nothing
    # Effects: Sets global variables, exits if a file can't be opened,
    #  creates files in the file system

    global logDiagFile, logCurFile, qcFile, impcFile, alleleFile, noteloadFile
    global jNumber, createdBy, inHeritMode, alleleStatus
    global transmissionState, alleleCollection, strainList
    global colonyToAlleleDict, alleleBySymbolDict, labCodeDict, markerDict
    global colonyDict, host, alleleTypeTransDict, impcAlleleTypeList
    global impcSubTypeList, calcAlleleDict

    db.useOneConnection(1)

    logDiagFile = os.getenv('LOG_DIAG')
    logCurFile = os.getenv('LOG_CUR')
    qcFile = os.getenv('QC_FILE')
    impcFile = os.getenv('SOURCE_COPY_INPUT_FILE')
    alleleFile = os.getenv('ALLELE_FILE')
    noteloadFile = os.getenv('CID_NOTE_FILE')
    jNumber = os.getenv('JNUMBER')
    createdBy = os.getenv('CREATEDBY')
    inHeritMode = os.getenv('INHERIT_MODE')
    alleleStatus = os.getenv('ALLELE_STATUS')
    transmissionState = os.getenv('TRANSMISSION_STATE')
    alleleCollection = os.getenv('ALLELE_COLLECTION')
    host = os.getenv('HOST')
    
    impcAlleleTypeList = str.split(os.getenv('IMPC_ALLELETYPES'), '|')
    impcSubTypeList = str.split(os.getenv('IMPC_SUBTYPES'), '|')
    #print('impcAlleleTypeList: %s' % impcAlleleTypeList)
    #print('impcSubTypeList: %s' % impcSubTypeList)

    alleleTypeTransString = os.getenv('ALLELE_TYPE_TRANS')
    alleleTypeTransDict = dict(x.split('=') for x in alleleTypeTransString.split('\n'))
    #print('alleleTypeTransString: %s' % alleleTypeTransString)
    #print('alleleTypeTransDict: %s' % alleleTypeTransDict)

    if openFiles() != 0:
        sys.exit(1)

    # Query for IKMC Allele Colony Name - there are multi per allele
    results = db.sql('''select distinct n.note as cidNote, 
            a.symbol as alleleSymbol, t.term as alleleStatus, 
            t2.term as alleleType, m.symbol as markerSymbol, m._Marker_key,
            a1.accid as alleleID, a2.accid as markerID, 
            a1.preferred as allelePref, a2.preferred as markerPref
        from MGI_Note n, ALL_Allele a, MRK_Marker m, ACC_Accession a1, ACC_Accession a2, 
            VOC_Term t, VOC_Term t2
        where n._NoteType_key = 1041
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
        colonyIDString = str.strip(r['cidNote'])
        alleleSymbol = r['alleleSymbol']
        alleleStatus = r['alleleStatus']
        alleleType = r['alleleType']
        alleleID  = r['alleleID']
        markerSymbol = r['markerSymbol']
        markerKey = r['_Marker_key']
        markerID = r['markerID']
        # create allele object
        allele = Allele(alleleID, alleleSymbol, alleleStatus, alleleType, markerID, markerSymbol, markerKey, colonyIDString)
        colonyIDList = str.split(colonyIDString, '|')
        
        # map the allele to each colony ID and create lookup
        for c in colonyIDList:
            cLower = str.lower(c)
            if cLower not in colonyToAlleleDict:
                colonyToAlleleDict[cLower] = []
                colonyToAlleleDict[cLower].append(allele)
            else: # this colony ID is assoc w/>1 allele (or there is a dupe)
                symbolList = []
                alleles =  colonyToAlleleDict[cLower]
                for a in alleles:
                    symbolList.append(a.asym)
                if alleleSymbol not in symbolList: # don't add duplicate alleles
                    colonyToAlleleDict[cLower].append(allele)

    # Query for alleles with colony IDs
    results = db.sql('''select n._Object_key as alleleKey, n.note
        from MGI_Note n
        where n._NoteType_key = 1041''', 'auto')
    for r in results:
        colonyDict[r['alleleKey']] =  r['note']
 
    # Query for alleles and create lookup
    results = db.sql('''select a._Allele_key, a.symbol as alleleSymbol, 
            t.term as alleleStatus, t2.term as alleleType, a1.accid as alleleID, 
            a2.accid as markerID, m.symbol as markerSymbol, m._Marker_key
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
        markerKey = r['_Marker_key']
        colonyID = ''
        if alleleKey in colonyDict:
            colonyID = colonyDict[alleleKey]
        # create allele object
        allele = Allele(alleleID, alleleSymbol, alleleStatus, alleleType, markerID, markerSymbol, markerKey, colonyID)
        alleleBySymbolDict[alleleSymbol] = allele
        alleleByIDDict[alleleID] = allele

    # Query for lab codes and create lookup
    results = db.sql('''select term, abbreviation from VOC_Term
        where _Vocab_key = 71''', 'auto')
    for r in results:
        labCodeDict[r['abbreviation']] = str.strip(r['term'])
    
    # Query for markers and create lookup
    results = db.sql('''select a.accid, m.symbol, m.name
        from MRK_Marker m, ACC_Accession a
        where m._Marker_Status_key = 1
        and m._Marker_Type_key in (1, 7)
        and m._Marker_key = a._Object_key
        and a._MGIType_key = 2
        and a._LogicalDB_key = 1
        and a.prefixPart = 'MGI:' ''', 'auto')
    for r in results:
        markerDict[r['accid']] = '%s|%s' % (r['name'], r['symbol'])

    # Query for strains
    results = db.sql('''select strain from PRB_Strain
        where private = 0''', 'auto')
    for r in results:
        strainList.append(r['strain'])

    return 0

def openFiles():
    # Purpose: Open input/output files.
    # Returns: 1 if error, else 0
    # Assumes: Nothing
    # Effects: Sets global variables, exits if a file can't be opened,
    #  creates files in the file system

    global fpLogDiag, fpLogCur, fpQC
    global fpIMPC, fpAllele, fpNoteload

    #
    # Open the Log Diag file; append to existing file
    #
    try:
        fpLogDiag = open(logDiagFile, 'a+')
    except:
        print('Cannot open file: ' + logDiagFile)
        return 1

    #
    # Open the Log Cur file; append to existing file
    #
    try:
        fpLogCur = open(logCurFile, 'a+')
    except:
        print('Cannot open file: ' + logCurFile)
        return 1

    #
    # Open the QC file
    #
    try:
        fpQC = open(qcFile, 'w')
    except:
        print('Cannot open file: ' + qcFile)
        return 1

    #
    # Open the IMPC file
    #
    try:
        fpIMPC = open(impcFile, 'r')
    except:
        print('Cannot open file: ' + impcFile)
        return 1

    #
    # Open the IMPC file with genotype sequence #
    #
    try:
        fpAllele = open(alleleFile, 'w')
    except:
        print('Cannot open file: ' + alleleFile)
        return 1

    #
    # Open Colony ID noteload file
    #
    try:
        fpNoteload = open(noteloadFile, 'w')
    except:
        print('Cannot open file: ' + noteloadFile)
        return 1

    return 0


#
# Purpose: Close files.
#
def closeFiles():
    # Purpose: Close all file descriptors
    # Returns: 1 if error, else 0
    # Assumes: all file descriptors were initialized
    # Effects: Nothing
    # Throws: Nothing

    try:
        fpLogDiag.close()
        fpLogCur.close()
        fpQC.close()
        fpIMPC.close()
        fpAllele.close()
        fpNoteload.close()
    except:
        return 1
    return 0

def queryMGIType(id): # An MGI Accession ID
    # Purpose: Find the MGI Type of an MGI Accession ID
    # Returns: '' if "id" is not in database; else tableName
    #	of the MGI Type
    # Assumes:  connection to a database
    # Effects: Nothing
    # Throws: Nothing

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
        return str.join(', ', typeList)

def findLabCode(allele): # an IMPC allele symbol
    # Purpose: Finds the labcode in an allele subscript
    # Returns: '' if no labCode found in "allele; else the labCode
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    # match ')' literally then get group where token(s) in 
    # [a-zA-Z0-9_] (\w) and at least one token (+)
    labCode = ''
    match = labCodeFinder.search(allele)
    if match:
        labCode = match.group(1)
    return labCode

def findAlleleBySymbol(symbol): # an MGI allele symbol
    # Purpose: query for an allele by its symbol
    # Returns: the result set from the query; may be empty
    # Assumes:  db connection
    # Effects: Nothing
    # Throws: Nothing

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
        
def createAlleleFile():
    # Purpose: Read the IMPC file and QC. Create a Allele input file
    # Returns: 1 if error,  else 0
    # Assumes: file descriptors have been initialized
    # Effects: writes to the file system
    # Throws: Nothing

    global missingRequiredValueList, cidMatchToMultiList
    global cidMatchMarkerIdMismatchList, cidMatchAlleleStatusDiscrepList
    global cidMatchAlleleSSMismatchList, badNomenList
    global labCodeNotInMgiList, markerIdNotInMgiList, alleleIdNotInMGIList
    global strainNotInMgiList, unknownAlleleClassList, unknownAlleleTypeList
    global unknonhwnSubTypeList, alleleIdMatchAlleleStatusDiscrepList
    global alleleIdMatchMarkerIdMismatchList, alleleIdMatchAlleleSSMismatchList
    global alleleIdMatchColonyIDMismatchList 
    global alleleIdMatchColonyIdMatchToMultiList
    global alleleIdMatchColonyIdMatchToDiffAlleleList
    global symbolMatchAlleleStatusDiscrepList, symbolMatchColonyIdMismatchList
    global symbolMatchMultiAlleleList, calcAlleleDict, atTransKeyNotInMgiList
    global linesSkippedCt, linesLoadedCt, allelesFoundCt, lineNum

    header = fpIMPC.readline()
    lineNum = 1 # ignoring header
    for line in fpIMPC.readlines(): 
        lineNum += 1
        hasError = 0
        alleleFound = 0
        tokens = list(map(str.strip, line[:-1].split('\t')))
        #print('#### Split input line: %s' % tokens)
        
        # tokens[0] -  marker symbol, not used by the load
        markerID = tokens[1]
        colonyID = tokens[2]
        strain = str.strip(tokens[3]) # colony background strain
        alleleClass = tokens[4] 
        alleleType = tokens[5] 
        alleleSubType = tokens[6] 
        alleleSymbol = tokens[7] # full symbol, was just superscript
        alleleID = tokens[8] # can be blank, if present allele has already been created
        # Translate colony background strain; 2 cases
        if strain == 'C57BL/6NTac/Den':
            strain = 'C57BL/6NTac'
        elif strain == 'C57BL/6NTac/USA':
            strain = 'C57BL/6NTac'

        # check if in the database, if not load allele with Not Specified strain
        # but still report 11/8/22
        if strain not in strainList:
            strainNotInMgiList.append('%s%s%s' % (lineNum, TAB, line))
            strain = 'Not Specified'

        # Requirement 7.2A1 Missing or Rejected Values for Required Fields
        missingDataList = []
        # report missing required values
        if markerID == '':
            missingDataList.append('Marker ID')
        if colonyID == '':
            missingDataList.append('Colony ID')
        #if strain == '': # 11/8/22  strain now set above
        #    missingDataList.append('Strain') 
        if alleleClass == '':
            missingDataList.append('Allele Class (type)')
        if alleleType == '':
            missingDataList.append('Allele (mutation) Type')
        # empty alleleDescription now means to not create molecular note
        #if alleleDescription == '': 11/22 not loading molecular note now
        #    missingDataList.append('Allele Description')
        if alleleSymbol == '':
            missingDataList.append('Allele Symbol')

        if len(missingDataList):
            missingRequiredValueList.append('%s%s%s%s%s' % (lineNum, TAB, str.join(', ', missingDataList), TAB, line))
            #print('  ### missing fields in input file, skip remaining QC')
            linesSkippedCt += 1
            continue	# If missing fields skip remainder of QC

        # Requirement 7.2A1 col2
        if markerID not in markerDict:  
            markerIdNotInMgiList.append('%s%s%s' % (lineNum, TAB, line))
            hasError = 1

         # Requirement 7.2A1 col8
        if str.lower(alleleClass) != 'endonuclease-mediated':
            unknownAlleleClassList.append('%s%s%s' % (lineNum, TAB, line))
            hasError = 1 
        else:
            alleleClass = 'Endonuclease-mediated' # not capitalized in the file, cap in DB

        # Requirement 7.2A1 col9
        if str.lower(alleleType) not in impcAlleleTypeList:
            unknownAlleleTypeList.append('%s%s%s' % (lineNum, TAB, line))
            hasError = 1
        # Requirement 7.2A1 col10
        if alleleSubType != '' and str.lower(alleleSubType) not in impcSubTypeList:
            unknownSubTypeList.append('%s%s%s' % (lineNum, TAB, line))
            hasError = 1
        if hasError: # skip to next line if any of the above checks fails
            #print('  ### unexpected data in input file, skip remaining QC')
            linesSkippedCt += 1
            continue

        #
        # get the marker name and symbol and calculate the allele symbol
        #
        marker = markerDict[markerID] # we've checked that markerID is in DB above
        markerName, markerSymbol = str.split(marker, '|')

        calcAlleleSymbol = alleleSymbol # '%s<%s>' % (markerSymbol, alleleSymbol)

        # BEGIN ALLELE ID PRESENT IN INPUT
        # Requirement  7.2.C
        # if the allele MGI ID is present in the input, check if it is
        # 1. in the database
        # 2. if so, is it an approved allele
        # 3. if not, is it another object type
        if alleleID != '': # Allele ID present Flow Diagram Box B
            #print('  #### Allele ID Present: %s' % alleleID)
            if alleleID in alleleByIDDict:    # Requirement 7.2.D Allele ID in MGI check
                dbA = alleleByIDDict[alleleID]
                #print('dbA.asym: %s' % dbA.asym)
                # if not 'Approved', don't do any other checks.
                if dbA.ast != 'Approved':    # Requirement 7.2.D3 Allele ID status check
                    alleleIdMatchAlleleStatusDiscrepList.append('%s%s%s%s%s' % \
                        (lineNum, TAB, dbA.ast, TAB, line))
                    hasError = 1
                else:   
                    # Requirement 7.2.D1 Marker ID check, 2ndary OK
                    if markerID != dbA.mid: # check to see if 2ndary
                        #print('markerID: %s dbID: %s' % (markerID, dbA.mid))
                        results = db.sql('''select accid
                                from ACC_Accession
                                where _MGIType_key = 2
                                and _LogicalDB_key = 1
                                and _Object_key = %s ''' % dbA.mk, 'auto')
                        isSecondary = 0
                        for r in results:
                            #print('dbAccID: %s' %  r['accid'])
                            if markerID == r['accid']:
                                isSecondary = 1
                                break
                        if not isSecondary:
                            alleleIdMatchMarkerIdMismatchList.append( '%s%s%s%s%s%s%s%s%s%s%s' % \
                                (lineNum, TAB, alleleID, TAB, calcAlleleSymbol, TAB, dbA.mid, TAB, dbA.ms, TAB, line))
                            hasError = 1
                    # Requirement 7.2.D2 Allele symbol check
                    #print('alleleSymbol: %s' % alleleSymbol)
                    #print('dbAlleleSymbol: %s' % dbA.asym)
                    if str.find(dbA.asym, alleleSymbol) == -1:
                        alleleIdMatchAlleleSSMismatchList.append('%s%s%s%s%s%s%s' % \
                            (lineNum, TAB, dbA.aid, TAB, dbA.asym, TAB, line))
                        hasError = 1
                    # Requirement 7.2.D4 Colony Name/ID check
                    # From the set of cid(s) (0..n) associated with allele ID in the 
                    # db, incoming cid must be in this set (if any exist)
                    # if incoming cid not in the set - add it to the db
                    # Additional cids in db is OK as long as input cid in db
                    
                    #print('dbA.cid: %s' % dbA.cid)
                    dbColonyIDList = []
                    cidError = 0 # assume there's no error
                    if dbA.cid != '':
                        dbColonyIDList.append(str.lower(dbA.cid)) # for lower case compare
                    #print('dbColonyIDList: %s' % dbColonyIDList)
                    #print('IncColonyID: %s' % colonyID)

                    # Requirement 7.2.D4a Allele ID match, Colony ID Mismatch
                    if dbColonyIDList != [] and str.lower(colonyID) not in dbColonyIDList:
                        alleleIdMatchColonyIDMismatchList.append('%s%s%s%s%s%s%s%s%s' % (lineNum, TAB, alleleID, TAB, dbA.asym, TAB, dbA.cid, TAB, line))
                        hasError = 1
                        cidError = 1

                    # Requirement 7.2.D4b Colony ID matches MULTIPLE  alleles in the database
                    # get the set of allele(s) (0..n) associated with incoming cid
                    # if there are multiple alleles in the set report
                    if str.lower(colonyID) in colonyToAlleleDict:
                        #print('%s in colonyToAlleleDict' % colonyID)
                        allelesByCidList = colonyToAlleleDict[str.lower(colonyID)]
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
                    if hasError == 0 and cidError == 0:
                        alleleFound = 1
                        if dbColonyIDList == []:
                            # Requirement 7.2.D4 if no error and no cid in the database, add a 
                            # new note to the allele
                            fpNoteload.write('%s%s%s%s' % (alleleID, TAB, colonyID, CRT))

            else: # Requirement 7.2.C1 Allele ID not in MGI OR matches different object type
                #print('Allele ID not in MGI OR matches a different object type')
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
            #print('  #### Allele ID not in input, colonyID is: %s' % colonyID)
            if str.lower(colonyID) in colonyToAlleleDict: 
                alleleList = colonyToAlleleDict[str.lower(colonyID)]
                # Requirement 7.2.F1 Colony ID Matches Multiple Alleles in MGI
                if len(alleleList) > 1:
                    for dbA in alleleList:
                        # report multiple alleles for a colony ID
                        cidMatchToMultiList.append('%s%s%s%s%s%s%s%s' % (lineNum, TAB, dbA.aid, TAB, dbA.asym, TAB, line, CRT))
                        hasError = 1
                        #print('  ###  multiple alleles for colony ID, skip remaining checks')
                        linesSkippedCt += 1
                        continue # multiple alleles for colony ID, don't do further checks

                dbA = alleleList[0] # there is only one
                # Requirement 7.2.F3 allele Status Check
                if dbA.ast != 'Approved':  
                    cidMatchAlleleStatusDiscrepList.append('%s%s%s%s%s%s%s%s%s' % \
                        (lineNum, TAB,  dbA.aid, TAB, dbA.asym, TAB, dbA.ast, TAB, line))
                    hasError = 1
                else:
                    # The following two checks could be replaced with a 
                    # calculated allele symbol match
                    # Requirement 7.2.F2 Marker ID check
                    if markerID != dbA.mid:
                        cidMatchMarkerIdMismatchList.append( '%s%s%s%s%s%s%s' % \
                            (lineNum, TAB, dbA.aid, TAB, dbA.asym, TAB, line))
                        hasError = 1
                    # Requirement 7.2.F2 Allele symbol check
                    #print('cid match, alleleSymbol: %s' % alleleSymbol)
                    #print('cid match dbAlleleSymbol: %s' % dbA.asym)
                    if str.find(dbA.asym, alleleSymbol) == -1:
                        cidMatchAlleleSSMismatchList.append('%s%s%s%s%s%s%s' % \
                            (lineNum, TAB, dbA.aid, TAB, dbA.asym, TAB, line))
                        hasError = 1
                if hasError == 0:
                    alleleFound = 1
                # END COLONY ID MATCH
                
            # BEGIN NO COLONY ID MATCH Requirement 7.2.G
            else:
                #print('  #### No allele ID, no cid match, calculate allele symbol and check in DB')

                symbolError = 0 # default to no error

                # if match found in database, these will be assigned.
                aID = '' 
                symbol = ''

                # Requirement 7.2.H Allele symbol check
                results = findAlleleBySymbol(calcAlleleSymbol)

                # Requirement 7.2.H No CID Match, Allele Symbol Match
                #elif len(results) == 1:
                if len(results) == 1:
                    #print('no cid match, but symbol match: %s' % results)
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
                        # if there is a cid for the symbol it has to be a 
                        # mismatch with the inc cid
                        if allele.cid != '': 
                            symbolMatchColonyIdMismatchList.append('%s%s%s%s%s%s%s%s%s' % \
                                (lineNum, TAB, aID, TAB, symbol, TAB, allele.cid, TAB, line))
                            symbolError = 1
                            hasError = 1
                    # Requirement 7.2.H4 No CID Match, Symbol match, and no errors
                    # add new note to the allele
                    if symbolError == 0 and hasError == 0:
                        alleleFound = 1
                        fpNoteload.write('%s%s%s%s' % (aID, TAB, colonyID, CRT))

                # Requirement 7.2.H3 check for multiple (duplicate) alleles in the database
                #else: # len(results) > 1:
                elif len(results) > 1:
                    #print('no cid match, symbol match to dupe alleles in database: %s' % results)
                    for r in results:
                        #print(r)
                        aID = r['accid']
                        symbol = r['symbol']
                        symbolMatchMultiAlleleList.append('%s%s%s%s%s%s%s' % \
                            (lineNum, TAB, aID, TAB, symbol, TAB, line))
                        symbolError = 1
                        hasError = 1

                # else no results, so create allele if labCode in DB

        # if we've found the allele in MGI count it and go to next line
        if alleleFound:
            allelesFoundCt +=1
            continue

        # if we have not found the allele in MGI, but we have errors, count
        # and continue
        elif hasError:
            linesSkippedCt += 1
            continue  

        # Requirement 7.2.I So, we have a new allele; check the lab code in 
        # the symbol to make sure it is in the Cell Line Lab Code vocab
        labName = ''
        #print('  #### checking allele nomenclature')
        if len(re.findall(pattern1, alleleSymbol)) > 1 or len(re.findall(pattern2, alleleSymbol)) > 1:
            #print('  #### bad allele nomen, not creating allele')
            badNomenList.append('%s%s%s%s%s' % (lineNum, TAB, alleleSymbol, TAB, line))
            hasError = 1
        #print('  #### checking lab code')
        labCode = findLabCode(alleleSymbol)

        if labCode not in labCodeDict:
            labCodeNotInMgiList.append('%s%s%s%s%s' % (lineNum, TAB, labCode, TAB, line))
            #print('  #### bad lab code, not creating allele')
            hasError = 1

        # translate allele type. The key is the pipe-delim IMPC alleleType
        # and subType, value is pipe-delim MGI alleleType and subType
        # impc key and mgi value may not have a subtype - thefore no pipe
        # mgi alleleType and subType may be multi-valued ';' delimited

        # Requirement 7.2.A1g  check that the allele type/subtype 'key' has a translation
        if alleleSubType == '':
            impcKey = str.lower(alleleType)
        else:
            impcKey = '%s|%s' % (str.lower(alleleType), str.lower(alleleSubType))
        if impcKey not in alleleTypeTransDict:
            atTransKeyNotInMgiList.append('%s%s%s%s%s' % (lineNum, TAB, impcKey, TAB, line))
            #print('  #### alleleType/subType combo not in translation')
            hasError = 1
        else:
            mgiValue = alleleTypeTransDict[impcKey]
            #print('impcKey: %s mgiValue: %s' % (impcKey, mgiValue))
            # The case where there is no subtype
            mgiAlleleType = mgiValue
            mgiSubType = ''

            # The case where there is a subtype
            if str.find(mgiValue, '|') != -1:
                mgiAlleleType, mgiSubType = str.split(mgiValue, '|')
            #print('mgiAlleletype: %s mgiSubType: %s' % ( mgiAlleleType, mgiSubType))

            # both mgiAlleleType and mgiSubType can be multivalued
            alleleTypes = str.split(mgiAlleleType, ';')
            subTypes = str.split(mgiSubType, ';')
            #print('alleleTypes: %s subTypes: %s' % (alleleTypes, subTypes))

                
        #
        # If no errors write out to dictionary of lines
        #
        if hasError == 1: # error in the lab code
            linesSkippedCt += 1
        else:
            #print('  #### No allele identified in DB and no errors; translate stuff and create allele')
            # get the sequencNum from the allele
            seqNumFinder = re.compile ( '<em(.*)\(' )
            match = seqNumFinder.search(alleleSymbol)
            sequenceNum = match.group(1)

            # get the lab name from the lab code
            labName = labCodeDict[labCode]

            # calculate allele name
            alleleName = alleleNameTemplate % (sequenceNum, labName)
            if calcAlleleSymbol not in calcAlleleDict:
                calcAlleleDict[calcAlleleSymbol] = []
            
            # the line we want to write to the allele file if no dupes
            alleleLine = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (markerID, TAB, markerSymbol, TAB, mgiAlleleType, TAB, alleleDescription, TAB, colonyID, TAB, strain, TAB, calcAlleleSymbol, TAB, alleleName, TAB, inHeritMode, TAB, alleleClass, TAB, mgiSubType, TAB, alleleStatus, TAB, transmissionState, TAB, alleleCollection, TAB, jNumber, TAB, createdBy, CRT)
            calcAlleleDict[calcAlleleSymbol].append([alleleLine, lineNum, line])

    for key in calcAlleleDict:
        if len(calcAlleleDict[key]) > 1: # dupe in input
            #print('  ### Dupe alleles in input')
            #print(calcAlleleDict[key])
            for l in calcAlleleDict[key]:
                dupeAlleleInInputList.append('%s%s%s%s' % (l[1], TAB, l[2], CRT))
        else:
            linesLoadedCt += 1
            fpAllele.write(calcAlleleDict[key][0][0])

    return 0

def writeQCReport():
    # Purpose: write all QC errors to the QC report file
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: writes to the file system
    # Throws: Nothing

    fpQC.write('Total lines in the input file (%s:%s) including header: %s%s%s' % (host, impcFile, lineNum, CRT, CRT))
    fpQC.write('Total alleles found in the DB: %s%s%s' % (allelesFoundCt, CRT, CRT))
    fpQC.write('Total lines from the input file loaded: %s%s%s' % (linesLoadedCt, CRT, CRT))
    fpQC.write('Total lines from the input file skipped: %s%s%s' % (linesSkippedCt, CRT, CRT))

    fpQC.write('%s%s7.2.A.1 Required Value Missing or Invalid%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%s Missing Value(s)%sInput Line%s' % (TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(missingRequiredValueList):
        fpQC.write(str.join(CRT, missingRequiredValueList))
    fpQC.write('Total: %s' % len(missingRequiredValueList))

    fpQC.write('%s%s7.2.A1 MGI Marker ID not in MGI%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Line%s' % (TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(markerIdNotInMgiList):
         fpQC.write(str.join(CRT, markerIdNotInMgiList))
    fpQC.write('Total: %s' % len(markerIdNotInMgiList))

    fpQC.write('%s%s7.2.A1 Colony Background Strain not in MGI%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Line%s' % (TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(strainNotInMgiList):
         fpQC.write(str.join(CRT, strainNotInMgiList))
    fpQC.write('Total: %s' % len(strainNotInMgiList))

    fpQC.write('%s%s7.2.A1 Allele Class not Endonuclease-mediated%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Line%s' % (TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(unknownAlleleClassList):
         fpQC.write(str.join(CRT, unknownAlleleClassList))
    fpQC.write('Total: %s' % len(unknownAlleleClassList))
   
    fpQC.write('%s%s7.2.A1 Allele (mutation) Type not in Translated Set%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Line%s' % (TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(unknownAlleleTypeList):
         fpQC.write(str.join(CRT, unknownAlleleTypeList))
    fpQC.write('Total: %s' % len(unknownAlleleTypeList)) 

    fpQC.write('%s%s7.2.A1 Allele Subtype not in Translated Set%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Line%s' % (TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(unknownSubTypeList):
         fpQC.write(str.join(CRT, unknownSubTypeList))
    fpQC.write('Total: %s' % len(unknownSubTypeList))

    fpQC.write('%s%s7.2.C1 MGI Allele ID present, No MGI Allele Match%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sObjectType%sInput Line%s' % (TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdNotInMGIList):
         fpQC.write(str.join(CRT, alleleIdNotInMGIList))
    fpQC.write('Total: %s' % len(alleleIdNotInMGIList))

    fpQC.write('%s%s7.2.D3 Allele ID Match, Allele Status Discrepancy%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sAllele Status%sInput Line%s' % (TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdMatchAlleleStatusDiscrepList):
         fpQC.write(str.join(CRT, alleleIdMatchAlleleStatusDiscrepList))
    fpQC.write('Total: %s' % len(alleleIdMatchAlleleStatusDiscrepList))

    fpQC.write('%s%s7.2.D1 Allele ID Match, Marker ID Mismatch%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sAllele ID%sAllele Symbol%sDB Marker ID%sDB Marker Symbol%sInput Line%s' % (TAB, TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdMatchMarkerIdMismatchList):
         fpQC.write(str.join(CRT, alleleIdMatchMarkerIdMismatchList))
    fpQC.write('Total: %s' % len(alleleIdMatchMarkerIdMismatchList))

    fpQC.write('%s%s7.2.D2 Allele ID Match, Allele Symbol  Mismatch%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sInput Line%s' % (TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdMatchAlleleSSMismatchList):
         fpQC.write(str.join(CRT, alleleIdMatchAlleleSSMismatchList))
    fpQC.write('Total: %s' % len(alleleIdMatchAlleleSSMismatchList))

    fpQC.write('%s%s7.2.D4a Allele ID Match, Colony ID Mismatch%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sAllele ID%sDB Allele Symbol%sDB CID%sInput Line%s' % (TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdMatchColonyIDMismatchList):
         fpQC.write(str.join(CRT, alleleIdMatchColonyIDMismatchList))
    fpQC.write('Total: %s' % len(alleleIdMatchColonyIDMismatchList))

    fpQC.write('%s%s7.2.D4b Allele ID match, Colony ID Match to Multi MGI Alleles%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Allele ID%sInput Allele Symbol%sInput Colony ID%sDB Allele ID%sDB Allele Symbol%sDB Allele Type%sDB Colony ID%sInput Line%s' % (TAB, TAB, TAB, TAB, TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdMatchColonyIdMatchToMultiList):
         fpQC.write(str.join(CRT, alleleIdMatchColonyIdMatchToMultiList))
    fpQC.write('Total: %s' % len(alleleIdMatchColonyIdMatchToMultiList))

    fpQC.write('%s%s7.2.D4b Allele ID match, Colony ID Match to Different Allele%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Allele ID%sInput Allele Symbol%sInput Colony ID%sDB Allele ID%sDB Allele Symbol%sDB Allele Type%sDB Colony ID%sInput Line%s' % (TAB, TAB, TAB, TAB, TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(alleleIdMatchColonyIdMatchToDiffAlleleList):
         fpQC.write(str.join(CRT, alleleIdMatchColonyIdMatchToDiffAlleleList))
    fpQC.write('Total: %s' % len(alleleIdMatchColonyIdMatchToDiffAlleleList))

    fpQC.write('%s%s7.2.F1 Colony ID Matches Multiple Alleles%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sInput Line%s' % (TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(cidMatchToMultiList):
         fpQC.write(str.join(CRT, cidMatchToMultiList))
    fpQC.write('Total: %s' % len(cidMatchToMultiList))

    fpQC.write('%s%s7.2.F2a Colony ID Match, Marker ID Mismatch%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sInput Line%s' % (TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(cidMatchMarkerIdMismatchList):
         fpQC.write(str.join(CRT, cidMatchMarkerIdMismatchList))
    fpQC.write('Total: %s' % len(cidMatchMarkerIdMismatchList))

    fpQC.write('%s%s7.2.F2b Colony ID Match, Allele Symbol Mismatch%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sInput Line%s' % (TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(cidMatchAlleleSSMismatchList):
         fpQC.write(str.join(CRT, cidMatchAlleleSSMismatchList))
    fpQC.write('Total: %s' % len(cidMatchAlleleSSMismatchList))

    fpQC.write('%s%s7.2.F3 Colony ID Match, Allele Status Discrepancy%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sDB Allele Status%sInput Line%s' % (TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(cidMatchAlleleStatusDiscrepList):
         fpQC.write(str.join(CRT, cidMatchAlleleStatusDiscrepList))
    fpQC.write('Total: %s' % len(cidMatchAlleleStatusDiscrepList))

    fpQC.write('%s%s7.2.H1 Allele Symbol Match, Allele Status Discrepancy%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sDB Allele Status%sInput Line%s' % (TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(symbolMatchAlleleStatusDiscrepList):
         fpQC.write(str.join(CRT, symbolMatchAlleleStatusDiscrepList))
    fpQC.write('Total: %s' % len(symbolMatchAlleleStatusDiscrepList))

    fpQC.write('%s%s7.2.H2 Allele Symbol Match, Colony ID Mismatch%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sDB Allele CID%sInput Line%s' % (TAB, TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(symbolMatchColonyIdMismatchList):
         fpQC.write(str.join(CRT, symbolMatchColonyIdMismatchList))
    fpQC.write('Total: %s' % len(symbolMatchColonyIdMismatchList))

    fpQC.write('%s%s7.2.H3 Allele Symbol Match to Multiple Alleles%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sDB Allele ID%sDB Allele Symbol%sInput Line%s' % (TAB, TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(symbolMatchMultiAlleleList):
         fpQC.write(str.join(CRT, symbolMatchMultiAlleleList))
    fpQC.write('Total: %s' % len(symbolMatchMultiAlleleList))

    fpQC.write('%s%sNew check: Allele Symbol has incorrect nomenclature%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sAllele Symbol%sInput Line%s' % (TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(badNomenList):
        fpQC.write(str.join(CRT, badNomenList))
    fpQC.write('Total: %s' % len(badNomenList))

    fpQC.write('%s%s7.2.I No Allele Match, Lab Code not Present%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sLab Code%sInput Line%s' % (TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)

    if len(labCodeNotInMgiList):
        fpQC.write(str.join(CRT, labCodeNotInMgiList))
    fpQC.write('Total: %s' % len(labCodeNotInMgiList))
   
    fpQC.write('%s%s7.2.A1g Allele (mutation) Type/Allele Subtype combination not in Translated Set%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sIMPC alleleType|subType%sInput Line%s' % (TAB, TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(atTransKeyNotInMgiList):
        fpQC.write(str.join(CRT, atTransKeyNotInMgiList))
    fpQC.write('Total: %s' % len(atTransKeyNotInMgiList))
 
    fpQC.write('%s%sDuplicate Allele in Input%s%s' % (CRT, CRT, CRT, CRT))
    fpQC.write('Line#%sInput Line%s' % (TAB, CRT))
    fpQC.write('_____________________________________________________________%s' % CRT)
    if len(dupeAlleleInInputList):
        fpQC.write(str.join(CRT, dupeAlleleInInputList))
    fpQC.write('Total: %s' % len(dupeAlleleInInputList))

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

if closeFiles() != 0:
    sys.exit(1)

db.useOneConnection(0)

sys.exit(0)
