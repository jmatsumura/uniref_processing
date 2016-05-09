#!/usr/bin/python

# The purpose of this script is to build a custom UniRef100 fasta dataset
# to run RapSsearch2 against. The target subset is dependent on the input
# files fed to this program. The data needs to be filtered so that it 
# consists of a list of relevant UniProt IDs in the second tab-delimited
# column of the file. 
#
# These IDs will then be used to subset the UniRef100 file.  
# relevant UniProt IDs is is just bacterial proteins. In 
# order to accomplish this, the most recent sprot file for just bacterial 
# proteins needs to be obtained to isolate only those UniRef clusters that
# have a bacterial representative. The script has three stages:
# 
# 1) Build a set representing all the UniRef100 related accessions present 
# that have some experimental evidence denoted by their GO annotation
# 2) map these IDs to their current UniRef representative
# 3) use these accessions to build a subset for those uniref hits with
# this level of annotation
#
# HOWTO:
# ./build_goset_uniref100.py path_to_go_accs_file path_to_uniref_file path_to_map_file
#
# Author: James Matsumura

import sys, os, re, gzip

goFile =  str(sys.argv[1]) 
unirefFile =  str(sys.argv[2]) # important that this is same version # as map file
mapFile =  str(sys.argv[3]) 

theGoFile = open(goFile, 'r') 
origUnirefFile = gzip.open(unirefFile, 'rb') 
mappingFile = gzip.open(mapFile, 'rb') 
relevantUnirefFile = open('./go_to_uniref_with_evidence.txt', 'w')
outFile = gzip.open('./custom_goev_uniref100.fasta.gz', 'wb')

footerFound = False
accessionFound = False
relevantUnirefEntry = False 

regexForUnirefAccession = r"^>UniRef100\_(\w+)\s+.*"
regexForMappedAccession = r"UniRef100\_(\w+)"
uniqueIds = set()
uniqueUnirefIds = set()

print "stage 1"
# 1)
# This file has already been preprocessed using bash/vim so that the second column
# contains the relevant UniProt IDs linked to some GO annotation that had 
# experimental evidence to back it.
for line in theGoFile:

	line = line.strip('\n')
	extractUniprot = line.split('\t')
	print extractUniprot[1]
	uniqueIds = uniqueIds | {extractUniprot[1]}

print "stage 2"
# 2) 
# Must map each UniProt entry to its corresponding current UniRef representative.
for line in mappingFile:
	elements = line.split('\t')
	if elements[0] in uniqueIds and elements[7] not in uniqueUnirefIds:
		mappedId = re.search(regexForMappedAccession, elements[7])
		finalId = mappedId.group(1)
		uniqueUnirefIds = uniqueUnirefIds | {finalId}
		relevantUnirefFile.write(finalId+'\n')
	
print "stage 3"	
# 3) 
# Each UniRef entry is denoted with the UniRef identity level followed by
# the UniProt accession cluster representative like so:
# UniRef100_Q6GZX4. This will have been generated from Step 2.
for line in origUnirefFile:

	if(line.startswith('>')):
		relevantUnirefEntry = False # must be reset each entry
		# Some odd formatting in the UniRef file? need to 
		# actually check to make sure it's in proper format
		findEntry = re.search(regexForUnirefAccession, line)
		if(findEntry):
			foundEntry = findEntry.group(1)
			# Perhaps the accession wasn't included in the map file. This ideally
			# should have no impact if the map file was made correctly but I'm
			# adding it just in case. 
			if(foundEntry in uniqueUnirefIds or foundEntry in uniqueUnirefIds):
				relevantUnirefEntry = True
				outFile.write(line)	

	elif(relevantUnirefEntry == True):
		outFile.write(line)	
