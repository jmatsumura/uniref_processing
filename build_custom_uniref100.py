#!/usr/bin/python

# The purpose of this script is to build a custom UniRef100 fasta dataset
# to run RapSsearch2 against. The target is just bacterial proteins. In 
# order to accomplish this, the most recent sprot file for just bacterial 
# proteins needs to be obtained to isolate only those UniRef clusters that
# have a bacterial representative. The script has three stages:
# 
# 1) isolate protein accessions that are annotated with some level of 
# experimental evidence
# 2) map these IDs to their current UniRef representative
# 3) use these accessions to refine any uniref hits that match a protein with
# this level of annotation
#
# HOWTO:
# ./build_custom_uniref100.py path_to_sprot_file path_to_uniref_file path_to_map_file
#
# Author: James Matsumura

import sys, os, re, gzip, time
from Bio import SwissProt

sprotFile =  str(sys.argv[1]) 
unirefFile =  str(sys.argv[2]) # important that this is same version as map file
mapFile =  str(sys.argv[3]) 

bacteriaOnlyFile = gzip.open(sprotFile, 'rb') # large files, use compression
origUnirefFile = gzip.open(unirefFile, 'rb') 
mappingFile = gzip.open(mapFile, 'rb') 
relevantBacteriaFile = open('./bacteria_with_evidence.txt', 'w')
relevantUnirefFile = open('./uniref_with_evidence.txt', 'w')
outFile = gzip.open('./custom_uniref100.fasta.gz', 'wb')

# Only want to find those with experimental evidence backing the annotation.
# Use: http://bioportal.bioontology.org/ontologies/ECO/?p=classes&conceptid=root
evidenceCodes = ('ECO:0000269','ECO:0000006','ECO:0000179', \
     	         'ECO:0000360','ECO:0005606','ECO:0000325', \
		 	 	 'ECO:0000180','ECO:0005604','ECO:0000002', \
		 		 'ECO:0005605','ECO:0000073','ECO:0000059', \
		 		 'ECO:0000008','ECO:0001094','ECO:0005516', \
		 		 'ECO:0000021','ECO:0000340','ECO:0000220', \
		 		 'ECO:0005504','ECO:0005031')

footerFound = False
accessionFound = False
relevantUnirefEntry = False 

regexForAccession = r"^AC\s+(.*);"
regexForFooter = r"^\/\/"
regexForUnirefAccession = r"^>UniRef100\_(\w+)\s+.*"
regexForMappedAccession = r"UniRef100\_(\w+)"
uniqueIds = set()
uniqueUnirefIds = set()

# 1)
# The annotations here are messy. Evidence codes are often in the comments (CC)
# or other tags like feature table (FT) or reference comments (RC). Thus, need
# to check multiple attributes of each entry for any trace of evidence.
print "stage 1"
for line in bacteriaOnlyFile:

	if (footerFound == True): # reinitialize values for next record
		accessionFound = False
		footerFound = False

	elif(accessionFound == True):
		if(re.search(regexForFooter, line)):
			footerFound = True

		else:
			if any(x in line for x in evidenceCodes):
				# It appears that each accession tag can have
				# multiple accessions tied to it. These all go
				# to the same representative in the UniProt site,
				# but, going to include them all as if they were
				# separate entities in case of some timing discrepancies
				# for when the UniRef100 dataset was constructed. 
				if(';' in foundAccession):
					multiAccessions = foundAccession.split('; ')
					for x in multiAccessions: # iterate over this ~2-3 len list
						if x not in uniqueIds:
							uniqueIds = uniqueIds | {x}
							relevantBacteriaFile.write(x+'\n')

				elif(foundAccession not in uniqueIds): 
					uniqueIds = uniqueIds | {foundAccession}
					relevantBacteriaFile.write(foundAccession+'\n')

	else:
		findAccession = re.search(regexForAccession, line)
		if(findAccession):
			foundAccession = findAccession.group(1)
			accessionFound = True

time.sleep(20) # using this to do a check
# 2) 
# Must map each UniProt entry to its corresponding current UniRef representative.
print "stage 2"
for line in mappingFile:
	elements = line.split('\t')
	if elements[0] in uniqueIds and elements[7] not in uniqueUnirefIds:
		mappedId = re.search(regexForMappedAccession, elements[7])
		finalId = mappedId.group(1)
		uniqueUnirefIds = uniqueUnirefIds | {finalId}
		relevantUnirefFile.write(finalId+'\n')
	
time.sleep(20)
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
