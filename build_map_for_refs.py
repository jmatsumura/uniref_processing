#!/usr/bin/python

# The purpose of this script is to build a map for the relevant references
# for each hit against the custom UniRef100 dataset. This means capturing
# the following:
# 1) UniProt accession
# 2) UniRef representative
# 3) SwissProt evidence
# 4) SwissProt reference
# 5) GO evidence
# 6) GO reference
#
# Note that since we only care about experimental evidence from UniProt, 
# UniProt accessions with evidence will only be derived from SwissProt 
# and denoted by the evidence code ECO:0000269.
# 
# This script will have multiple stages in order to build a complete map file.
# It requires a number of input files, namely:
# 1) UniProt to UniRef map
# 2A) SwissProt.dat complete datasets 
# 2B) Trembl.dat complete dataset 
# 3) TSV for GO annotations with evidence code, accession, and ref ID
# 4) Map file for relevant GO accessions to UniProt accessions
#
# HOWTO:  testing 
# ./build_map_for_refs.py /path_to_uniprot_uniref_map /swissprot.dat /trembl.dat /go_refs.tsv /go_to_uniprot_map
#
# Author: James Matsumura

import sys, os, re, gzip

uniprot_uniref_map =  str(sys.argv[1]) # important to get the same dated versions of all UniProt files 
#sprot_dat =  str(sys.argv[2]) 
#go_tsv =  str(sys.argv[4]) 
#go_uniprot_map =  str(sys.argv[5]) 

prot_ref_map_file = gzip.open(uniprot_uniref_map, 'rb') 
#sprot_file = gzip.open(sprot_dat, 'rb') 
#go_tsv_file = open(go_tsv, 'r') 
#go_prot_map_file = open(go_uniprot_map, 'r') 
#outFile = gzip.open('./m.fasta.gz', 'wb')

# This object will house the seven attributes noted in the comments above. 
class Entry1:
	def __init__(self, prot_acc, ref_acc, go_terms):
		self.prot_acc = prot_acc
		self.ref_acc = ref_acc
		#self.sprot_ev_code = None
		#self.sprot_ref = None
		#self.go_ev_code = None
		#self.go_ref = None
		self.go_terms = go_terms

	#def setSprotValues(self, sprot_ev_code_value, sprot_ref_value):
	#	self.sprot_ev_code = sprot_ev_code_value
	#	self.sprot_ref = sprot_ref_value

	#def setGOvalues(self, go_ev_code_value, go_ref_value):
	#	self.go_ev_code = go_ev_code
	#	self.go_ref = go_ref_value

footerFound = False
accessionFound = False
relevantUnirefEntry = False 

regexForAccession = r"^AC\s+(.*);"
regexForFooter = r"^\/\/$"
regexForMappedAccession = r"UniRef100\_(\w+)"
regexForSprotReferences = r"ECO:0000269\|(PubMed:\d+)"

uniqueIds = set()
uniquePMIds = set()
uniqueUnirefIds = set()
entryList = []
sprotData = {}

print 'stage1'
# Begin building the list of Entry objects
for line in prot_ref_map_file:
	mappings = line.split('\t')
	# appears that not all entries have an UniRef100 ID
	if('UniRef100' in mappings[7]):
		uniref_acc = re.search(regexForMappedAccession, mappings[7]).group(1)
	else:
		uniref_acc = None
	entryList.append(Entry1(prot_acc=mappings[0],ref_acc=uniref_acc, go_terms=mappings[6]))
print 'stage2'
# Just gather the data from the sprot file, add these values to their objects later
for line in sprot_file:
	if (footerFound == True): # reinitialize values for next record
		accessionFound = False
		footerFound = False
		uniquePMIds.clear()
	elif(accessionFound == True):
		if(re.search(regexForFooter, line)):
			if(';' in foundAccession):
				multiAccessions = foundAccession.split('; ')
				for x in multiAccessions: # iterate over this ~2-3 len list
					sprotData[x] = '|'.join(uniquePMIds)
					print('|'.join(uniquePMIds))
			footerFound = True
		else:
			if('ECO:0000269' in line):
				pmid = re.search(regexForSprotReferences, line).group(1)
				if pmid not in uniquePMIds:
					uniquePMIds = uniquePMIds | {pmid}
	else:
		findAccession = re.search(regexForAccession, line)
		if(findAccession):
			foundAccession = findAccession.group(1)
			accessionFound = True
# Same as the previous loop, just gather the GO data here
#for line in go_tsv_file:

# Again, just linking data here not completing the objects yet
#for line in go_tsv_file:

# Fill out the rest of the object attributes
#for x in entryList:

print(entry_list[1].ref_acc)
