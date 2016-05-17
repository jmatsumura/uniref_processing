#!/usr/bin/python
#
# The purpose of this script is to build a map for the relevant references
# for each hit against the custom UniRef100 dataset. This means capturing
# the following:
# 1) UniProt accession
# 2) UniRef representative
# 3) UniRef GO evidence
# 4) SwissProt evidence
# 5) SwissProt reference
# 6) GO evidence (from GO database)
# 7) GO reference
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
# EXAMPLE TAB-DELIMITED OUTPUT FILE:
# -----------------------------------------------------------------------------------------------------------------------------------------------
# | DB:ACC     | GO EV CODE | PM ID  (GO)   | GO term (GO) | UniProt acc | UniRef acc | Go term (UniProt) | PM ID (UniProt) | PM ID (UniRef100) |
# -----------------------------------------------------------------------------------------------------------------------------------------------
# | SGD:S000911| IGI        | PMID:21315072 | GO:0005634   | Q7Z8M2      | C8VAU1     | GO:0043940        |  PMID:20870878  | PMID: 2091923912  |
# -----------------------------------------------------------------------------------------------------------------------------------------------
#
# This file is a map for tying any UniRef match to whatever evidence is tied to any members of the cluster
# whether this evidence is derived from GO (columns 1-4) or UniProt (7-9) databases. Description of each column:
# 1. DB:ACC - This is an accession particular to a database. This could be NCBI, UniProt, FlyBase, etc. GO stores 
# all of its evidence linked to these non-uniform IDs so they were mapped to their respective UniProt acc when possible.
# 2. GO EV CODE - One  of the six experimental evidence codes: EXP / IDA / IPI / IMP / IGI / IEP
# 3. PM ID (GO) - PubMed ID for this entry noted by the GO database
# 4. GO term (GO) - GO term noted by the GO database
# 5. UniProt acc - This was obtained by mapping column 1 to UniProt IDs or, if column 1 is absent, it was pulled
# directly from SwissProt
# 6. UniRef acc - Maps column 5 to its UniRef cluster representative (***this is dynamic with each release***) 
# 7. GO term (UniProt) - GO terms found via the UniProt-UniRef map file from their site
# 8. PM ID (UniProt) - Evidence tied to the manually curated entries indicated by ECO:0000269 in SwissProt
# 9. PM ID (UniRef) - Same technique as 8. This is also included to give the user a choice for whether they want
# to link references that are derived from the cluster members or the cluster representative of the UniRef100 hit. 
#
# Author: James Matsumura

import sys, os, re, gzip
from collections import defaultdict

uniprot_uniref_map =  str(sys.argv[1]) # important to get the same dated versions of all UniProt files 
sprot_dat =  str(sys.argv[2]) 
go_uniprot_map =  str(sys.argv[3]) 
go_tsv =  str(sys.argv[4]) 

prot_ref_map_file = gzip.open(uniprot_uniref_map, 'rb') 
sprot_file = gzip.open(sprot_dat, 'rb') 
go_tsv_file = open(go_tsv, 'r') 
go_prot_map_file = open(go_uniprot_map, 'r') 
outFile1 = open('./map_file.v1.tsv', 'w')
outFile2 = './map_file.v2.tsv'
outFile3 = './map_file.v3.tsv'
outFileFinal = './map_file.final.tsv'

# This object will house the first three attributes noted in the comments above. 
class Entry1:
	def __init__(self, prot_acc, ref_acc, go_terms):
		self.prot_acc = prot_acc
		self.ref_acc = ref_acc
		self.go_terms = go_terms

# This object will house GO data. This means the GO noted accession which can come from a variety of 
# sources like ZFIN, UniProt, RefSeq, etc. as well as the evidence type and reference/source ID. 
class Entry2:
	def __init__(self, go_noted_acc, evidence_type, reference_id, go_term):
		self.go_noted_acc = go_noted_acc
		self.evidence_type = evidence_type
		self.reference_id = reference_id
		self.go_term = go_term

footerFound = False
accessionFound = False
relevantUnirefEntry = False 

regexForAccession = r"^AC\s+(.*);"
regexForFooter = r"^\/\/$"
regexForMappedAccession = r"UniRef100\_(\w+)"
regexForSprotReferences = r"ECO:0000269\|PubMed:(\d+)"
regexForFBgnIds = r"[A-Z]+[a-z]*(\d+)"
regexForGOid = r":(.*)"

uniqueSprotIds = set()
uniquePMIds = set()
uniqueUnirefIds = set()
entry1List = []
entry2List = []
sprotData = {}
goData = {}

print 'stage1'
# Begin building the list of Entry objects
for line in prot_ref_map_file:
	mappings = line.split('\t')
	# appears that not all entries have an UniRef100 ID
	if 'UniRef100' in mappings[7]:
		uniref_acc = re.search(regexForMappedAccession, mappings[7]).group(1)
	else:
		uniref_acc = None
	entry1List.append(Entry1(prot_acc=mappings[0],ref_acc=uniref_acc, go_terms=mappings[6]))

print 'stage2'
# Just gather the data from the sprot file, add these values to their objects later. Note
# that only a hash/dict is needed here as there are only two data points to store. 
for line in sprot_file:
	if footerFound == True: # reinitialize values for next record
		accessionFound = False
		footerFound = False
		uniquePMIds.clear()
	elif accessionFound == True:
		if re.search(regexForFooter, line):
			if ';' in foundAccession:
				multiAccessions = foundAccession.split('; ')
				for x in multiAccessions: # iterate over this ~2-3 len list
					sprotData[x] = '|'.join(uniquePMIds)
			else:
				sprotData[foundAccession] = '|'.join(uniquePMIds)
			footerFound = True
		else:
			if 'ECO:0000269|PubMed' in line:
				pmid = re.search(regexForSprotReferences, line).group(1)
				if pmid not in uniquePMIds:
					uniquePMIds = uniquePMIds | {pmid}
	else:
		findAccession = re.search(regexForAccession, line)
		if findAccession:
			foundAccession = findAccession.group(1)
			accessionFound = True

print 'stage3'
# Want to start with this since it'd be a waste of time to find the evidence
# for those GO entries that don't have a corresponding UniRef entity.
for line in go_prot_map_file:
	mappings = line.split('\t')
	# Ideally there would be no duplicate GO noted IDs and no duplicate
	# UniProts. Since this is not the case, simply append every UniProt
	# acc linked to a particular GO noted ID. 
	goData.setdefault(mappings[0], []).append(mappings[1])

print 'stage4'
# Just gathering the relevant data, not building the final file yet. 
for line in go_tsv_file:
	elements = line.split('\t')

	if elements[1] == '':
		continue
	else:
		ref_id = ':'.join([elements[2], elements[3]])
		entry2List.append(Entry2(go_noted_acc=elements[1],evidence_type=elements[0],reference_id=ref_id,go_term=elements[4].strip(' ')))

print 'stage5'
# Now all the data has been gathered, build the final map file.
for x in entry2List:
	outFile1.write('\t'.join([x.go_noted_acc, x.evidence_type, x.reference_id, x.go_term]))

outFile1.close() # switch to reading it

# These next chunks are inefficient, but I only need to run this code once to get the output and
# the inconsistent formatting of the data pulled made it too much a hassle to reformat only to 
# make this step more efficient. Broken up in such a way via intermediary files used to avoid 
# hitting >n^2 complexity. 

print 'stage6'
# First, add in the UniProt accs related to the noted GO ID
with open('./map_file.v1.tsv', 'r') as input_file, open(outFile2, 'w') as output_file:
	relevant = False
	go_to_uniprot = ''
	for line in input_file:
		line = line.replace('\n','')
		elements = line.split('\t')
		if ':' in elements[0]:
			if 'FB:' in elements[0]:
				elements[0] = 'FBGN' + re.search(regexForFBgnIds, elements[0]).group(1)	
			else:
				elements[0] = re.search(regexForGOid, elements[0]).group(1)	
		for k,v in goData.iteritems():
			if k == elements[0]:
				relevant = True
				go_to_uniprot = ','.join(v)
				break
			else:
				relevant = False
		if relevant:
			output_file.write(line + '\t' + go_to_uniprot.replace('\n','') + '\n')
		else:
			output_file.write(line + '\t' + '\n') # keep the tabs consistent

print 'stage7'
# Now that UniProt accs are present, map to UniRef accs. It is important to 
# remember the size of the objects for this loop. Only want to iterate over
# the largest object (entry1List) if we absolutely have to and always be sure
# to break at the point of finding the relevant info. 
with open(outFile2, 'r') as input_file, open(outFile3, 'w') as output_file:
	for line in input_file:
		relevant = False
		ref_to_prot = ''
		prot_go_dat = ''
		line = line.replace('\n','')
		elements = line.split('\t')
		if elements[4]=='':
			# Can uncomment out if ALL go data is needed, but since these don't
			# map to a UniProt acc I don't see the point of including currently
			# as we can't map these from a UniRef100 match.
			#output_file.write(line + '\t' + '\t' + '\n') # keep the tabs consistent
			continue
		elif ref_to_prot == '': # need to perform this block first
			# Do a check for those with multiple accs and those with only one.
			# In order to have proper mapping, need to iterate over these
			# individually
			if ',' in elements[4]:
				individualAccs = elements[4].split(',')
				for j in individualAccs:
					for x in entry1List:
						if x.prot_acc == j:
							relevant = True
							if ref_to_prot != '': # append commas
								ref_to_prot += ','
								prot_go_dat += ','
							ref_to_prot += x.ref_acc
							prot_go_dat += x.go_terms
							break
					if not relevant:
						if ref_to_prot != '': # append commas
							ref_to_prot += ','
							prot_go_dat += ','
						ref_to_prot += 'NONE'
						prot_go_dat += 'NONE'
			else:
				for x in entry1List:
					if x.prot_acc in elements[4]:
						relevant = True
						ref_to_prot += x.ref_acc
						prot_go_dat += x.go_terms
						break
		if relevant:
			output_file.write(line + '\t' + ref_to_prot + '\t' + prot_go_dat + '\n')

print 'stage8'
# Finally, append the SwissProt data and all the references associated with each UniProt acc
with open(outFile3, 'r') as input_file, open(outFileFinal, 'w') as output_file:
	for line in input_file:
		line = line.replace('\n','')
		elements = line.split('\t')
		uniprot_refs = '' # PM IDs linked to the accs
		uniref_refs = ''
		if elements[4]=='' and elements[5]=='': # no UniRef/UniProt
			continue
		# Should assume that if there's a UniRef, there's a UniProt
		elif not elements[4]=='' and not elements[5]=='': 
			uniqueSprotIds = uniqueSprotIds | {elements[4]} # Track which SwissProt already present
			if ',' in elements[4]:
				uniprots = elements[4].split(',')
				for x in uniprots:
					if uniprot_refs == '':
						uniprot_refs += sprotData[x]
					else:
						uniprot_refs += ';'+sprotData[x]
			else:
				uniprot_refs += sprotData[elements[4]]
			if ',' in elements[5]:
				unirefs = elements[5].split(',') # possible to have multiple here
				for x in unirefs:
					if uniref_refs == '':
						uniref_refs += 'PMID:'+sprotData[x]
					else:
						uniref_refs += ';PMID:'+sprotData[x]
			else:
				uniref_refs = sprotData[elements[5]]
		output_file.write(line + '\t' + uniprot_refs + '\t' + uniref_refs + '\n')	

print 'stage9'
# Up til now, building on the assumption that a GO noted accession is present. However,
# need to be able to map those entries which only were found to have evidence through
# SwissProt. This will then exclude columns 2-4. 
with open(outFileFinal, 'w') as output_file:
	for key in sprotData:
		if key not in uniqueSprotIds:
			unirefId = ''
			goData = ''
			for x in entry1List:
					if x.prot_acc == key:
						unirefId = x.ref_acc
						goData = x.go_terms
						break
			output_file.write('UniProtKB:'+key+'\t'+'\t'+'\t'+key+'\t'+unirefId+'\t'+goData+'\t'+sprotData[x]+'\t'+sprotData[unirefId]+'\n')
