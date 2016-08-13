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
# This set of scripts script will have multiple stages in order to build a 
# complete map file. These requires a number of input files, namely:
# 1) UniProt to UniRef map
# 2A) SwissProt.dat complete datasets 
# 2B) Trembl.dat complete dataset 
# 3) TSV for GO annotations with evidence code, accession, and ref ID
# 4) Map file for relevant GO accessions to UniProt accessions
#
# HOWTO:  
# ./build_map_phase_1.py /path_to_uniprot_uniref_map /path_to_go_uniprot_map /path_to_go_data_tsv
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
# These scripts were broken up into phases as memory requirements (input files are huge) were limiting the ability
# to use cheap hash lookups and forcing iteration over huge lists. With each step happening in phases grabbing
# just bits of the data, the overall time for completion is much less. Scripts are named with respect to the order
# they need to occur in (build_map_phase_1.py --> build_map_phase_2.py --> etc.)
#
# Author: James Matsumura

import sys, os, re, gzip

uniprot_uniref_map =  str(sys.argv[1]) # important to get the same dated versions of all UniProt files 
go_uniprot_map =  str(sys.argv[2]) 
go_tsv =  str(sys.argv[3]) 

prot_ref_map_file = gzip.open(uniprot_uniref_map, 'rb') 
go_tsv_file = open(go_tsv, 'r') 
go_prot_map_file = open(go_uniprot_map, 'r') 
outFile1 = open('./map_file.v1.tsv', 'w')
outFile2 = './phase_1.tsv'

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

regexForMappedAccession = r"UniRef100\_(\w+)"
regexForFBgnIds = r"[A-Z]+[a-z]*(\d+)"
regexForGOid = r":(.*)"

uniquePMIds = set()
uniqueUnirefIds = set()
entry1List = []
entry2List = []
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
# Want to start with this since it'd be a waste of time to find the evidence
# for those GO entries that don't have a corresponding UniRef entity.
for line in go_prot_map_file:
	mappings = line.split('\t')
	# Ideally there would be no duplicate GO noted IDs and no duplicate
	# UniProts. Since this is not the case, simply append every UniProt
	# acc linked to a particular GO noted ID. 
	goData.setdefault(mappings[0], []).append(mappings[1])

print 'stage3'
# Just gathering the relevant data, not building the final file yet. 
for line in go_tsv_file:
	elements = line.split('\t')

	if elements[1] == '':
		continue
	else:
		ref_id = ':'.join([elements[2], elements[3]])
		entry2List.append(Entry2(go_noted_acc=elements[1],evidence_type=elements[0],reference_id=ref_id,go_term=elements[4].strip(' ')))

print 'stage4'
# Now all the data has been gathered, build the final map file.
for x in entry2List:
	outFile1.write('\t'.join([x.go_noted_acc, x.evidence_type, x.reference_id, x.go_term]))

outFile1.close() # switch to reading it

# This next step isn't ideal, but only need to run this step once.

print 'stage5'
# First, add in the UniProt accs related to the noted GO ID
with open('./map_file.v1.tsv', 'r') as input_file, open(outFile2, 'w') as output_file:
	relevant = False
	go_to_uniprot = ''
	for line in input_file:
		line = line.replace('\n','')
		elements = line.split('\t')
		if ':' in elements[0]:
			if 'FB:' in elements[0]:
				elements[0] = 'FBGN' + re.match(regexForFBgnIds, elements[0]).group(1)	
			else:
				elements[0] = re.match(regexForGOid, elements[0]).group(1)	
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
