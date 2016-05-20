#!/usr/bin/python
#
# This script follows phase 3 and simply requires the same input file.
#
# Author: James Matsumura

import sys, os, re, gzip

uniprot_uniref_map =  str(sys.argv[1]) # important to get the same dated versions of all UniProt files 

prot_ref_map_file = gzip.open(uniprot_uniref_map, 'rb') 
outFile = './phase_3.5.tsv'

regexForMappedAccession = r"UniRef100\_(\w+)"
protData = {}

print 'stage1'
# Begin building the list of Entry objects
for line in prot_ref_map_file:
	mappings = line.split('\t')
	protData[mappings[0]] = mappings[6]

print 'stage2'
# Now that UniProt accs are present, map to UniRef accs. It is important to 
# remember the size of the objects for this loop. Only want to iterate over
# the largest object (entry1List) if we absolutely have to and always be sure
# to break at the point of finding the relevant info. 
with open('./phase_3.tsv', 'r') as input_file, open(outFile, 'w') as output_file:
	for line in input_file:
		relevant = False
		prot_to_go = ''
		line = line.replace('\n','')
		elements = line.split('\t')
		if elements[4]=='':
			# Can uncomment out if ALL go data is needed, but since these don't
			# map to a UniProt acc I don't see the point of including currently
			# as we can't map these from a UniRef100 match.
			#output_file.write(line + '\t' + '\t' + '\n') # keep the tabs consistent
			continue
		elif prot_to_go == '': # need to perform this block first
			# Do a check for those with multiple accs and those with only one.
			# In order to have proper mapping, need to iterate over these
			# individually
			if ',' in elements[4]:
				individualAccs = elements[4].split(',')
				for j in individualAccs:
					if j in protData:
						go_data = protData[j]
					else:
						relevant = False
					if go_data != None or go_data != '':
						relevant = True
						if prot_to_go != '': # append commas
							prot_to_go += ','
						prot_to_go += go_data
					if not relevant:
						if prot_to_go != '': # append commas
							prot_to_go += ','
						prot_to_go += 'NONE'
			else:
				if elements[4] in protData:
					go_data = protData[elements[4]]
				if go_data != None or go_data != '':
					relevant = True
					prot_to_go += go_data
		if relevant:
			output_file.write(line + '\t' + prot_to_go + '\n')
