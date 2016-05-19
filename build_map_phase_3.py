#!/usr/bin/python
#
# This script follows phase 2 and actually needs to be run twice. The first time it needs to be given
# the number '7' as its second option and the second time it requires the number '6'. This is simply
# pulling different data fields from the UniProt mapping file. Again, as noted in the phase 1 script,
# this is only being done sequentially to overcome the memory limitations due to the large input files. 
#
# Author: James Matsumura

import sys, os, re, gzip

uniprot_uniref_map =  str(sys.argv[1]) # important to get the same dated versions of all UniProt files 
col =  int(sys.argv[2]) # important to get the same dated versions of all UniProt files 

prot_ref_map_file = gzip.open(uniprot_uniref_map, 'rb') 
outFile = './phase_3.tsv'

regexForMappedAccession = r"UniRef100\_(\w+)"
protData = {}

print 'stage1'
# Begin building the list of Entry objects
for line in prot_ref_map_file:
	mappings = line.split('\t')
	# appears that not all entries have an UniRef100 ID
	if 'UniRef100' in mappings[col]:
		uniref_acc = re.search(regexForMappedAccession, mappings[col]).group(1)
	else:
		uniref_acc = None
	protData[mappings[0]] = uniref_acc

print 'stage2'
# Now that UniProt accs are present, map to UniRef accs. It is important to 
# remember the size of the objects for this loop. Only want to iterate over
# the largest object (entry1List) if we absolutely have to and always be sure
# to break at the point of finding the relevant info. 
with open('./phase_2.tsv', 'r') as input_file, open(outFile, 'w') as output_file:
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
					if j in protData:
						uniref = protData[j]
					else:
						relevant = False
					if uniref != '':
						relevant = True
						if ref_to_prot != '': # append commas
							ref_to_prot += ','
						ref_to_prot += uniref
					if not relevant:
						if ref_to_prot != '': # append commas
							ref_to_prot += ','
						ref_to_prot += 'NONE'
			else:
				if elements[4] in protData:
					uniref = protData[elements[4]]
				if uniref != '':
					relevant = True
					ref_to_prot += uniref
		if relevant:
			output_file.write(line + '\t' + ref_to_prot + '\n')
