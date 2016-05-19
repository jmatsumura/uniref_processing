#!/usr/bin/python
#
# Intermediate phase to add/append the SwissProt annotations.
#
# Author: James Matsumura

import sys, os, re, gzip

uniprot_uniref_map =  str(sys.argv[1])
sprot_dat =  str(sys.argv[2]) 

prot_ref_map_file = gzip.open(uniprot_uniref_map, 'rb')
sprot_file = gzip.open(sprot_dat, 'rb') 
outFile = './phase_2.tsv'

uniqueMapIds = set()
uniqueSprotIds = set()

print 'stage1'
for line in prot_ref_map_file:
	mappings = line.split('\t')
	uniqueMapIds = uniqueMapIds | {mappings[0]}

print 'stage2'
# Append the SwissProt data that wasn't detected by GO
with open('./phase_1.tsv', 'r') as input_file, open(outFile, 'w') as output_file:
	for line in input_file:
		line = line.replace('\n','')
		elements = line.split('\t')
		if elements[4]=='' and elements[5]=='': # no UniRef/UniProt
			continue
		# Should assume that if there's a UniRef, there's a UniProt
		elif not elements[4]=='': 
			if ',' in elements[4]:
				accs = elements[4].split(',')
				for x in accs:
					uniqueSprotIds = uniqueSprotIds | {x} # Track which SwissProt already present
			else:
				uniqueSprotIds = uniqueSprotIds | {elements[4]} # Track which SwissProt already present
		output_file.write(line + '\n')	

print 'stage3'
# Up til now, building on the assumption that a GO noted accession is present. However,
# need to be able to map those entries which only were found to have evidence through
# SwissProt. This will then exclude columns 2-4. 
with open(outFile, 'a') as output_file:
	for x in uniqueMapIds:
		if x not in uniqueSprotIds:
			output_file.write('UniProtKB:'+protAcc+'\t'+'\t'+'\t'+'\t'+protAcc)
