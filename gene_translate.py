'''
TO DO: CREATE HANDLER FOR USER-FRIENDLINESS!
CREATE CONSERVED DOMAIN DICTIONARY
'''

#allows reading from subdirectories (intention to keep directory together, but for others to use as well)
from os import path
wd = path.dirname(__file__)

#allows dictionary to be read from txt file
from ast import literal_eval

gene = raw_input("gene name: ")

#start of 5' utr at
utr_start = literal_eval(open(wd + gene + "/basics.txt","r").read())

#exons saved as list of tuples
gDNA_exons = literal_eval(open(wd + gene + "/exons.txt","r").read())

#codon table as dictionary in external document
codon_table = literal_eval(open("codon_table.txt","r").read())

#genomic DNA sequence as FASTA in external document
gDNA_sequence = open(wd + gene + "/FASTA.txt","r").read()

if path.isfile(wd + gene + "/regions.txt"):
	regions = literal_eval(open(wd + gene +"/regions.txt","r").read())

plus_strand = gDNA_exons[0][0] < gDNA_exons[0][1]

#creates list of cDNA splice boundaries from exon list
def create_cDNA_exons(gDNA_exons):
	cDNA_exons = [1]
	placeholder = 1
	for i in gDNA_exons:
		difference = (i[0] - i[1]) + 1
		cDNA_exons.append(difference + placeholder)
		placeholder += difference
	return cDNA_exons

#saves splice boundaries to variable
cDNA_exons = create_cDNA_exons(gDNA_exons)

#creates cDNA sequence on the fly from gDNA sequence and exons (just like a real polymerase)
def create_cDNA_sequence(sequence):
	cDNA_sequence = ""
	for i in range(len(gDNA_exons)):
		cDNA_sequence += gDNA_sequence[(norm(gDNA_exons[i][0])):(norm(gDNA_exons[i][1]))+1]
	return cDNA_sequence

#normalizes chromosomal base pair number to start at zero (used by create_cDNA_sequence function)
def norm(chromosome_number):
	if plus_strand:
		return chromosome_number - utr_start
	else:
		return utr_start - chromosome_number

#saves cDNA sequence to variable
cDNA_sequence = create_cDNA_sequence(gDNA_sequence)




#returns complement strand to given string (no real use at the moment)
def get_complement(sequence):
	complement = []
	complement_string = ""
	for bp in sequence:
		if bp == "A":
			complement.append("T")
		elif bp == "T":
			complement.append("A")
		elif bp == "G":
			complement.append("C")
		elif bp == "C":
			complement.append("G")
	for i in complement:
		complement_string += i
	return complement_string


#takes chromosomal base pair number and returns exon number (used by gDNA_to_cDNA function)
##actually returns exon number minus one because of zero indexing
def bp_to_exon(number):
	if plus_strand:
		for i in range(len(gDNA_exons)):
			if number >= gDNA_exons[i][0] and number <= gDNA_exons[i][1]:
				return i
	else:
		for i in range(len(gDNA_exons)):
			if number <= gDNA_exons[i][0] and number >= gDNA_exons[i][1]:
				return i

#takes chromosomal base pair number and returns cDNA base pair number
def gDNA_to_cDNA(number):
	exon_number = bp_to_exon(number)
	if plus_strand:
		return (gDNA_exons[exon_number][0] + number) + cDNA_exons[exon_number], exon_number + 1
	else:
		return (gDNA_exons[exon_number][0] - number) + cDNA_exons[exon_number], exon_number + 1



#takes cDNA base pair number and returns exon number (used by cDNA_to_gDNA function)
def cDNAbp_to_exon(number):
	for i in range(len(cDNA_exons)):
		if number >= cDNA_exons[i] and number < cDNA_exons[i+1]:
			return i

#takes cDNA base pair number and returns chromosomal base pair number
def cDNA_to_gDNA(number):
	exon_number = cDNAbp_to_exon(number)
	return gDNA_exons[exon_number][0] - (number - cDNA_exons[exon_number]), exon_number + 1


#takes amino acid number and returns region
def aa_to_region(number):
	for key in regions:
		if number >= regions[key][0] and number < regions[key][1]:
			return key

#print aa_to_region(int(raw_input("enter amino acid number: ")))

#finds which frame the base pair is in, based on cDNA transcript (used by substitution effect function)
def find_frame(number):
	codon = ""
	if number % 3 == 0:
		codon += cDNA_sequence[number-3]
		codon += cDNA_sequence[number-2]
		codon += cDNA_sequence[number-1]
	elif number % 3 == 1:
		codon += cDNA_sequence[number-1]
		codon += cDNA_sequence[number]
		codon += cDNA_sequence[number+1]
	elif number % 3 == 2:
		codon += cDNA_sequence[number-2]
		codon += cDNA_sequence[number-1]
		codon += cDNA_sequence[number]
	return codon, number % 3

def create_aa_sequence(sequence):
	codon_list = []
	codon = ""
	aa_sequence = ""
	for i in range(len(sequence)):
		if i % 3 == 0 or i % 3 == 1:
			codon += sequence[i]
		elif i % 3 == 2:
			codon += sequence[i]
			codon_list.append(codon)
			codon = ""
	for j in codon_list:
		aa_sequence += codon_table[j]
	return aa_sequence

#takes cDNA base pair number, and letter substituted and gives change
def substitution_effect(number, replacement):
	codon = find_frame(number)
	new_codon = ""
	if codon[1] == 0:
		new_codon += codon[0][0]
		new_codon += codon[0][1]
		new_codon += replacement
	elif codon[1] == 1:
		new_codon += replacement
		new_codon += codon[0][1]
		new_codon += codon[0][2]
	elif codon[1] == 2:
		new_codon += codon[0][0]
		new_codon += replacement
		new_codon += codon[0][2]
	return "original: " + codon[0] + " = " + codon_table[codon[0]] + " | new: " + new_codon + " = " + codon_table[new_codon]


'''
RUNS MISSENSE MUTATION ANALYZER
'''
#print substitution_effect(int(raw_input("enter location: ")), raw_input("enter replacement: "))

'''
TRANSLATES FROM CDNA NUMBER TO CHROMOSOME NUMBER
'''
#print cDNA_to_gDNA(int(raw_input("enter cDNA: ")))

'''
TRANSLATES FROM CHROMOSOME NUMBER TO CDNA NUMBER
'''
#print gDNA_to_cDNA(int(raw_input("enter gDNA: ")))
