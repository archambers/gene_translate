from ast import literal_eval
#start and end of each exon defined as numerical variable
#negative to account for reverse strand

'''
TO DO: PUT EXON START AND END POINTS IN EXTERNAL documents
PUT FASTAS IN EXTERNAL DOCUMENTS 
GENERALIZE FUNCTIONS
MAKE INTO REPOSITORY!!!

'''

#start and end of utr at beginning of gDNA sequence
five_prime_utr_start = -32224067
five_prime_utr_end = -32223929

#exon start points
exon_01_start = -32223928
exon_02_start = -32223086
exon_03_start = -32222806
exon_04_start = -32221325
exon_05_start = -32220878
exon_06_start = -32220641
exon_07_start = -32220284
exon_08_start = -32219786
exon_09_start = -32218108
exon_10_start = -32217266
exon_11_start = -32217067
exon_12_start = -32215385
exon_13_start = -32214255
exon_14_start = -32213840
exon_15_start = -32213252
exon_16_start = -32212911
exon_17_start = -32212627
exon_18_start = -32210936
exon_19_start = -32204389
exon_20_start = -32203882
exon_21_start = -32202599
exon_22_start = -32201500
exon_23_start = -32201006
exon_24_start = -32199145
exon_25_start = -32198730
exon_26_start = -32198559
exon_27_start = -32197594
exon_28_start = -32197072
exon_29_start = -32196421
exon_30_start = -32196150

gDNA_exons = [
	[32223928, 32223856],
	[32223086, 32223005],
	[32222806, 32222511],
	[32221325, 32220978],
	[32220878, 32220756],
	[32220641, 32220405],
	[32220284, 32220129],
	[32219786, 32219592],
	[32218108, 32217995],
	[32217266, 32217153],
	[32217067, 32216945],
	[32215385, 32215226],
	[32214255, 32214110],
	[32213840, 32213688],
	[32213252, 32213135],
	[32212911, 32212824],
	[32212627, 32212474],
	[32210936, 32210752],
	[32204389, 32204137],
	[32203882, 32203770],
	[32202599, 32202076],
	[32201500, 32201117],
	[32201006, 32200831],
	[32199145, 32198926],
	[32198730, 32198649],
	[32198559, 32198421],
	[32197594, 32197299],
	[32197072, 32196925],
	[32196421, 32196324],
	[32196150, 32195437]
]


#exon end points
exon_01_end = -32223856
exon_02_end = -32223005
exon_03_end = -32222511
exon_04_end = -32220978
exon_05_end = -32220756
exon_06_end = -32220405
exon_07_end = -32220129
exon_08_end = -32219592
exon_09_end = -32217995
exon_10_end = -32217153
exon_11_end = -32216945
exon_12_end = -32215226
exon_13_end = -32214110
exon_14_end = -32213688
exon_15_end = -32213135
exon_16_end = -32212824
exon_17_end = -32212474
exon_18_end = -32210752
exon_19_end = -32204137
exon_20_end = -32203770
exon_21_end = -32202076
exon_22_end = -32201117
exon_23_end = -32200831
exon_24_end = -32198926
exon_25_end = -32198649
exon_26_end = -32198421
exon_27_end = -32197299
exon_28_end = -32196925
exon_29_end = -32196324
exon_30_end = -32195437


#cDNA exon cutoffs
c1_start = 1
c1_end = 73
c2_start = 74
c2_end = 155
c3_start = 156
c3_end = 451
c4_start = 452
c4_end = 799
c5_start = 800
c5_end = 922
c6_start = 923
c6_end = 1159
c7_start = 1160
c7_end = 1315
c8_start = 1316
c8_end = 1510
c9_start = 1511
c9_end = 1624
c10_start = 1625
c10_end = 1738
c11_start = 1739
c11_end = 1861
c12_start = 1862
c12_end = 2021
c13_start = 2022
c13_end = 2167
c14_start = 2168
c14_end = 2320
c15_start = 2321
c15_end = 2438
c16_start = 2439
c16_end = 2526
c17_start = 2527
c17_end = 2680
c18_start = 2681
c18_end = 2865
c19_start = 2866
c19_end = 3118
c20_start = 3119
c20_end = 3231
c21_start = 3232
c21_end = 3755
c22_start = 3756
c22_end = 4139
c23_start = 4140
c23_end = 4315
c24_start = 4316
c24_end = 4535
c25_start = 4536
c25_end = 4617
c26_start = 4618
c26_end = 4756
c27_start = 4757
c27_end = 5052
c28_start = 5053
c28_end = 5200
c29_start = 5201
c29_end = 5298
c30_start = 5299
c30_end = 6012

codon_table = literal_eval(open("codon_table.txt").read())

#print codon_table
#cDNA and genomic DNA sequences as FASTA
Notch4_cDNA_sequence = open("Notch4_cDNA.txt","r").read()
Notch4_gDNA_sequence = open("Notch4_gDNA.txt","r").read()


'''
def translate(sequence):
	peptide = []
	for i in sequence:

'''

def create_cDNA_sequence(sequence):
	cDNA_interpolated = ""
	for i in range(len(gDNA_exons)):
		cDNA_interpolated += Notch4_gDNA_sequence[(norm(gDNA_exons[i][0])):(norm(gDNA_exons[i][1]))+1]
	return cDNA_interpolated



def norm(raw):
	return 32224067 - raw

#returns complement strand to given string
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

#takes gDNA base pair as positive number and returns corresponding exon
def bp_to_exon(number):
	if (number <= -exon_01_start and number >= -exon_01_end):
		return "exon_01"
	elif (number <= -exon_02_start and number >= -exon_02_end):
		return "exon_02"
	elif (number <= -exon_03_start and number >= -exon_03_end):
		return "exon_03"
	elif (number <= -exon_04_start and number >= -exon_04_end):
		return "exon_04"
	elif (number <= -exon_05_start and number >= -exon_05_end):
		return "exon_05"
	elif (number <= -exon_06_start and number >= -exon_06_end):
		return "exon_06"
	elif (number <= -exon_07_start and number >= -exon_07_end):
		return "exon_07"
	elif (number <= -exon_08_start and number >= -exon_08_end):
		return "exon_08"
	elif (number <= -exon_09_start and number >= -exon_09_end):
		return "exon_09"
	elif (number <= -exon_10_start and number >= -exon_10_end):
		return "exon_10"
	elif (number <= -exon_11_start and number >= -exon_11_end):
		return "exon_11"
	elif (number <= -exon_12_start and number >= -exon_12_end):
		return "exon_12"
	elif (number <= -exon_13_start and number >= -exon_13_end):
		return "exon_13"
	elif (number <= -exon_14_start and number >= -exon_14_end):
		return "exon_14"
	elif (number <= -exon_15_start and number >= -exon_15_end):
		return "exon_15"
	elif (number <= -exon_16_start and number >= -exon_16_end):
		return "exon_16"
	elif (number <= -exon_17_start and number >= -exon_17_end):
		return "exon_17"
	elif (number <= -exon_18_start and number >= -exon_18_end):
		return "exon_18"
	elif (number <= -exon_19_start and number >= -exon_19_end):
		return "exon_19" 
	elif (number <= -exon_20_start and number >= -exon_20_end):
		return "exon_20" 
	elif (number <= -exon_21_start and number >= -exon_21_end):
		return "exon_21" 
	elif (number <= -exon_22_start and number >= -exon_22_end):
		return "exon_22" 
	elif (number <= -exon_23_start and number >= -exon_23_end):
		return "exon_23" 
	elif (number <= -exon_24_start and number >= -exon_24_end):
		return "exon_24" 
	elif (number <= -exon_25_start and number >= -exon_25_end):
		return "exon_25" 
	elif (number <= -exon_26_start and number >= -exon_26_end):
		return "exon_26" 
	elif (number <= -exon_27_start and number >= -exon_27_end):
		return "exon_27" 
	elif (number <= -exon_28_start and number >= -exon_28_end):
		return "exon_28" 
	elif (number <= -exon_29_start and number >= -exon_29_end):
		return "exon_29" 
	elif (number <= -exon_30_start and number >= -exon_30_end):
		return "exon_30" 

#takes gDNA base pair as positive number and returns corresponding cDNA base pair and exon as tuple
#depends on bp_to_exon
def gDNA_to_cDNA(number):
	if bp_to_exon(number) == "exon_01":
		return (-exon_01_start - number), bp_to_exon(number)
	elif bp_to_exon(number) == "exon_02":
		return (-exon_02_start - number) + c2_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_03":
		return (-exon_03_start - number) + c3_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_04":
		return (-exon_04_start - number) + c4_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_05":
		return (-exon_05_start - number) + c5_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_06":
		return (-exon_06_start - number) + c6_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_07":
		return (-exon_07_start - number) + c7_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_08":
		return (-exon_08_start - number) + c8_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_09":
		return (-exon_09_start - number) + c9_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_10":
		return (-exon_10_start - number) + c10_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_11":
		return (-exon_11_start - number) + c11_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_12":
		return (-exon_12_start - number) + c12_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_13":
		return (-exon_13_start - number) + c13_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_14":
		return (-exon_14_start - number) + c14_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_15":
		return (-exon_15_start - number) + c15_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_16":
		return (-exon_16_start - number) + c16_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_17":
		return (-exon_17_start - number) + c17_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_18":
		return (-exon_18_start - number) + c18_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_19":
		return (-exon_19_start - number) + c19_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_20":
		return (-exon_20_start - number) + c20_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_21":
		return (-exon_21_start - number) + c21_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_22":
		return (-exon_22_start - number) + c22_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_23":
		return (-exon_23_start - number) + c23_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_24":
		return (-exon_24_start - number) + c24_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_25":
		return (-exon_25_start - number) + c25_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_26":
		return (-exon_26_start - number) + c26_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_27":
		return (-exon_27_start - number) + c27_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_28":
		return (-exon_28_start - number) + c28_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_29":
		return (-exon_29_start - number) + c29_start , bp_to_exon(number)
	elif bp_to_exon(number) == "exon_30":
		return (-exon_30_start - number) + c30_start , bp_to_exon(number)

#takes cDNA base pair and returns exon number
def cDNAbp_to_exon(number):
	if number >= c1_start and number < c2_start:
		return "exon_01"
	elif number >= c2_start and number < c3_start:
		return "exon_02"
	elif number >= c3_start and number < c4_start:
		return "exon_03"
	elif number >= c4_start and number < c5_start:
		return "exon_04"
	elif number >= c5_start and number < c6_start:
		return "exon_05"
	elif number >= c6_start and number < c7_start:
		return "exon_06"
	elif number >= c7_start and number < c8_start:
		return "exon_07"
	elif number >= c8_start and number < c9_start:
		return "exon_08"
	elif number >= c9_start and number < c10_start:
		return "exon_09"
	elif number >= c10_start and number < c11_start:
		return "exon_10"
	elif number >= c11_start and number < c12_start:
		return "exon_11"
	elif number >= c12_start and number < c13_start:
		return "exon_12"
	elif number >= c13_start and number < c14_start:
		return "exon_13"
	elif number >= c14_start and number < c15_start:
		return "exon_14"
	elif number >= c15_start and number < c16_start:
		return "exon_15"
	elif number >= c16_start and number < c17_start:
		return "exon_16"
	elif number >= c17_start and number < c18_start:
		return "exon_17"
	elif number >= c18_start and number < c19_start:
		return "exon_18"
	elif number >= c19_start and number < c20_start:
		return "exon_19"
	elif number >= c20_start and number < c21_start:
		return "exon_20"
	elif number >= c21_start and number < c22_start:
		return "exon_21"
	elif number >= c22_start and number < c23_start:
		return "exon_22"
	elif number >= c23_start and number < c24_start:
		return "exon_23"
	elif number >= c24_start and number < c25_start:
		return "exon_24"
	elif number >= c25_start and number < c26_start:
		return "exon_25"
	elif number >= c26_start and number < c27_start:
		return "exon_26"
	elif number >= c27_start and number < c28_start:
		return "exon_27"
	elif number >= c28_start and number < c29_start:
		return "exon_28"
	elif number >= c29_start and number < c30_start:
		return "exon_29"
	elif number >= c30_start and number <= c30_end:
		return "exon_30"

#takes cDNA base pair and returns gDNA and exon number as tuple
def cDNA_to_gDNA(number):
	if cDNAbp_to_exon(number) == "exon_01":
		return exon_01_start + number - c1_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_02":
		return exon_02_start + number - c2_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_03":
		return exon_03_start + number - c3_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_04":
		return exon_04_start + number - c4_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_05":
		return exon_05_start + number - c5_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_06":
		return exon_06_start + number - c6_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_07":
		return exon_07_start + number - c7_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_08":
		return exon_08_start + number - c8_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_09":
		return exon_09_start + number - c9_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_10":
		return exon_10_start + number - c10_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_11":
		return exon_11_start + number - c11_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_12":
		return exon_12_start + number - c12_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_13":
		return exon_13_start + number - c13_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_14":
		return exon_14_start + number - c14_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_15":
		return exon_15_start + number - c15_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_16":
		return exon_16_start + number - c16_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_17":
		return exon_17_start + number - c17_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_18":
		return exon_18_start + number - c18_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_19":
		return exon_19_start + number - c19_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_20":
		return exon_20_start + number - c20_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_21":
		return exon_21_start + number - c21_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_22":
		return exon_22_start + number - c22_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_23":
		return exon_23_start + number - c23_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_24":
		return exon_24_start + number - c24_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_25":
		return exon_25_start + number - c25_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_26":
		return exon_26_start + number - c26_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_27":
		return exon_27_start + number - c27_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_28":
		return exon_28_start + number - c28_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_29":
		return exon_29_start + number - c29_start, cDNAbp_to_exon(number)
	elif cDNAbp_to_exon(number) == "exon_30":
		return exon_30_start + number - c30_start, cDNAbp_to_exon(number)

#finds which frame the base pair is in, based on cDNA transcript
def find_frame(number):
	codon = ""
	if number % 3 == 0:
		codon += Notch4_cDNA_sequence[number-3]
		codon += Notch4_cDNA_sequence[number-2]
		codon += Notch4_cDNA_sequence[number-1]
	elif number % 3 == 1:
		codon += Notch4_cDNA_sequence[number-1]
		codon += Notch4_cDNA_sequence[number]
		codon += Notch4_cDNA_sequence[number+1]
	elif number % 3 == 2:
		codon += Notch4_cDNA_sequence[number-2]
		codon += Notch4_cDNA_sequence[number-1]
		codon += Notch4_cDNA_sequence[number]
	return codon, number % 3

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

#print substitution_effect(int(raw_input("enter location: ")), raw_input("enter replacement: "))

print create_cDNA_sequence(Notch4_gDNA_sequence) == Notch4_cDNA_sequence

#print gDNA_exons[0][1]
#print cDNA_to_gDNA(int(raw_input("enter cDNA: ")))

