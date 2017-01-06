import pandas as pd
from collections import OrderedDict, Counter

#PART 1: READ THE FILE

#Step 1: Read the file, drop the unmatched column
def File_Reader(file):
	df = pd.read_csv(file, delimiter = '\t')
	df = df.drop(df.columns[-1], axis = 1)
	column_labels = []
	for entry in df.columns:
		column_labels.append(entry[4:17])
	df.columns = column_labels
	return df

#Step 2: Relabel the indexes (time intensive)
def Location_Index_Relabeller(df):
	location_tups = [(df.ix[i,0], df.ix[i,1]) for i in range(0,len(df.index))]
	index = []
	for tup in location_tups:
		location = str(tup[0]) + '_' + str(tup[1])
		if location not in index:
			index.append(location)
		else:
			location += '_1'
			index.append(location)
	df.index = index
	df = df.drop(df.columns[0:3], axis=1)
	return df

#PART 2: READ THE OFFSPRING DATA

#Step 1: Count the alleles
def Allele_Counter(df):
	df = df.ix[:, :-4]
	data_dict = OrderedDict()
	for row in df.iterrows():
		key = row[0]
		allele_list = row[1][:]
		allele_count = dict(Counter(allele_list))
		data_dict[key] = allele_count
	return data_dict

#Step 2: Determine the top two maximum count reads for both 
def Min_Max_Analyzer(data_dict):
	allele_analysis = OrderedDict()
	location_errors = []

	for location in data_dict:
		frequencies = data_dict[location]
		#Filter out './.' entries; alleles into dictionary with format allelle: frequency
		for allele in frequencies:
			if allele == './.':
				del data_dict[location][allele]

		#If more than two frequencies are present, determine the top two frequencies
		while len(frequencies.keys()) > 2:
			counts = list(frequencies.values())
			allele = list(frequencies.keys())
			min_allele = allele[counts.index(min(counts))]
			del frequencies[min_allele]

		#If less than two frequencies are present, return error with location identifier			
		if len(frequencies.keys()) < 2:
			location_errors.append(location)
			print("Error at location %s; not enough alleles" %location)
		allele_analysis[location] = frequencies		
	
	return allele_analysis

#Step 3: Write the allele counts to a separate file
def Allele_Counter_File_Writer(data_dict, file_name):
	with open(file_name, 'w') as f:
		for key in data_dict:
			f.write(key)
			f.write(',')
			
			for allele in data_dict[key]:
				f.write(allele)
				f.write(',')
			
			f.write('\n')
			f.write(' ,')
			
			for allele_numb in data_dict[key]:
				f.write(str(data_dict[key][allele_numb]))
				f.write(',')
			
			f.write('\n')

#PART 3: READ THE PARENT INFO

#Step 1: Read the parent reads and determine if location is SNP, insertion/deletion, non-unique, or different reads
def Parent_Info_Stripper(parents):
	#Relabel the parent columns for easy access; will need to generalize
	parent_labels = ['P1-1', 'P1-2', 'P2-1', 'P2-2']
	parents.columns = parent_labels
	
	#Set empty dictionaries
	unique_sites = OrderedDict()
	insertion_sites = OrderedDict()
	non_matching_locs = []

	#Determine if location is SNP, insertion/deletion, non-unique, or different reads
	for row in parents.iterrows():
		#Assign parent reads to variable to help clean up code
		location = row[0]
		P1_1 = row[1]['P1-1']
		P1_2 = row[1]['P1-2']
		P2_1 = row[1]['P2-1']
		P2_2 = row[1]['P2-2']
		
		#If one read is blank, but the other isn't, use the non-blank read
		if P1_1 == './.' and P1_2 != './.':
			P1_1 = P1_2
		elif P1_2 == './.' and P1_1 != './.':
			P1_2 = P1_1
		if P2_1 == './.' and P2_2 != './.':
			P2_2 = P2_2
		elif P2_2 == './.' and P2_1 != './.':
			P2_1 = P2_1

		#Deals with case where both reads are blank for a parent
		if P1_1 == './.' and P1_2 == './.':
			P1_1 = 'NoRead'
		if P2_1 == './.' and P2_2 == './.':
			P2_1 = 'NoRead'
		
		#Do the sorting into unique and insertion/deletion sites
		if P1_1 == P1_2 and P2_1 == P2_2:
			if P1_1 != P2_1:
				if len(P1_1) == len(P2_1):
					single_sites[location] = (P1_1, P2_1)
				else:
					insertion_sites[location] = (P1_1, P2_1)
		else:
			non_matching_locs.append(location)
	return single_sites, insertion_sites, non_matching_locs


file = 'test.tab'
df = File_Reader(file)
df = Location_Index_Relabeller(df)
allele_count = Allele_Counter(df)
print(allele_count)