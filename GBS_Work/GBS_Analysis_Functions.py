'''
GBS Analysis Functions
Written by J. Hayes, Last updated 1/4/2017

Description:
A set of functions that are useful for analyzing the results of GBS sequence data
Currently includes an allele counter and SNP stripper

Usage:
No set usage yet, just a (hopefully) useful set of scripts

'''


import pandas as pd
from collections import OrderedDict, Counter

#Reads in a tab file and returns it as a DataFrame object; relabels columns as well to remove trash
#(Currently relabeller is set for Chickpea format data; need to make more general in future)
def File_Reader(file):
	df = pd.read_csv(file, delimiter = '\t')
	column_labels = []
	for entry in df.columns:
		column_labels.append(entry[4:17])
	df.columns = column_labels
	return df

#Takes a DataFrame object from File_Reader, and pulls out only the SNP sites
#Returns SNP info in an ordered dictionary
def SNP_Striper(df):
	data_dict = OrderedDict()
	for row in df.iterrows():
		skip = False
		key = str(row[1][0]) + '_' + str(row[1][1])
		
		#Determine if the reference read is an SNP
		if len(row[1][2]) > 1:
			skip = True
		
		#Determine if the other reads are SNPs (omits insertion polymorphisms)
		#Typically these reads are in the form: './.', hence 3 character limit
		for i in range(3,len(row[1])):
			if len(row[1][i]) > 3:
				skip = True

		#If skip has been triggered, do not add read to the data dictionary
		if skip:
			pass
		else:
			data_dict[key] = row[1]

	return data_dict

#Writes results (in ordered dictionary form) to a csv file
def File_Writer(data_dict, file_name):
	df_out = pd.DataFrame.from_dict(data_dict, orient = "index")
	df_out.to_csv(file_name)

#Counts the frequency of each allele to allow for parent analysis; assumes 
def Allele_Counter(df):
	df = df.ix[:, :-5]
	data_dict = OrderedDict()
	for row in df.iterrows():
		key = str(row[1][0]) + '_' + str(row[1][1])
		allele_list = row[1][3:]
		allele_count = dict(Counter(allele_list))
		data_dict[key] = allele_count
	return data_dict

#Writes the frequency count to a separate csv file; format is:
# Location_A Allelle A1 Allele A2 Allele A3
#                 #        #        #
# Location_B Allelle B1 Allele B2 Allele B3
#                 #        #        #
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

#Reads a file, searches for parents (non-general right now), and returns an abridged dataframe
def Parent_Info_Stripper(file):
	df = File_Reader(file)
	locations = df.ix[:, 0:2]
	parents = df.ix[:,-5:-1]
	parents = pd.concat([locations,parents], axis = 1)
	parent_labels = ['Chromosome', 'Pos', 'P1-1', 'P1-2', 'P2-1', 'P2-2']
	parents.columns = parent_labels

	return parents

#Takes a dataframe containing the parent data (from Parent_Info_Stripper) and calculates how many unique site, how many non-unique sites, and how many bad site reads
#are present
def Parent_Stat_Reader(parents):
	#set counters equal to zero
	i=0
	j=0
	h=0
	k=0
	l=0
	m=0

	#read through data
	for row in parents.iterrows():
		k += 1
		#If one parent read is blank and the other has a read; use the parent with a read value
		#Parent 1
		if row[1]['P1-1'] == './.' and row[1]['P1-2'] != './.':
			row[1]['P1-1'] = row[1]['P1-2']
		elif row[1]['P1-2'] == './.' and row[1]['P1-1'] != './.':
			row[1]['P1-2'] = row[1]['P1-1']
		#Parent 2
		if row[1]['P2-1'] == './.' and row[1]['P2-2'] != './.':
			row[1]['P2-1'] = row[1]['P2-2']
		elif row[1]['P2-2'] == './.' and row[1]['P2-1'] != './.':
			row[1]['P2-2'] = row[1]['P2-1']

		#If repeat parent reads do not match, tally here
		if row[1]['P1-1'] != row[1]['P1-2']:
			l += 1
		elif row[1]['P2-1'] != row[1]['P2-2']:
			m += 1

		#If both parents match, start tally here
		if row[1]['P1-1'] == row[1]['P1-2'] and row[1]['P2-1'] == row[1]['P2-2']:
			#If the parent reads are different, start tally
			if row[1]['P1-1'] != row[1]['P2-1']:
				#Parent reads are the same (usually SNP)
				if len(row[1]['P1-1']) == len(row[1]['P2-1']):
					#tup = (row[1]['P1-1'], row[1]['P2-1'])
					#parent_list.append(tup)
					i += 1
				#Parent reads are different (usually insertion/deletion event)
				else:
					#tup = (row[1]['P1-1'], row[1]['P2-1'])
					#insertion_list.append(tup)
					j+=1
			#Parent reads are the same
			else:
				h+=1

	#Print out the results of the analysis
	print('Unique; same length')
	print(i)
	print('Unique; insertion/deletion')
	print(j)
	print('Identical')
	print(h)
	print('Non-identical P1 Reads')
	print(l)
	print('Non-identical P2 Reads')
	print(m)
	print('Total')
	print(k)
	
	#Print a check sum; should equal the value under Total, otherwise something went wrong
	print("Check Sum")
	print(i+j+h+l+m) 

#Reads a data frame containing parent data, and puts the unique sites into two dictionaries: One where they are the same length (generally SNPs), and one where they are
#different lengths (generally insertion/deletion sites); returns two dictionaries, where keys are the location, w/format Chromsome_Pos
def Unique_Site_Stripper(parents):
	unique_sites = OrderedDict()
	insertion_sites = OrderedDict()
	for row in parents.iterrows():
		if row[1]['P1-1'] == './.' and row[1]['P1-2'] != './.':
			row[1]['P1-1'] = row[1]['P1-2']
		elif row[1]['P1-2'] == './.' and row[1]['P1-1'] != './.':
			row[1]['P1-2'] = row[1]['P1-1']

		if row[1]['P1-1'] == row[1]['P1-2'] and row[1]['P2-1'] == row[1]['P2-2']:
			key = str(row[1]['Chromosome']) + '_' + str(row[1]['Pos']) 
			if row[1]['P1-1'] != row[1]['P2-1']:
				#Parent reads are the same (usually SNP)
				if len(row[1]['P1-1']) == len(row[1]['P2-1']):
					tup = (row[1]['P1-1'], row[1]['P2-1'])
					unique_sites[key] = tup
				else:
					tup = (row[1]['P1-1'], row[1]['P2-1'])
					insertion_sites[key] = tup
	return unique_sites, insertion_sites

#Reads the offspring identifiers, identifies the top two maximum values, and then
#returns a dictionary in format: {Location: {Allele 1: #, Allele 2: #}}
def Min_Max_Analyzer(data_dict):
	allele_analysis = {}
	for key in data_dict:
		frequencies = {}
		
		#Filter out './.' entries; alleles into dictionary with format allelle: frequency
		for allele in data_dict[key]:
			if allele != './.':
				if allele not in frequencies.keys():
					frequencies[allele] = data_dict[key][allele]
				else:
					frequencies[allele].append(data_dict[key][allele])

		#If more than two frequencies are present, determine the top two frequencies
		while len(frequencies.keys()) > 2:
			counts = list(frequencies.values())
			allele = list(frequencies.keys())
			min_allele = allele[counts.index(min(counts))]
			del frequencies[min_allele]

		#If less than two frequencies are present, return error with location identifier			
		if len(frequencies.keys()) < 2:
			print("Error at location %s; not enough alleles" %key)
		allele_analysis[key] = frequencies		
	
	return allele_analysis

#Reads in a list of parent alleles and the counter list output frm Min_Max_Analyzer
#Assigns parentage to each based on this
def Allele_Assigner(allele_dict, parents_dict):
	assigned_allele = {}
	for key in parents_dict:
		holder = {}
		#Set allele values for convienent referencing
		allele_1 = list(allele_dict[key].keys())[0]
		allele_2 = list(allele_dict[key].keys())[1]
		
		#Compare and assign allele 1
		if allele_1 == parents_dict[key][0]:
			holder[allele_1] = 'A'
		elif allele_1 == parents_dict[key][1]:
			holder[allele_1] = 'B'

		#Compare and assign allele 2
		if allele_2 == parents_dict[key][0]:
			holder[allele_2] = 'A'
		elif allele_2 == parents_dict[key][1]:
 			holder[allele_2] = 'B'

 		#In case an allele can't be assigned (shouldn't happen, but just in case)
		if allele_1 != parents_dict[key][0] and allele_1 != parents_dict[key][1]:
			print("Whoops")
			print(key)
		if allele_1 != parents_dict[key][0] and allele_1 != parents_dict[key][1]:
			print("whoops")
			print(key)
		#Add results to output dictionary
		assigned_allele[key] = holder	
	return assigned_allele

def Parent_Assignment(df, assigned_alleles):
	print(len(df.index))
	print(len(df.columns))


#Run Functions here
file = './test.tab'

df = File_Reader(file)
data_dict = Allele_Counter(df)
allele_dict = Min_Max_Analyzer(data_dict)
unique_dict, insertion_dict = Unique_Site_Stripper(Parent_Info_Stripper(file))
assigned_alleles = Allele_Assigner(allele_dict, unique_dict)
Parent_Assignment(df, assigned_alleles)

