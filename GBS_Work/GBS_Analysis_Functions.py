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
import sys
import time

#Reads in a tab file and returns it as a DataFrame object; relabels columns as well to remove trash
#(Currently relabeller is set for Chickpea format data; need to make more general in future)
def File_Reader(file):
	print("Reading File")

	df = pd.read_csv(file, delimiter = '\t')
	column_labels = []
	for entry in df.columns:
		column_labels.append(entry[4:17])
	df.columns = column_labels
	return df

#Renames the index of the dataframe to be the location in the format Chromosome_Position; also removes the location and reference columns
#This needs to be made more efficienty (currently responsible for 1/2 run time)
def Location_Relabeller(df):
	print("relabelling sites")
	index = []
	for i in range(0, len(df.index)):
		location = str(df.ix[i,0]) + '_' + str(df.ix[i,1])
		if location not in index:
			index.append(location)
		else:
			location += '_1'
			index.append(location)
	df.index = index
	df = df.drop(df.columns[0:3], axis=1)
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
	print("Counting Alleles")
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
def Parent_Info_Stripper(df):
	print("Stripping Parents...not as dirty as it sounds")
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
	print("Here's the stats on the parent reads:")
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
	print("Picking the unique sites")
	unique_sites = OrderedDict()
	insertion_sites = OrderedDict()
	for row in parents.iterrows():
		#If one read is blank, but the other isn't, use the non-blank read
		if row[1]['P1-1'] == './.' and row[1]['P1-2'] != './.':
			row[1]['P1-1'] = row[1]['P1-2']
		elif row[1]['P1-2'] == './.' and row[1]['P1-1'] != './.':
			row[1]['P1-2'] = row[1]['P1-1']
		if row[1]['P2-1'] == './.' and row[1]['P2-2'] != './.':
			row[1]['P2-1'] = row[1]['P2-2']
		elif row[1]['P2-2'] == './.' and row[1]['P2-1'] != './.':
			row[1]['P2-2'] = row[1]['P2-1']

		#Deals with case where both reads are blank for a parent (thus making a false non-read)
		if row[1]['P1-1'] == './.' and row[1]['P1-2'] == './.':
			row[1]['P1-1'] = 'NoRead'
		if row[1]['P2-1'] == './.' and row[1]['P2-2'] == './.':
			row[1]['P2-1'] = 'NoRead'
		
		#Do the sorting into unique and insertion/deletion sites (may not be necessary?)
		if row[1]['P1-1'] == row[1]['P1-2'] and row[1]['P2-1'] == row[1]['P2-2']:
			key = str(row[1]['Chromosome']) + '_' + str(row[1]['Pos']) 
			if row[1]['P1-1'] != row[1]['P2-1']:
				#Parent reads are the same lenght (usually SNP)
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
	print("Min/Max Analysis")
	allele_analysis = OrderedDict()
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
	print("Parent Assignment: A Jerry Springer Episode")
	assigned_allele = OrderedDict()
	for key in parents_dict:
		holder = {}
		#Set allele values for convienent referencing
		allele_1 = list(allele_dict[key].keys())[0]
		allele_2 = list(allele_dict[key].keys())[1]
		
		#Compare and assign allele 1
		if allele_1 == parents_dict[key][0]:
			holder['A'] = allele_1
		elif allele_1 == parents_dict[key][1]:
			holder['B'] = allele_1

		#Compare and assign allele 2
		if allele_2 == parents_dict[key][0]:
			holder['A'] = allele_2
		elif allele_2 == parents_dict[key][1]:
 			holder['B'] = allele_2

 		#In case an allele can't be assigned (shouldn't happen, but just in case)
		if allele_1 != parents_dict[key][0] and allele_1 != parents_dict[key][1]:
			print("Whoops")
			print(key)
		if allele_2 != parents_dict[key][0] and allele_2 != parents_dict[key][1]:
			print("whoops")
			print(key)
		#Add results to output dictionary
		
		if len(holder.keys()) == 2:
			assigned_allele[key] = holder
		else:
			print("Error, either too many or too few keys")
			print("Location: %s, Alleles: %s and %s" %(key, allele_1, allele_2))				
	return assigned_allele

#Assigns a parent to the read from each offspring at each location (i.e., Loc 1, Offspring 1 was C/G, now is A, to show its from parent A)
#Returns the result as a dataframe
def Parent_Assignment(df, assigned_alleles):
	print("And the father is...Assigning parents to each read")
	#Remove any locations where there is no corresponding assigned allele (likely because parents didn't match)
	for location in df.index:
		if location not in assigned_alleles.keys():
			df = df.drop(location, axis = 0)
	
	#Add counters for tracking program completion
	total = len(assigned_alleles.keys())
	i = 1

	#Assign to parent; if it doesn't match either parent option, it is set to unassigned; skips cases where there was no initial read
	for key in assigned_alleles:
		#Allow for monitoring when processing large files
		print("Assigning reads from %s, location %s of %s" % (key, i, total))
		print(key)
		i += 1

		#The heart of the program; goes by rows
		for j in range(0, len(df.columns)):
			allele = df.ix[key,j]
			df.loc[key,df.columns[j]]
			if assigned_alleles[key]['A'] == allele:
				df.ix[key,j] = 'A'
			elif assigned_alleles[key]['B'] == allele:
				df.ix[key,j] = 'B'
			elif allele == './.':
				pass
			else:
				df.ix[key,j] = 'Unassigned'
	return df

#Takes two dictionaries and combines them into a third dictionary; this third dictionary is returned
#NOTE: Only works with Python 3 or greater (syntax not valid in Python 2.7)
def Dictionary_Combiner(d1, d2):
	print("Combining into one unique dictionary")
	total_dict = {**d1, **d2}
	return total_dict



#Run and test functions here
file = './maf01_dp6_mm08.tab'
df = File_Reader(file)
time.clock()
df1 = Location_Relabeller(df)
print(time.clock())

data_dict = Allele_Counter(df)
allele_dict = Min_Max_Analyzer(data_dict)
parents = Parent_Info_Stripper(df)
Parent_Stat_Reader(parents)
unique_dict, insertion_dict = Unique_Site_Stripper(parents)
all_sites = Dictionary_Combiner(unique_dict, insertion_dict)
assigned_alleles = Allele_Assigner(allele_dict, all_sites)
results = Parent_Assignment(df1, assigned_alleles)

print("Saving results")
results.to_csv('test.csv')
