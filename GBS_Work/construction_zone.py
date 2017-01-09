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

#Step 2: Relabel the indexes
def Location_Index_Relabeller(df):
	index = []
	for row in df.iterrows():
		location = row[1][0] + '_' + str(row[1][1])
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

	Allele_Counter_File_Writer(data_dict, "Offsrping_Allele_Counts.csv")
	return data_dict

#Step 2: Determine the top two maximum count reads for both; return in format {LocA: {Allele1: #, Allele2: #}, LocB: {Allele 1: #, Allele 2: #},...}
def Min_Max_Analyzer(data_dict):
	allele_analysis = OrderedDict()
	short_alleles = []

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
			short_alleles.append(key)
		else:
			allele_analysis[key] = frequencies		

	Allele_Counter_File_Writer(allele_analysis, 'Parsed_Allele_Counts.csv')
	
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
def Parent_Info_Stripper(df):
	#Relabel the parent columns for easy access; will need to generalize
	parents = df.ix[:, -4:]
	parent_labels = ['P1-1', 'P1-2', 'P2-1', 'P2-2']
	parents.columns = parent_labels
	
	#Set empty dictionaries
	unique_sites = OrderedDict()
	insertion_sites = OrderedDict()
	identical_sites = OrderedDict()
	total_sites = OrderedDict()
	non_matching_locs = OrderedDict()

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
					unique_sites[location] = (P1_1, P2_1)
					total_sites[location] = (P1_1, P2_1)
				else:
					insertion_sites[location] = (P1_1, P2_1)
					total_sites[location] = (P1_1, P2_1)
			else:
				identical_sites[location] = (P1_1, P2_1)
		else:
			non_matching_locs[location] = ((P1_1, P1_2), (P2_1, P2_2))

	#Write unique sites, insertion sites, and non-matching sites to csv files
	Site_File_Writer(non_matching_locs, "Non_Matching.csv")
	Site_File_Writer(identical_sites, "Identical_Sites.csv")
	Site_File_Writer(unique_sites, 'Unique_Sites.csv') 
	Site_File_Writer(insertion_sites, 'Insertion_Sites.csv')

	return total_sites

#Step 2: Write the parent read info into separated csv files for further by-hand analysis
def Site_File_Writer(site_dict, file_name):
	with open(file_name, 'w') as f:
		for loc in site_dict:
			f.write(loc + ' ,' + str(site_dict[loc][0]).strip('(').strip(')') + ' ,' + str(site_dict[loc][1]).strip('(').strip(')'))
			f.write('\n')

#Step 3: Print out some statistics about the parent reads
def Parent_Stat_Reader(df):
	parents = df.ix[:, -4:]
	parent_labels = ['P1-1', 'P1-2', 'P2-1', 'P2-2']
	parents.columns = parent_labels
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
	print('Unique; same length: %s' %i)
	print('Unique; insertion/deletion: %s' %j)
	print('Identical: %s' %h)
	print('Non-identical P1 Reads: %s' %l)
	print('Non-identical P2 Reads: %s' %m)
	print('Total: %s' %k)
	
	#Print a check sum; should equal the value under Total, otherwise something went wrong
	check_sum = i + j + h + l + m
	print("Check Sum: %s" %check_sum)

#Step 4: Combine the unique and insertion sites into a single dictonary
def Dictionary_Combiner(d1, d2):
	total_dict = {**d1, **d2}
	return total_dict

#PART 4: ASSIGN THE READS TO PARENTS
#Step 1: Assign the top two frequency reads to either parent A or parent B
def Allele_Assigner(allele_dict, parents_dict):
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

#Step 2: Assign the offspring reads to either parent A or parent B
def Read_Assignment(df, assigned_alleles):
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




file = 'test.tab'
df = File_Reader(file)
df = Location_Index_Relabeller(df)
allele_count = Allele_Counter(df)
min_max = Min_Max_Analyzer(allele_count)
total_sites = Parent_Info_Stripper(df)
#Parent_Stat_Reader(df)
#assigned_alleles = Allele_Assigner(min_max, total_sites)
#df = Read_Assignment(df, assigned_alleles)
#print(df)


#print(single_sites)
#print()