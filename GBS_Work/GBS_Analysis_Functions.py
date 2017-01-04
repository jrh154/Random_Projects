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
	data_dict = OrderedDict()
	for row in df.iterrows():
		key = str(row[1][0]) + '_' + str(row[1][1])
		allele_list = row[1][3:]
		allele_count = dict(Counter(allele_list))
		data_dict[key] = allele_count

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