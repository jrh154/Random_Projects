import pandas as pd
from collections import OrderedDict
import GBS_Analysis_Functions

file = './maf01_dp6_mm08.tab'

#Reads a file, searches for parents (non-general right now), and returns the read stats (i.e., unique, non-unique, etc)
def Parent_Stat_Reader(file):
	#Read the file with the parent reads in, parse the parents out, and then relabel the columns to useful values
	df = GBS_Analysis_Functions.File_Reader(file)
	locations = df.ix[:, 0:2]
	parents = df.ix[:,-5:-1]
	parents = pd.concat([locations,parents], axis = 1)
	parent_labels = ['Chromosome', 'Pos', 'P1-1', 'P1-2', 'P2-1', 'P2-2']
	parents.columns = parent_labels
	print(parents)
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


df = GBS_Analysis_Functions.File_Reader(file)
df1 = df.ix[:, 0:2]
df2 = df.ix[:, -5:-1]
df3 = pd.concat([df1, df2], axis = 1)

print(df3)