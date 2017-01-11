'''
Parse_SNP
Written by J. Hayes, Last Editted 1/11/2017

Purpose: The purpose of this script is to take a list of reads from a GBS analysis
and parse the file so only SNPs remain (no insertions or deletions)

Input: .tab file containing the read information (directly put out by following GBS analysis
directions in Processes.md)

Output: .tab file containing only SNPs; will be output in the same directory as the input file

Usage: python Parse_SNP.py <input tab file> [output file name]
(If no output file name is given, the default will be the input file name + _SNP.tab)

'''

from collections import OrderedDict
from os import mkdir
from os.path import dirname, isfile, isdir
import pandas as pd
import sys

def Parse_SNP(file_in, file_out):
	#Open files
	df = pd.read_csv(file_in, delimiter = '\t')
	data_dict = OrderedDict()

	#Parse through the data file
	for row in df.iterrows():
		key = str(row[1][0]) + '_' + str(row[1][1])
		if len(row[1][2]) == 1:
			SNP = True
		else:
			SNP = False

		if SNP:
			for i in range(3,len(row[1])):
				if len(row[1][i]) != 3:
					SNP = False
					break

		if SNP:
			data_dict[key] = row[1]

	df_out = pd.DataFrame.from_dict(data_dict, orient = "index")
	df_out.to_csv(file_out)


if len(sys.argv) == 3:
	if not isdir(dirname(sys.argv[2])):
		print("The directory %s does not exist" %dirname(sys.argv[2]))
		create_dir = input("Would you like me to create it for you? [Y/N]: ")
		if create_dir == 'Y' or create_dir == 'y' or create_dir == 'Yes' or create_dir == 'yes':
			mkdir(dirname(sys.argv[2]))
		else:
			sys.exit("ERROR: Please try again with a different save name")

	if isfile(sys.argv[1]) and isdir(dirname(sys.argv[2])):
		Parse_SNP(sys.argv[1], sys.argv[2])
	elif not isfile(sys.argv[1]):
		sys.exit("Sorry, I can't find the .tab file to input; please try again")

elif len(sys.argv) == 2:
	if isfile(sys.argv[1]):
		file_out = sys.argv[1][:-4] + '_SNP.tab'
		Parse_SNP(sys.argv[1], file_out)
	else:
		sys.exit("Sorry, I can't find the .tab file to input; please try again")
else:
	print("I'm sorry, I don't understand; I take up to 2 arguments")
	print("Input File (required), output file name (optional)")

