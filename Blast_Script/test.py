import pandas as pd
import numpy as np

file = "C:/Users/John/Desktop/test_set.csv"
reference = "C:/Users/John/Desktop/test_file.fa"
output = "./test_out.csv"

df = pd.read_csv(file, header = None)
location_list = []
for snp_loc in df[0]:
	location = snp_loc.split('_')[0]
	location_list.append(location)
location_list = list(set(location_list))

sequence_dict = {}
with open(reference, 'r') as f:
	record = False
	for line in f:
		line = line.strip('\n')
		if '>' in line:
			record = False
			if line.strip('>') in location_list:
				location = line.strip('>')
				sequence_dict[location] = ''
				record = True
		if record and location not in line:
			sequence_dict[location] += line

with open(output, 'w') as f:
	for location in sequence_dict:
		f.write(location + ', ' + sequence_dict[location] + '\n')