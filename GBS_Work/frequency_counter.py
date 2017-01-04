from collections import Counter, OrderedDict
import pandas as pd

file = './maf01_dp6_mm08.tab'

df = pd.read_csv(file, delimiter = '\t')
data_dict = OrderedDict()

for row in df.iterrows():
	key = str(row[1][0]) + '_' + str(row[1][1])
	allele_list = row[1][3:]
	allele_count = dict(Counter(allele_list))
	data_dict[key] = allele_count




with open('DP6_Allele_Count.csv', 'w') as f:
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

