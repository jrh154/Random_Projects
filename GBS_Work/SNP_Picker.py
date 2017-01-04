import pandas as pd
from collections import OrderedDict

file = './maf01_dp4_mm08.tab'
#file = 'test.tab'

df = pd.read_csv(file, delimiter = '\t')
data_dict = OrderedDict()

for row in df.iterrows():
	skip = False
	key = str(row[1][0]) + '_' + str(row[1][1])
	if len(row[1][2]) > 1:
		skip = True
	
	for i in range(3,len(row[1])):
		if len(row[1][i]) > 3:
			skip = True
	
	if skip:
		pass
	else:
		data_dict[key] = row[1]

df_out = pd.DataFrame.from_dict(data_dict, orient = "index")
df_out.to_csv('dp4_SNP.csv')