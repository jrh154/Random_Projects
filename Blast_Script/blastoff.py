from Bio.Blast import NCBIWWW
import pandas as pd
from os.path import join

seq_list = "C:/Users/John/Desktop/scaffold_seqs.csv"
save_dir = "C:/Users/John/Desktop/Blast_Results"

df = pd.read_csv(seq_list, header = None)

for row in df.iterrows():
	seq_loc = row[1][0].replace(":", "_")
	sequence = row[1][1]
	print("Blasting region %s..." %seq_loc)
	blast_result = NCBIWWW.qblast('blastn', 'nt', sequence)
	print("Saving sequence %s..." %seq_loc)
	save_file = open(join(save_dir, seq_loc+'.xml'), 'w')
	save_file.write(blast_result.read())
	save_file.close()
	blast_result.close()