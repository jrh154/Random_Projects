file = "C:/Users/John/Desktop/CDC_Togo.fa"
out_file = "C:/Users/John/Desktop/test_file.fa"

i = 0
with open(out_file, 'w') as f1:
	with open(file, 'r') as f2:
		for line in f2:
			if i < 500:
				f1.write(line)
			else:
				break
			i+=1