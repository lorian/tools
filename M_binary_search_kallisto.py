import os

def test_kallisto(f):
	return os.system('time kallisto index -i ensembl_cdna_test_index {}'.format(' '.join(f)))

def split_files(f):
	split_A = f[0:len(f)/2]
	split_B = f[len(f)/2:]
	return split_A,split_B

bad_files = [f for f in os.listdir('.') if f.endswith('.cdna.all.fa')]

while len(bad_files) > 1:
	split_A,split_B = split_files(bad_files)
	if not test_kallisto(split_A): # runs without error
		bad_files = split_B
	else:
		bad_files = split_A

print bad_files
test_kallisto(bad_files)
