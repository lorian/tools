import os

def test_kallisto(f):
	return os.system('time /home/pmelsted/kallisto/build/src/kallisto mini-index -i ensembl_cdna_test_index {}'.format(' '.join(f)))

def split_files(f):
	split_A = f[0:len(f)/2]
	split_B = f[len(f)/2:]
	return split_A,split_B

bad_files = [f for f in os.listdir('.') if f.endswith('.cdna.all.fa')]
#bad_files.reverse()

cause_errors = []
for f in bad_files:
	if test_kallisto([f]):
		cause_errors.append(f)

print cause_errors


