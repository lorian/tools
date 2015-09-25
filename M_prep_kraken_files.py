import os

mfa = open('i100.fasta','r')
new_mfa = open('i100_kraken.fasta','w')
for line in mfa:
	if line[0] == '>': # fasta name
		new_mfa.write('>gi|' + line.partition('gi|')[2]) # remove the strain label at the front of the header because kraken is stupid
	else:
		new_mfa.write(line)

new_mfa.close()
