# Replace spaces in mfa names to _ so entire name is retained

import sys

filename = ""
iterarg = iter(sys.argv)
next(iterarg) #skip name of function
for arg in iterarg:
	if filename == "":
		filename = arg #avoid extra space at beginning
	else:
		filename = filename + " " + arg

mfa = open(filename,'r')
text = ""

for line in mfa:
	if line.startswith('>'):
		text = text + line.replace (" ", "_")
	else:
		text = text + line

fa = open(filename[:filename.rfind('.')] + '_named.fa', 'w')
fa.seek(0)
fa.write(text)
fa.truncate()
fa.close()
