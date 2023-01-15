import filecmp
import os
import sys

fname1 = sys.argv[1]
fname2 = sys.argv[2]
if len(sys.argv) > 3:
	verbose = 1
else:
	verbose = 0
f1 = open(fname1)
f2 = open(fname2)

result = filecmp.cmp(fname1, fname2, shallow=False)
if (result):
	print("pass")
else:
	f1_lines = f1.readlines()
	f2_lines = f2.readlines()

	i = 0

	for i in range (0, len(f1_lines)):
		if f1_lines[i] == f2_lines[i]:
			print("Line ", i, "is identical")
		else:
			print("Line ", i, ":")
			if verbose:
				for x in f1_lines[i]:
					print(x, ord(x))
				print("---")
				for x in f2_lines[i]:
					print(x, ord(x))
			else:
				print(f1_lines[i])
				print(f2_lines[i])

	f1.close()
	f2.close()
	
