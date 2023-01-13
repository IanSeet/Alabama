import os
import math
import glob
import shutil
import sys
scriptDir = os.path.dirname(__file__)
genrNum = int(sys.argv[1])
popSize = int(sys.argv[2])
fitNum = int(sys.argv[3])

folderName = "genr" + str(genrNum)
parentPath = os.path.join(folderName, folderName + "parents")
with open(parentPath, 'r') as f:
	parentLines = f.readlines()
index = 0

protected = []
for line in parentLines:
	words = line.split(" ")
	if words[0] == "index:":
		protected.append(int(words[1]))
protectedSet = set(protected)

#print(protectedSet)
	
for j in range(0, popSize):
	if j not in protectedSet:
		fileName = folderName + "/" + folderName + "_" + str(j)
		for k in range(0, fitNum):
			simName = fileName + "_a" + str(k)
			trajFileName = simName + ".xyz"
			os.remove(trajFileName)
