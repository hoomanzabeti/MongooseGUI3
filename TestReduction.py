# This file contains functions for testing the model reduction process
# Created by: Leonid Chindelevitch
# Last modified: January 30, 2017

TEST_PATH = "/Users/admin/Documents/Mongoose/Repository/trunk/TestFiles"
from MetaMerge import *
os.chdir(TEST_PATH)
Network = ReadASCIIMatrix('TestMatrixExt.txt')
Network = [[Fraction(x) for x in y] for y in Network]
print("There are " + str(getSize(Network)[0]) + " metabolites and " + str(getSize(Network)[1]) + " reactions")
Irrev = ReadASCIIMatrix('TestIrrevExt.txt')[0]
print("There are " + str(len(Irrev)) + " irreversible reactions among them")
output = reduceMatrix(Network, Irrev, Filename = "AlternativeReduction.txt")
res1 = subprocess.call(['diff', 'AlternativeReduction.txt', 'TestReductionExt.txt'])
res2 = subprocess.call(['diff', 'AlternativeReductionFull.txt', 'TestReductionExtFull.txt'])
### If no output is produced by the last two commands, the test is successful.