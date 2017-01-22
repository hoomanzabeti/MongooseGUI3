# This file contains functions for processing the metabolite and reaction match matrices
# Created by: Leonid Chindelevitch
# Last modified: January 13, 2012

import codecs
from OutputProcessing import CreateSMatrix
from ReactionMatching import *
from Utilities import *

def printMatches(comparedFeatures, lowerBd, upperBd, Reacts0, Reacts1, Metabs0, Metabs1, Filename = ''):
    # Prints out the reactions corresponding to pairs whose comparison value is between the given bounds.
    # Uses the output of compareFeatures and the desired cutoff.
    Inds0, Inds1, Mat = comparedFeatures[0], comparedFeatures[1], comparedFeatures[-1]
    good = [(Inds0[x], Inds1[y]) for x in range(len(Mat)) for y in range(len(Mat[0])) if Mat[x][y] >= lowerBd and Mat[x][y] <= upperBd]
    output= ''
    for matchedPair in good:
        output += '\n'.join([str(x) + ') ' + ReconstructReaction(Reacts0[x], Metabs0, 0) for x in matchedPair[0]])
        output += '\n\n'
        output += '\n'.join([str(y) + ') ' + ReconstructReaction(Reacts1[y], Metabs1, 1) for y in matchedPair[1]])
        output += '\n\n' + '-'*100 + '\n\n'
    if Filename:
        f = codecs.open(Filename, encoding='utf-8', mode='w')
        f.write(output)
        f.close()
    else:
        print(output)
    return good

def printMetabMatches(Mat, lowerBd, upperBd, Metabs0, Metabs1, Filename = ''):
    # Prints out the reactions corresponding to pairs whose comparison value is between the given bounds.
    # Uses the output of compareFeatures and the desired cutoff.
    good = [(x, y) for x in range(len(Mat)) for y in range(len(Mat[0])) if Mat[x][y] >= lowerBd and Mat[x][y] <= upperBd]
    output= ''
    for matchedPair in good:
        (x,y) = matchedPair
        output += '\n' + str(x) + ') ' + Metabs0[x] + " <==> " + str(y) + ') ' + Metabs1[y]
        output += '\n\n'
    if Filename:
        f = codecs.open(Filename, encoding='utf-8', mode='w')
        f.write(output)
        f.close()
    else:
        print(output)
    return good

def expandMatch(Reacts0, Reacts1, Metabs0, Metabs1, curReactMatch, curMetabMatch, Filename = ''):
    # Expands a match by producing a list of candidates to be taken into consideration
    matched0 = [x for x in range(len(Reacts0)) if [z for z in curReactMatch[x] if z > 0]]
    matched1 = [x for x in range(len(Reacts1)) if [z for z in [curReactMatch[y][x] for y in range(len(Reacts0))] if z > 0]]
    toMatch0 = [x for x in range(len(Reacts0)) if x not in matched0]
    toMatch1 = [x for x in range(len(Reacts1)) if x not in matched1]
    MetabMap01 = matrixToDict(curMetabMatch)
    correctMatch = []
    perfectMatch = []
    noMatch = []
    output = ''
    problems = []
    consider = []
    for i in toMatch0:
        for j in toMatch1:
            if curReactMatch[i][j] >= 0: # ignore the "forbidden" negative values!
                A = compareReact(Reacts0[i], Reacts1[j], MetabMap01, [])
                if type(A) != type(0):
                    if len(A[0]) <= 1 and len(A[1]) <= 1:
                            output += (str(i) + ') ' + ReconstructReaction(Reacts0[i], Metabs0, 0, A[0]) + '\n')
                            output += (str(j) + ') ' + ReconstructReaction(Reacts1[j], Metabs1, 1, A[1]) + '\n')
                            output += '\n\n' + '-'*100 + '\n\n'
                            consider.append([i,j])
                else:
                    perfectMatch.append([i,j])
    if Filename:
        f = codecs.open(Filename, encoding='utf-8', mode='w')
        f.write(output)
        f.close()
    else:
        print(output)
    return (perfectMatch, consider)

def maxMatching(MetabMat, Metabs0, Metabs1):
    # Greedily computes the best one-to-one metabolite matching based on a metabolite matching matrix
    # Reorders the metabolites in the second set according to the numbering of the first one.
    # Also returns two lists of unmatched metabolites (possibly empty) - avoids matches with score 0!
    score = 0
    freeMetabs0, freeMetabs1 = [x for x in Metabs0], [x for x in Metabs1]
    pairs = []
    while(freeMetabs0 and freeMetabs1):
        (bestx, besty, bests) = max([(x,y,MetabMat[x][y]) for x in freeMetabs0 for y in freeMetabs1], key = lambda x:x[2])
        score += bests
        if not bests:
            break
        pairs.append([Metabs0.index(bestx), Metabs1.index(besty)])
        freeMetabs0.remove(bestx)
        freeMetabs1.remove(besty)
    pairs = sorted(pairs, key = lambda x:x[0])
    rMetabs0, rMetabs1 = [Metabs0[x[0]] for x in pairs], [Metabs1[x[1]] for x in pairs]
    return (rMetabs0, rMetabs1, freeMetabs0, freeMetabs1, score)

def checkTransformable(MetabMatch, Reacts0, Reacts1):
    # Checks whether any transformable metabolites are part of the same grouping given by a metabolite matching matrix
    # Returns the index in the cover of the violating groups with bad reactions and the network they belong to (0 or 1)
    m0, m1 = getSize(MetabMatch)
    Mat0, Mat1 = CreateSMatrix(Reacts0, m0), CreateSMatrix(Reacts1, m1)
    (Heights, Widths) = RectangleCover(MetabMatch)
    OK = []
    for i in range(len(Heights)):
        group0, group1 = Heights[i], Widths[i]
        if len(group0) > 1:
            allReags0 = sum([[y for y in range(len(Reacts0)) if Mat0[x][y] < 0] for x in group0],[])
            allProds0 = sum([[y for y in range(len(Reacts0)) if Mat0[x][y] > 0] for x in group0],[])
            overlap0  = sorted(list(set(allReags0).intersection(allProds0)))
            if len(overlap0):
                OK.append((i, overlap0, 0))
        if len(group1) > 1:
            allReags1 = sum([[y for y in range(len(Reacts1)) if Mat1[x][y] < 0] for x in group1],[])
            allProds1 = sum([[y for y in range(len(Reacts1)) if Mat1[x][y] > 0] for x in group1],[])
            overlap1  = sorted(list(set(allReags1).intersection(allProds1)))
            if len(overlap1):
                OK.append((i, overlap1, 1))
    return (OK, Heights, Widths)
