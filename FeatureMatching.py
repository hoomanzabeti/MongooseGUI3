# This file contains functions for identifying matching species or reactions in MetaMerge
# Created by: Leonid Chindelevitch
# Last modified: January 13, 2012

from Utilities import *
import difflib

def matchStrings(Strings1, Strings2, cutoff = 0.6, maxNum = 10):
    # Finds pairs of strings which are above a given similarity threshold,
    # at most maxNum for each string in the first group; returns a matrix
    # whose i,j-th entry contains 1 if and only if strings i and j match.
    L1 = len(Strings1)
    L2 = len(Strings2)
    Matches = [0] * L1
    count = 0
    for ind in range(L1):
        Matches[ind] = [0] * L2
        curString1 = Strings1[ind]
        if count % 10 == 0:
            print(('Processed ' + str(count) + ' metabolites so far'))
        curMatches = set(difflib.get_close_matches(curString1, Strings2, maxNum, cutoff))
        if curMatches:
            for index in range(L2):
                Matches[ind][index] += int(Strings2[index] in curMatches)
        count = count + 1
    return Matches

def compareFeature(List0, List1):
    # Compares two lists of features, assumed to be strings, numbers or lists of the same type.
    # Returns a matrix containing the common features of the lists in the appropriate entries.
    L0 = len(List0)
    L1 = len(List1)
    M = [0]*L0
    for i in range(L0):
        M[i] = [0]*L1
        for j in range(L1):
            M[i][j] = []
    active0, active1 = [x for x in range(L0) if List0[x]], [x for x in range(L1) if List1[x]]
    redList0, redList1 = [_f for _f in List0 if _f], [_f for _f in List1 if _f]
    first0, first1 = redList0[0], redList1[0]
    Numbers, Strings, Lists, Dicts = False, False, False, False
    if type(first0) != type(first1) and not (isinstance(first0, str), isinstance(first1, str)):
        print(('Error: non-matching feature types!' + str(first0) + str(first1)))
        return
    elif type(first0) == type(0):
        Numbers = True
    elif isinstance(first0, str):
        Strings = True
    elif type(first0) == type([]):
        Lists = True
    elif type(first0) == type({}):
        Dicts = True
    else:
        print(('Error: unrecognized feature type!' + str(first0) + str(first1)))
        return
    if Strings: # strings need to be split into words for easier matching
        redList0, redList1 = [x.lower().split() for x in redList0], [x.lower().split() for x in redList1]
        Lists = True
    if Lists: # list features need a special sorting procedure to facilitate finding common elements
        sortList0 = sorted(sum([[(y, i) for y in redList0[i]] for i in range(len(redList0))],[]))
        sortList1 = sorted(sum([[(y, i) for y in redList1[i]] for i in range(len(redList1))],[]))       
    else:
        sortList0 = sorted([(redList0[x], x) for x in range(len(redList0))])
        sortList1 = sorted([(redList1[x], x) for x in range(len(redList1))])
    perm0, perm1 = [x[1] for x in sortList0], [x[1] for x in sortList1]
    active0, active1 = [active0[x] for x in perm0], [active1[x] for x in perm1]
    sortList0, sortList1 = [x[0] for x in sortList0], [x[0] for x in sortList1]
    pairs = pointerWalk(sortList0, sortList1)
    for pair in pairs:
        x, y = pair[0], pair[1]
        M[active0[x]][active1[y]].append((sortList0[x], sortList1[y]))
    return M

def supplementFeatures(items, features, default = ['']):
    # Adds the default value to a list of feature dictionaries
    # The default may be specified for each feature separately
    L = len(features)
    if len(default) == 1:
        default = default * L
    if len(default) != L:
        print('Error: the number of default feature values does not match!')
    for ind in range(L):
        curFeature = features[ind]
        curDefault = default[ind]
        for item in items:
            if curFeature not in item:
                item[curFeature] = curDefault
    return items

def compareFeatures(items1, items2, features1, features2, minMatching = 1):
    # Compares the specified lists of features from two models
    # items1 and items2 are the lists of feature dictionaries
    # features1 and features2 are the strings of feature names
    # minMatching is a vector specifying the smallest number of
    # matches needed to count a pair of features as matching;
    # it is typically best to use 1 or 2 for string features.
    L1, L2 = len(features1), len(features2)
    if L1 != L2:
        print('Error: the numbers of features do not match!')
        return
    if type(minMatching) == type(0):
        minMatching = [minMatching] * L1 
    if len(minMatching) != L1:
        print('Error: the number of cutoffs does not match!')
        return
    M1, M2 = len(items1), len(items2)
    allComparisons = []
    for ind in range(L1):
        print(('Processing features ' + features1[ind] + ' and ' + features2[ind]))
        curMin = minMatching[ind]
        curFeature1 = [items1[x][features1[ind]] for x in range(M1)]
        curFeature2 = [items2[x][features2[ind]] for x in range(M2)]
        curComparedFeatures = compareFeature(curFeature1, curFeature2)
        curComparison = [0] * M1
        for x in range(M1):
            curComparison[x] = [0] * M2
            for y in range(M2):
                if curComparedFeatures[x][y] and len(curComparedFeatures[x][y]) >= curMin:
                    curComparison[x][y] = 1
        allComparisons.append(curComparison)
    # Finally, collect the global information
    M = reconcileMatches(allComparisons)
    return M

def reconcileMatches(Matches, scoreVector = 1):
    # Creates an overall match matrix from a collection of match matrices.
    # Adds a score from a specified score vector for each existing match.
    # By default, the scoreVector has value 1 for each matching feature.
    L = len(Matches)
    m = len(Matches[0])
    n = len(Matches[0][0])
    if type(scoreVector) == type(0):
        scoreVector = [scoreVector]*L
    Overall = [[]]*m
    for i in range(m):
        Overall[i] = [0]*n
        for j in range(n):
            for k in range(L):
                if Matches[k] and Matches[k][i][j]:
                    Overall[i][j] += scoreVector[k]
    return Overall

def addMatrices(Matrix1, Matrix2):
    # Adds two matrices assumed to be rectangular of the same size
    Matrix = [[Matrix1[x][y] + Matrix2[x][y] for y in range(len(Matrix1[0]))] for x in range(len(Matrix1))]
    return Matrix
