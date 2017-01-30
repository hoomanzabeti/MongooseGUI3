# This file contains functions for identifying matching species or reactions in MetaMerge
# Created by: Leonid Chindelevitch
# Last modified: January 13, 2012

from Utilities import *

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
    if type(first0) != type(first1):
        print(('Error: non-matching feature types!' + str(first0) + str(first1)))
        return
    elif type(first0) == type(0):
        Numbers = True
    elif type(first0) == type(''):
        Strings = True
    elif type(first0) == type([]):
        Lists = True
    elif type(first0) == type({}):
        Dicts = True
    else:
        print('Error: unrecognized feature type!')
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

def supplementFeatures(items, features, default = ''):
    # Adds the default value to a list of feature dictionaries
    # The default may be specified for each feature separately
    if type(default) != type([]):
        default = [default] * len(features)
    newItems = [[]]*len(items)
    for ind0, item in enumerate(items):
        newItems[ind0] = {}
        for ind1, feature in enumerate(features):
            if feature not in item:
                newItems[ind0][feature] = default[ind1]
            else:
                newItems[ind0][feature] = item[feature]
    return newItems

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
        print(('Processing feature number ' + ind))
        curMin = minMatching[ind]
        curFeature1 = [items1[x][features1[ind]] for x in range(M1)]
        curFeature2 = [items2[x][features2[ind]] for x in range(M1)]
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

def reconcileMatches(Matches, scoreVector = [1]):
    # Creates an overall match matrix from a collection of match matrices.
    # Adds a score from a specified score vector for each existing match.
    # By default, the scoreVector has value 1 for each matching feature.
    L = len(Matches)
    m = len(Matches[0])
    n = len(Matches[0][0])
    if len(scoreVector) == 1:
        scoreVector = [1]*L
    Overall = [[]]*m
    for i in range(m):
        Overall[i] = [0]*n
        for j in range(n):
            for k in range(L):
                if Matches[k] and Matches[k][i][j]:
                    Overall[i][j] += scoreVector[k]
    return Overall
