# This file contains miscellaneous auxiliary functions for the merging process
# Created by: Leonid Chindelevitch
# Last modified: May 31, 2013

import time, re
from fractions import *

def revDict(Dict):
    # returns a reverse mapping to the one defined by Dict
    rev = {}
    for i  in Dict:
        u = Dict[i]
        if u in rev:
            print(('Error: key ' + str(u) + ' already in mapping!'))
        else:   
            rev[Dict[i]] = i
    return rev

def mergeDict(Dict1, Dict2):
    # returns a dictionary containing the union of key-value pairs
    import copy
    Dict = copy.deepcopy(Dict1)
    for i in Dict2:
        if i not in Dict:
            Dict[i] = Dict2[i]
        else:
            if Dict[i] == Dict2[i]:
                pass
            else:
                print(('Contradiction: ' + str(i) + ' maps to ' + Dict[i] + ' in the first case, but to ' + Dict2[i] + ' in the second case'))
    return Dict

def findCommon(Dict1, Dict2):
    # finds all elements of Dict1 and Dict2 which match in value
    rev1 = revDict(Dict1)
    rev2 = revDict(Dict2)
    common = set(rev1.keys()).intersection(list(rev2.keys()))
    Match = {}
    for i in common:
        Match[rev1[i]] = rev2[i]
    return Match

def makeLists(List, split = False):
    # This function creates lists out of the true entries of a list
    # If split is set to true, splits each string element into words
    L = len(List)
    outList = [[]]*L
    for i in range(L):
        if List[i]:
            if split:
                outList[i] = List[i].split()
            else:
                outList[i] = [List[i]]
    return outList

def pointerWalk(List0, List1):
    # Performs a pointer walk over a pair of lists (sorted in the same direction)
    # Returns the set of all pairs of indices (x,y) such that List0[x] = List1[y]
    L0, L1 = len(List0), len(List1)
    breaks0 = [0] + [x for x in range(1,L0) if List0[x] != List0[x-1]] + [L0]
    breaks1 = [0] + [x for x in range(1,L1) if List1[x] != List1[x-1]] + [L1]
    i, j = 0, 0
    Pairs = []
    while (i < len(breaks0) - 1) and (j < len(breaks1) - 1):
        if List0[breaks0[i]] == List1[breaks1[j]]:
            Pairs += [(x,y) for x in range(breaks0[i], breaks0[i+1]) for y in range(breaks1[j], breaks1[j+1])]
            i += 1
            j += 1
            continue
        elif List0[breaks0[i]] < List1[breaks1[j]]:
            i += 1
            continue
        else: # List0[breaks0[i]] > List1[breaks1[j]]
            j += 1
            continue
    return Pairs

def reconcileLists(List1, List2):
    # Produces a consensus list of features from two incomplete lists (first one more trusted).
    # Returns the first value in positions the lists disagree, the agreed-upon value elsewhere.
    if len(List1) != len(List2):
        print('Error: the lists have incompatible lengths!')
    List = [0]*len(List1)
    for x in range(len(List1)):
        if List1[x]:
            List[x] = List1[x]
        else:
            List[x] = List2[x]
    return List

def getDegrees(Matrix):
    # This function computes the in- and out-degrees of each metabolite.
    m = len(Matrix)
    inDegrees  = [len([x for x in Matrix[i] if x < 0]) for i in range(m)]
    outDegrees = [len([x for x in Matrix[i] if x > 0]) for i in range(m)]
    return (inDegrees, outDegrees)

def matrixToDict(Matrix):
    # Creates a mapping (asymmetric: rows to columns) from a given "matching" matrix
    Dict = {}
    for i in range(len(Matrix)):
        matches = [j for j in range(len(Matrix[i])) if Matrix[i][j] > 0]
        if len(matches) == 1:
            Dict[i] = matches[0]
        elif len(matches) > 1:
            Dict[i] = matches
    return Dict

def createOrdering(duplicates, length):
    # Creates an ordering for a list of given length such that all the duplicates go to the end
    # The assumption is that the duplicate always follows the first copy of the same element (!)
    # For instance, running it on [5,7,9,12] with length 13 gives [0,1,2,3,4,6,8,10,11,5,7,9,12]
    ordering = [x for x in range(length) if x not in duplicates] + duplicates
    permutation = [ordering.index(x) for x in range(length)]
    return ordering, permutation

def recodeList(List, indices):
    # Recodes a list according to a map (vector of indices)
    List = [List[x] for x in indices]
    return List

def recodeMatrix(Matrix, indices1, indices2):
    # Recodes a matrix according to two maps (vectors of indices)
    L1, L2 = len(indices1), len(indices2)
    NewMatrix = [0] * L1
    for ind1 in range(L1):
        NewMatrix[ind1] = [0] * L2
        for ind2 in range(L2):
            NewMatrix[ind1][ind2] = Matrix[indices1[ind1]][indices2[ind2]]
    return NewMatrix

def recodeAndSortPairs(pairs, indices):
    # Recodes the pairs representing a reaction according to a map
    # (vector of indices) and sorts the result by the first element
    recodedPairs = [[indices[pair[0]], pair[1]] for pair in pairs]
    sortedPairs = sorted(recodedPairs, key = lambda x:x[0])
    return sortedPairs

def RectangleCover(Matrix):
    # Greedily covers a matrix (converted to a binary one) with rectangles; returns the resulting subsets in order
    m, n = len(Matrix), len(Matrix[0])
    BMatrix = [[bool(Matrix[x][y]) for y in range(n)] for x in range(m)]
    Widths, Heights = [], []
    S = sum(sum(BMatrix,[]))
    rowSums, colSums  = [len([_f for _f in BMatrix[i] if _f]) for i in range(m)], [len([_f for _f in [BMatrix[i][j] for i in range(m)] if _f]) for j in range(n)]
    while S:
        best = max([(x, rowSums[x]) for x in range(m)] + [(m + y, colSums[y]) for y in range(n)], key = lambda x:x[1])
        L = best[1]
        if best[0] < m: # a row has won
            # Get the width of the rectangle, then extend the height as much as possible
            width  = [x for x in range(n) if BMatrix[best[0]][x]]
            height = [y for y in range(m) if rowSums[y] == L and sum([BMatrix[y][x] for x in width]) == L]
        else: # best[0] >= m, a column has won
            # Get the height of the rectangle, then extend the width as much as possible
            height = [y for y in range(m) if BMatrix[y][best[0]-m]]
            width = [x for x in range(n) if colSums[x] == L and sum([BMatrix[y][x] for y in height]) == L]
        # update the cover information
        Widths.append(width)
        Heights.append(height)
        for x in height:
            for y in width:
                BMatrix[x][y] = False
                rowSums[x] -= 1
                colSums[y] -= 1
                S -= 1
    return (Heights, Widths)

def findRepeats(listOfLists):
    # Finds and returns entries shared by several lists in a given list of lists
    D = {}
    for ind, List in enumerate(listOfLists):
        for item in List:
            myAdd(D, item, ind)
    repeats = []
    for x in list(D.keys()):
        if len(D[x]) != 1:
            repeats.append((x, D[x]))
    return repeats

def end(n):
    # computes the cardinal affix appropriate for n
    if n > 10 and n < 20:
        return 'th'
    elif n % 10 == 1:
        return 'st'
    elif n % 10 == 2:
        return 'nd'
    elif n % 10 == 3:
        return 'rd'
    else:
        return 'th'

def checkTransitive(Matrix):
    # Checks whether the bipartite graph described by the matrix is transitive,
    # in other words, whether it can be covered by DISJOINT bipartite cliques.
    # If the returned set is empty, this means that everything is transitive.
    (Heights, Widths) = RectangleCover(Matrix)
    allHeights, allWidths = sum(Heights,[]), sum(Widths,[])
    M = [[x for x in allHeights if allHeights.count(x) > 1]] + [[y for y in allWidths if allWidths.count(y) > 1]]
    return M

def ConvertTime():
    # converts the current time to the format YYYY-MM-DDTHH:MM:SSZ
    T = time.localtime()
    return '-'.join([str(x) for x in [T.tm_year, T.tm_mon, T.tm_mday]]) + 'T' + ':'.join([str(x) for x in [T.tm_hour, T.tm_min, T.tm_sec]]) + 'Z'

def Pad(integer, length, char = '0'):
    # returns a padded encoding of a given integer
    encoding = str(integer)
    L = len(encoding)
    return char*(length-L) + encoding

def findExactlyOne(listOfStrings, substrings):
    # Returns the index of the unique string in a given list containing the specified substrings.
    # returns None if there are no matches or multiple matches. substrings can be a single string!
    if type(substrings) == type(''):
        substrings = [substrings]
    candidates = [x for x in listOfStrings if any([y.lower() in x.lower() for y in substrings])]
    option = processCandidates(candidates, ', '.join(substrings))
    if option is None:
        return
    else:
        return listOfStrings.index(option)

def processCandidates(candidates, identity):
    # Prompts the user to select a string from a list of candidates based on identity
    # If there is a single candidate or no candidates, no user input is required
    if len(candidates) == 0:
        print(('Warning: no candidate for ' + identity + ' found in the list'))
        return
    elif len(candidates) > 1:
        print(('Warning: multiple candidates for ' + identity + ' found in the list'))
        print((' '.join([str(ind) + ') ' + string for (ind, string) in enumerate(candidates)])))
        opt = int(input("Please select the correct string from the list above and enter its number, or -1 if none \n"))
        if opt == -1:
            return
    else:
        opt = 0
    return candidates[opt]

def DtoC(List):
    # Transforms a list of lists which is an AND of ORs into
    # another list of lists which is an OR of ANDs, or vice versa (De Morgan)!
    n = len(List)
    lens = [len(x) for x in List]
    product = 1
    for i in range(n):
        product = product * lens[i]
    indices = [0]*product
    for i in range(product):
        index = [0]*n
        N = i
        for j in range(n-1,0,-1):
            index[j] = N % lens[j]
            N = int(N / lens[j])
        index[0] = N
        indices[i] = index
    for i in range(product):
        index = indices[i]
        for j in range(len(index)):
            index[j] = List[j][index[j]]
        indices[i] = index
    return indices

def extractUnique(List):
    # This function extracts the unique elements of a given list
    if len(List):
        List = sorted(List)
        goodInds = [0] + [x for x in range(1, len(List)) if List[x] != List[x - 1]]
        List = [List[x] for x in goodInds]
    return List

def convertToFraction(element):
    # This function converts a numeric or a string fraction into a fraction representation.
    if type(element) in [type(Fraction(0)), type(0), type(0)]:
        return Fraction(element)
    if type(element) == type(0.0):
        element = str(element)  # convert to string first!
    if type(element) in [type(''), type('')]:
        string = element
        power = 0
        string = string.lower()
        if '/' in string:
            # NOTE: assume we cannot have both 'e' and '/'
            index = string.index('/')
            num = int(string[:index])
            den = int(string[index+1:])
        else:
            if 'e' in string:
                index = string.index('e')
                power = int(string[index+1:])
                string = string[:index]
            if '.' in string:
                index = string.index('.')
                ndec = len(string[index+1:])
                den = 10**ndec
                num = int(round(float(string)*den))
            else:
                num = int(string)
                den = 1
        if power >= 0:
            num = num*(10**power)
        else:
            den = den*(10**(-power))
        coeff = Fraction(num, den)
        return coeff
    else:
        print(('Unacceptable input type! ' + str(element)))
        return []

def convertTableToDictionary(table, keyColumn = 'name', valueColumn = -1, header = True):
    # Converts a table (possibly incomplete) into a dictionary
    # If header = True, the first row of the table is considered to be the header
    # The keyColumn is either a column name (if header = True) or a column index for keys.
    # The valueColumn is either a column name (if header = True) or a column index for values.
    Dict = {}
    if header:
        headings = table[0]
        table = table[1:]
        if type(keyColumn) == type(''):
            keyColumn = findExactlyOne(headings, keyColumn)
        if type(valueColumn) == type(''):
            valueColumn = findExactlyOne(headings, valueColumn)
        valueName = headings[valueColumn]
    for ind, row in enumerate(table):
        if ind % 100 == 0:
            print(('Processed ' + str(ind) + ' rows'))
        key = row[keyColumn]
        try:
            value = row[valueColumn]
        except:
            value = None
        if header:
            Dict[key] = {}
            Dict[key][valueName] = value
        else:
            Dict[key] = value
    return Dict
            
def augmentDictionaryByValue(Dict1, Dict2):
    # Augments a dictionary of dictionaries
    for key in list(Dict1.keys()):
        value = list(Dict1[key].values())[0]
        if value in Dict2:
            Dict1[key].update(Dict2[value])
    return Dict1

def cleanupTags(page):
    # Takes an HTML page and removes everything between tags
    # Also removes end-of-line and nbsp characters and splits
    # the page into understandable human-readable word units.
    page = page.replace('&nbsp;', ' ').replace('\n', ' ')
    opens  = [x for x in range(len(page)) if page[x] == '<']
    closes = [x for x in range(len(page)) if page[x] == '>']
    L = len(opens) - 1
    cleanPage = sum([page[closes[x] + 1: opens[x+1]].split() for x in range(L)], [])
    return cleanPage

def removeTags(page):
    # Alternative to cleanupTags from http://love-python.blogspot.com/2008/07/strip-html-tags-using-python.html
    page = page.replace('&nbsp;', ' ').replace('\n', ' ')
    p = re.compile(r'<.*?>')
    return p.sub('', page).strip()

def findTrueIndices(vector):
    # Finds the true indices in a given vector
    return [x for x in range(len(vector)) if vector[x]]

def findFalseIndices(vector):
    # Finds the false indices in a given vector
    return [x for x in range(len(vector)) if not vector[x]]

def mapList(items, List):
    # Same as recodeList, but in reverse order of arguments!
    return [List[x] for x in items]

def updateRecord(fullRecord, currentRecord, currentMarked, mark):
    # This function marks all the currently marked records in a given full record
    # by assigning the value of mark to it. WARNING: This function works in place!
    for item in currentMarked:
        fullRecord[currentRecord[item]] = mark
    return

def filterOut(List, removedItems):
    # This function removes the items indexed by the second list from the first one.
    newList = [List[x] for x in range(len(List)) if x not in removedItems]
    return newList

def getSize(Matrix):
    # This function returns the size of a matrix (number of rows, number of columns).
    m = len(Matrix)
    if m > 0:
        n = len(Matrix[0])
    else:
        n = 0
    return (m, n)

def myAdd(dictionary, key, value):
    # This function adds a value to the list corresponding to the key in a dictionary.
    if key in dictionary:
        dictionary[key].append(value)
    else:
        dictionary[key] = [value]
    return

def groupIdentical(List):
    # This function returns a list of lists of indices corresponding to identical elements.
    # For instance, the input [1,2,3,2,1,3,4] should return [[0,4],[1,3],[2,5],[6]] in order.
    F = sorted([(List[x],x) for x in range(len(List))])
    breaks = [0] + [x for x in range(1, len(F)) if F[x-1][0] != F[x][0]] + [len(F)]
    groups = [[F[x][1] for x in range(breaks[i], breaks[i+1])] for i in range(len(breaks)-1)]
    return groups

def compareMatrices(Mat1, Mat2, tol = 1e-10):
    # This function returns a list of differences between two matrices exceeding a tolerance.
    m,n = getSize(Mat1)
    m1,n1 = getSize(Mat2)
    if m != m1 or n != n1:
        print("Error: the matrices don't have the same size!")
        return
    diffs = []
    if Mat1 != Mat2:
        for i in range(m):
            cur1 = Mat1[i]
            cur2 = Mat2[i]
            curDiffs = [abs(cur1[j] - cur2[j]) for j in range(n)]
            diffs += [(i,j,cur1[j],cur2[j]) for j in range(n) if curDiffs[j] > tol]
    return diffs
