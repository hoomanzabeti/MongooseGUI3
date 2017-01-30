Fracs = True
from Utilities import *
from ModelProcessing import findFeasible
from math import gcd

def findUnidirectionalOld(N, Irrev, option = 'null', verbose = False):
    # This function finds all unidirectional (effectively irreversible) reactions.
    # NOTE: It assumes that all the reactions in the network can have nonzero flux;
    # otherwise may incorrectly classify blocked reactions as only negative.
    # The option can be 'null' for nullspace (default) or 'row' for rowspace.
    m, n = getSize(N)
    onlyPos, onlyNeg = [], []
    for i in range(n):
        if verbose and i % 100 == 0:
            print(('Processed ' + str(i) + ' reactions so far'))
        if i not in Irrev:
            # test for effective irreversibility of the reaction
            (val1, vec1) = findFeasible(N, i, Irrev, True,  'sub+' + str(i) + 'sets.lp', option = option)
            if (type(val1) == type([]) and len(val1) == 0): # infeasible
                onlyNeg.append(i)
            else:
                (val0, vec0) = findFeasible(N, i, Irrev, False, 'sub-' + str(i) + 'sets.lp', option = option)
                if (type(val0) == type([]) and len(val0) == 0): # infeasible
                    onlyPos.append(i)
    return (onlyPos, onlyNeg)

def extremeElement(vector, Max = False):
    # This function finds the smallest element of a vector.
    # If Max is True, finds the maximum element instead
    if Fracs:
        if Max:
            return max(vector)
        else:
            return min(vector)
    else:
        curmin = vector[0]
        for i in range(1, len(vector)):
            if (Max and compare(vector[i], curmin) == 1) or (not Max and compare(vector[i], curmin) == -1):
                curmin = vector[i]
    return curmin

def filterList(List):
    # This function removes the zero entries from a list of fractions
    zero = convertToFraction(0)
    myNone = lambda x:compare(x,zero)
    return list(filter(myNone, List))

def representable(frac):
    # This function decides whether a fraction is representable as an exact decimal
    den = getden(frac)
    while den % 2 == 0:
        den /= 2
    while den % 5 == 0:
        den /= 5
    if den == 1:
        return True
    else:
        return False

def getnum(frac):
    # Returns the numerator of a fraction
    if Fracs:
        return frac.numerator
    else:
        return frac[0]

def getden(frac):
    # Returns the denominator of a fraction
    if Fracs:
        return frac.denominator
    else:
        return frac[1]

def makeFloat(frac):
    # This function returns the floating point number corresponding to a fraction
    if Fracs:
        return float(frac)
    else:
        return float(frac[0])/frac[1]

def negate(frac):
    # This function returns the negative of a fractional number
    if Fracs:
        return -frac
    else:
        return [-frac[0], frac[1]]

def absValue(frac):
    # This function returns the absolute value of a fractional number
    if Fracs:
        return abs(frac)
    else:
        return [abs(frac[0]), abs(frac[1])]

def compare(frac0, frac1):
    # This function defines the comparison of two fractional numbers, frac0 & frac1.
    # Returns 1 if frac0 is bigger, -1 if frac1 is bigger, and 0 if they are equal.
    if Fracs:
        #return cmp(frac0, frac1) #replaced with ((frac0 > frac1) - (frac0 < frac1))
        #taken from Python 3 website: https://docs.python.org/3.0/whatsnew/3.0.html
        return ((frac0 > frac1) - (frac0 < frac1))
    else:
        diff = subtract(frac0, frac1)
        if diff[0]*diff[1] > 0:
            return 1
        elif diff[0]*diff[1] < 0:
            return -1
        else: # the numerator is 0
            return 0

def multiply(frac0, frac1):
    # This function defines the multiplication of two fractional numbers.
    if Fracs:
        return frac0*frac1
    else:
        [num0, den0] = frac0
        [num1, den1] = frac1
        num = num0*num1
        den = den0*den1
        gcd = GCD(num, den)
        return [num/gcd, den/gcd]

def add(frac0, frac1):
    # This function defines the addition of two fractional numbers.
    if Fracs:
        return frac0+frac1
    else:
        [num0, den0] = frac0
        [num1, den1] = frac1
        num = num0*den1 + num1*den0
        den = den0*den1
        gcd = GCD(num, den)
        return [num/gcd, den/gcd]

def divide(frac0, frac1):
    # This function defines the division of two fractional numbers.
    if Fracs:
        return frac0/frac1
    else:
        [num0, den0] = frac0
        [num1, den1] = frac1
        num = num0*den1
        den = num1*den0
        gcd = GCD(num, den)
        return [num/gcd, den/gcd]

def subtract(frac0, frac1):
    # This function defines the subtraction of two fractional numbers.
    if Fracs:
        return frac0-frac1
    else:
        [num0, den0] = frac0
        [num1, den1] = frac1
        num = num0*den1 - num1*den0
        den = den0*den1
        gcd = GCD(num, den)
        return [num/gcd, den/gcd]

def GCD(a,b):
    # This function computes the greatest common divisor of two integers
    if Fracs:
        return gcd(a,b)
    else:
        while b:
            a,b = b, a%b
        return a

def dotProduct(List1, List2):
    # This function computes a dot product between two lists of numbers
    if len(List1) != len(List2):
        print('Error: incompatible vectors!')
    else:
        return sum([List1[x] * List2[x] for x in range(len(List1))])

def correlationRank(values1, values2, option = 'Spearman'):
    # This function computes the correlation of the ranks between two variables.
    # The possible options are Spearman and Kendall; both give a value in [-1,1].
    # NOTE: we assume that there are no ties for the current implementation!
    m, n = len(values1), len(values2)
    if m != n:
        print('Error: the lists do not have the same dimension!')
        return
    if option == 'Spearman':
        ranks1, ranks2 = rankList(values1), rankList(values2)
        deltas = [ranks1[x] - ranks2[x] for x in range(n)]
        sumsqrs = float(sum([x*x for x in deltas]))
        return 1 - 6*sumsqrs/(n*(n*n-1))
    elif option == 'Kendall':
        Sum = 0
        for i in range(n):
            for j in range(i+1, n):
               Sum += (compare(values1[i], values1[j]) * compare(values2[i], values2[j]))
        return 2*float(Sum)/(n*(n-1))
    else:
        print('Error: unrecognized option!')
        return

def rankList(values, reverse = False):
    # This function returns the ranks of the values; the smallest gets rank 1.
    # Ties are broken in a stable way; if reverse is True, ranks are reversed.
    myList = [(values[x],x) for x in range(len(values))]
    myList = sorted(myList, reverse = reverse)
    return [x[1] for x in myList]