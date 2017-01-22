# This file contains functions for processing matching reactions in MetaMerge
# Created by: Leonid Chindelevitch
# Last modified: January 13, 2012

def compareReact(Re1, Re2, MetabMap, Exempt):
    # Compares two reactions based on a mapping between metabolites
    # returns 1 or -1 if they are completely matched by metabolites
    # (the sign indicates the relative direction of the reactions);
    # otherwise, returns a list of unmatched metabolites in each one
    # the metabolites with indices contained in Exempt are ignored!

    # Note: now accommodates the case when MetabMap is one-to-many!
    factor = 1
    Metabs1 = [x[0] for x in Re1]
    Coeffs1 = [x[1] for x in Re1]
    Metabs2 = [x[0] for x in Re2]
    Coeffs2 = [x[1] for x in Re2]
    Matched1 = []
    Matched2 = []
    Unmatched1 = []
    Unmatched2 = []
    for i in range(len(Metabs1)):
        item = Metabs1[i]
        if item in MetabMap and item not in Exempt:
            citem = MetabMap[item]
            found = 0
            if type(citem) == type(0):
                if citem in Metabs2:
                    found = 1
            else:   # a list instead of a single element
                citem = set(citem).intersection(Metabs2)
                if len(citem) > 1:
                    pass
                    # print('Contradiction: more than one matching reactant!')
                if len(citem):
                    found = 1
                    citem = (list(citem))[0]
            if found:
                j = Metabs2.index(citem)
                coeff1 = Coeffs1[i]
                coeff2 = Coeffs2[j]
                if factor == 1:
                    if coeff1 == coeff2:
                        Matched1.append(i)
                        Matched2.append(j)
                    elif coeff1 == -coeff2:
                        if not (Matched1 or Matched2):
                            factor = -1
                            Matched1.append(i)
                            Matched2.append(j)
                        else:
                            pass
                            # print('Contradiction: some reactants have opposite signs!')
                    else:
                        pass
                        # print('Contradiction: incompatible coefficients!')
                elif factor == -1:
                    if coeff1 == -coeff2:
                        Matched1.append(i)
                        Matched2.append(j)
                    elif coeff1 == coeff2:
                        pass
                        # print('Contradiction: some reactants have opposite signs!')
                    else:
                        pass
                        # print('Contradiction: incompatible coefficients!')
                else:
                    print('This should never happen during compareReact!')
            else:
                Unmatched1.append(i)
        elif item not in Exempt:
            Unmatched1.append(i)
        else:
            pass
    Unmatched2 = [x for x in range(len(Metabs2)) if x not in Matched2]
    if not (Unmatched1 or Unmatched2):
        return factor
    else:
        return (Unmatched1, Unmatched2)
        
def addReacts(ReactList):
    # creates a combined reaction representing the sum of those in a given list
    Metabs = []
    Coeffs = []
    Output = []
    for React in ReactList:
        for item in React:
            if item[0] in Metabs:
                ind = Metabs.index(item[0])
                Coeffs[ind] += item[1]
            else:
                Metabs.append(item[0])
                Coeffs.append(item[1])
    # postprocessing
    for i in range(len(Metabs)):
        # ignore 0 coefficients; append all the reactants first, then products
        if Coeffs[i] < 0:
            Output = [[Metabs[i], Coeffs[i]]] + Output
        elif Coeffs[i] > 0:
            Output = Output + [[Metabs[i], Coeffs[i]]]
    return Output

def mapReact(React, metabMap, exempt):
    # Maps a reaction through a metabolite map; skips exempt metabolites
    return [[metabMap[x[0]],x[1]] for x in React if x[0] not in exempt]

def ReconstructReaction(React, Metabs, option = 0, mark = [], number = False, numbering = []):
    # Outputs a reaction based on a global list of metabolites;
    # the option also indicates the style of the output (JP/GB)
    # The metabolites at positions in mark are highlighted.
    # If number = True, the metabolites are also numbered.
    # By default, the numbering is from 0 increasing by 1, but
    # a list of numbers can be supplied as an argument instead.
    if option:
        nRound = 3
        connect = ' == '
    else:
        nRound = 1
        connect = ' <==> '
    L = len(React)
    Mets   = [Metabs[x[0]] for x in React] 
    Coeffs = [round(x[1], nRound) for x in React]      
    fullReact = ['<< '*int(x in mark) + '('*option + str(abs(Coeffs[x])) + ')'*option + ' ' + Mets[x] + ' >>'*int(x in mark) for x in range(L)]
    if number:
        if not numbering:
            numbering = list(range(L))
        fullReact = ['[' + str(numbering[x]) + '] ' + fullReact[x] for x in range(len(fullReact))]
    Reactants = [fullReact[x] for x in range(L) if Coeffs[x] < 0]
    Products  = [fullReact[x] for x in range(L) if Coeffs[x] > 0]
    reaction = ' + '.join(Reactants) + connect + ' + '.join(Products)
    return reaction

def CombineReacts(react0, react1, opposite = False):
    # Takes the union of the metabolites from two reactions, producing a warning if some coefficients are mismatched
    # If opposite is True, the function assumes that the two reactions are oppositely directed; otherwise, co-directed
    if opposite: # reverse the second reaction before matching
        react1 = [[x[0], -x[1]] for x in react1]
    Metabs0, Coeffs0 = [x[0] for x in react0], [x[1] for x in react0]
    Metabs1, Coeffs1 = [x[0] for x in react1], [x[1] for x in react1]
    newReact = [x for x in react0]
    for i in range(len(Metabs1)):
        met1 = Metabs1[i]
        coeff1 = Coeffs1[i]
        if met1 in Metabs0:
            ind0 = Metabs0.index(met1)
            coeff0 = Coeffs0[ind0]
            if coeff0 != coeff1:
                print(('Warning: the coefficients for metabolite ' + str(met1) + ' differ between ' + str(react0) + ' and ' + str(react1)))
            # newReact[ind0][1] = coeff0 + coeff1
        else:
            newReact.append([met1, coeff1])
    return sorted(newReact, key = lambda x:x[1])
