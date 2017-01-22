# This file contains the functions necessary for processing genes in MetaMerge
# Created by: Leonid Chindelevitch
# Last modified: January 13, 2012

def findEssentialGenes(EssReacts, BoolGenes):
    # Finds the genes essential for at least one of the essential reactions
    Essential = []
    for i in EssReacts:
        if BoolGenes[i]:
            Essential += [x for x in BoolGenes[i][0] if all([x in BoolGenes[i][j] for j in range(1,len(BoolGenes[i]))])]
    return sorted(list(set(Essential)))

def OneOrTwo(Genes1, Genes2):
    # Creates the OR of two lists of genes in conjunctive normal form
    return getUnique(Genes1 + Genes2)

def OneAndTwo(Genes1, Genes2):
    # Creates the AND of two lists of genes in conjunctive normal form          
    return getUnique([OneOrTwo(x,y) for x in Genes1 for y in Genes2])

def findOverlap(Genes1, Genes2):
    # Finds the common genes among two lists of genes
    return sorted(list(set(Genes1).intersection(Genes2)))

def getUnique(List):
    # Returns the unique elements in a list
    if not List:
        return []
    else:
        F = sorted(List)
        return [F[0]] + [F[x] for x in range(1,len(F)) if F[x] != F[x-1]]
