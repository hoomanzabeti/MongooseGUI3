# This file contains functions for analyzing a metabolic network and is unrelated to MetaMerge
# Created by: Leonid Chindelevitch
# Last modified: November 30, 2016

from GeneProcessing import *
from functools import reduce

def classifyExchange(FullNetwork, externalMetabs, Irrev):
    # This function finds all the exchange reactions in a given full network
    # (with the external metabolites included) and classifies them into 3 types:
    # input only, output only, and mixed, and additionally by irreversibility.
    m = len(FullNetwork)
    n = len(FullNetwork[0])
    Exchange = [y for y in range(n) if [_f for _f in [FullNetwork[x][y] for x in externalMetabs] if _f]]
    InputIrr, InputRev, OutputIrr, OutputRev, MixedIrr, MixedRev = [], [], [], [], [], []
    for react in Exchange:
        In, Out = False, False
        for x in externalMetabs:
            if FullNetwork[x][react] < 0:
                In = True
            elif FullNetwork[x][react] > 0:
                Out = True
        if In and Out:
            if react in Irrev:
                MixedIrr.append(react)
            else:
                MixedRev.append(react)
        elif In and (not Out):
            if react in Irrev:
                InputIrr.append(react)
            else:
                InputRev.append(react)
        elif (not In) and Out:
            if react in Irrev:
                OutputIrr.append(react)
            else:
                OutputRev.append(react)
        else:
            print('This should never happen!')
    return (InputIrr, InputRev, OutputIrr, OutputRev, MixedIrr, MixedRev)

def findReachable(Network, Inputs, Allowed):
    # This function finds the reachable reactions and metabolites in a metabolic network
    # given a set of input reactions and a subset of allowed input reactions (all other
    # input reactions are assumed to be blocked). NOTE: Network must be RECONFIGURED!
    # Returns two lists - one for metabolites, one for reactions - with the numbers of
    # steps needed to reach it, or -1 if it cannot be reached at all (-2 if forbidden).
    m = len(Network)
    n = len(Network[0])
    Metabs = [-1 for x in range(m)]
    Reacts = [-1 for x in range(n)]
    Reagents = [set([y for y in range(m) if Network[y][x] < 0]) for x in range(n)]
    Products = [set([y for y in range(m) if Network[y][x] > 0]) for x in range(n)]
    Forbidden = [x for x in Allowed if x not in Inputs]
    for x in Forbidden:
        Reacts[x] = -2
    Iter = 0
    newReacts = Allowed
    while newReacts:
        newMetabs = []
        for x in newReacts:
            Reacts[x] = Iter
            newMetabs += [y for y in Products[x] if Metabs[y] == -1]
        newMetabs = list(set(newMetabs))
        Candidate = []
        for y in newMetabs:
            Metabs[y] = Iter
            Candidate += [x for x in range(n) if Reacts[x] == -1 and y in Reagents[x]]
        Candidate = list(set(Candidate))
        newReacts = [x for x in Candidate if all([(Metabs[y] >= 0) for y in Reagents[x]])]
        Iter += 1
    return (Metabs, Reacts)

def findEssentialTop(Network, Inputs, Allowed, target, option = 'basic'):
    # This function finds the essential metabolites in a metabolic network for a target
    # given a set of input reactions and a subset of allowed input reactions (all other
    # input reactions are assumed to be blocked). NOTE: Network must be RECONFIGURED!
    # The two possible options are 'basic' (return reactions essential for reachability)
    # and 'sum' (return reactions whose deletion increases the number of synthesis steps)
    (M0,R0) = findReachable(Network, Inputs, Allowed)
    Essential = []
    if R0[target] < 0:
        print('Error: the target reaction is not reachable under the given configuration!')
        return Essential
    candidates = [x for x in range(len(R0)) if R0[x] >= 0]
    candidates.remove(target)
    if option == 'sum':
        Comps = [y for y in range(len(M0)) if Network[y][target] < 0]
        Sum = sum([M0[y] for y in Comps])
    for cand in candidates:
        (M,R) = findReachable(Network, Inputs + [cand], Allowed)
        if option == 'basic':
            if R[target] < 0:
                Essential.append(cand)
        elif option == 'sum':
            if R[target] < 0 or sum([M0[y] for y in Comps]) > Sum:
                Essential.append(cand)
        else:
            print('Error: unrecognized option!')
    return Essential
        
def checkPredictions(EssReacts, Genes, GoldEss, GoldNonEss, merged = False, mapping = []):
    # This function computes the number of correct and wrong essentiality predictions
    # EssReacts contains the list of essential reactions, Genes, the corresponding sets
    # of genes in OR of ANDs form (i.e. an OR of sublists, with AND within each sublist)
    # GoldEss and GoldNonEss contain the standards for essential and non-essential genes
    # If merged is True, assumes that Genes is a list of size 2, one for each network;
    # in that case, a mapping of the reactions should be supplied as the last argument.
    if merged:
        EssReactsT = [mapping[x] for x in EssReacts]
        All=[[Genes[x[1]][x[0]] for x in y] for y in EssReactsT]
        AllT=[reduce(OneOrTwo,x) for x in All if x]
        EssGenes = findEssentialGenes(list(range(len(EssReacts))), AllT)
    else:
        EssGenes = findEssentialGenes(EssReacts, Genes)
    if merged: # preliminary transformation
        Genes = Genes[0] + Genes[1]
    allGenes = reduce(lambda u,v:u.union(v), [reduce(lambda y,z:y.union(z), x, set([])) for x in Genes])
    nonEssGenes = allGenes.difference(EssGenes)
    Pos = findOverlap(GoldEss, allGenes)
    Neg = findOverlap(GoldNonEss, allGenes)
    TP = findOverlap(Pos, EssGenes)
    FP = findOverlap(Neg, EssGenes)
    FN = findOverlap(Pos, nonEssGenes)
    TN = findOverlap(Neg, nonEssGenes)
    return (TP, FP, FN, TN)

def checkEssentials(EssGenes, AllGenes, GoldEss, GoldNonEss):
    # This function checks a list of essential genes (subset of all genes) against a gold standard
    nonEssGenes = AllGenes.difference(EssGenes)
    TP = findOverlap(GoldEss,     EssGenes)
    FP = findOverlap(GoldNonEss,  EssGenes)
    FN = findOverlap(GoldEss,     nonEssGenes)
    TN = findOverlap(GoldNonEss,  nonEssGenes)
    return (TP, FP, FN, TN)

def getStats(prediction):
    # This function computes statistics based on a set of predictions returned by checkPredictions
    tp, fp, fn, tn = [len(x) for x in prediction]
    sens = round(1.0*tp/(tp+fn),3)
    spec = round(1.0*tn/(tn+fp),3)
    corr = round(1.0*(tp+tn)/(tp+fn+fp+tn),3)
    return (sens, spec, corr)

def DeletionEffect(Genes, merged = False, mapping = []):
    # This function computes the sets of reactions that get disabled by a gene deletion
    # given each of the boolean expressions over genes catalyzing each network reaction
    # If merged is True, assumes that Genes is a list of size 2, one for each network;
    # in that case, a mapping of the reactions should be supplied as the last argument.
    Dico = {}
    if merged:
        All = [[Genes[x[1]][x[0]] for x in y] for y in mapping]
        Genes =[reduce(OneOrTwo,x,[]) for x in All if x]
    for x in range(len(Genes)):
        for y in findEssentialGenes([x], Genes):
            if y in Dico:
                Dico[y] += [x]
            else:
                Dico[y] = [x]
    return Dico

def findHardCoupled(Matrix):
    # This function finds hard-coupled reaction sets
    import networkx as nx
    m, n = len(Matrix), len(Matrix[0])
    metDegrees = [len([y for y in range(n) if Matrix[x][y]]) for x in range(m)]
    suspects = [x for x in range(m) if metDegrees[x] == 2]
    ranges = [tuple([y for y in range(n) if Matrix[x][y]]) for x in suspects]
    G = nx.Graph()
    G.add_nodes_from(list(range(n)))
    G.add_edges_from(ranges)
    return nx.connected_components(G)