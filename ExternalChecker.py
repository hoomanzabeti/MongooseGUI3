from ModelProcessing import *

def checkFluxMode(Mode, Network, Irreversible):
    val, vec = computeDistance(Network, Mode, norm = 'inf', Irrev = Irreversible)
    if not val:
        print('The mode is not in the flux cone as specified; returning the closest flux mode')
        return vec
    else:
        print('The mode is in the flux cone!')
        return True

def checkCutSet(Set, Network, Irreversible, Target):
    res = testCutSet(Set, Network, Target, rec = False, I = Irreversible)
    if not res:
        print('The set is not a cut set for the target reaction; computing an augmented set')
        (m, n) = getSize(Network)
        columnIndices = list(range(n))
        NetworkRed = [[Network[x][y] for y in range(n) if y not in Set] for x in range(m)]
        columnIndices = filterOut(columnIndices, Set)
        IrreversibleRed = [columnIndices.index(x) for x in Irreversible if x not in Set]
        TargetRed = columnIndices.index(Target)
        NetworkRec, Rev = reconfigureNetwork(NetworkRed, IrreversibleRed)
        columnIndices += [columnIndices[x] for x in Rev]
        (val, vec) = findMinAdded(NetworkRec, TargetRed, Filename = 'CutsetAugment.lp')
        cutset = [int(x[1:]) for x in list(vec.keys()) if x.startswith('T')]
        bestMCS = extractMinimal(cutset, testCutSet, [NetworkRed, TargetRed, 'CutsetAugment.lp'])
        actualMCS = [columnIndices[x] for x in bestMCS]
        return sorted(Set + actualMCS)
    else:
        print('The set is a cut set for the target reaction')
        return True