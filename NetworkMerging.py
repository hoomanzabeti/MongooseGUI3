# This file contains various approaches to merging two networks
# Created by: Leonid Chindelevitch
# Last modified: January 13, 2012

from Utilities import *
from ReactionMatching import *

def mergeNetworks(Reacts0, Reacts1, ReactEquiv, del0, del1, rev0, rev1, MetabEquiv, nMetab0, nMetab1, ext0, ext1, exempt0, exempt1):
    # This function merges two networks by using a map between their reactions
    (metabMapNew0, metabMapNew1, metabFullMap, externals) = prepareNumbering(MetabEquiv, nMetab0, nMetab1, ext0, ext1, exempt0, exempt1)
    reacts0 = set(range(len(Reacts0)))
    reacts1 = set(range(len(Reacts1)))   
    reactsMapped0 = [x for x in list(ReactEquiv.keys()) if type(x)==type(0)]
    reactsMapped0 += sum([[y for y in x] for x in list(ReactEquiv.keys()) if type(x)!=type(0)],[])
    reactsMapped1 = [x for x in list(ReactEquiv.values()) if type(x)==type(0)]
    reactsMapped1 += sum([[y for y in x] for x in list(ReactEquiv.values()) if type(x)!=type(0)],[])
    reactsUnmapped0 = sorted(list(reacts0.difference(reactsMapped0 + del0)))
    reactsUnmapped1 = sorted(list(reacts1.difference(reactsMapped1 + del1)))
    nEquiv = len(ReactEquiv)
    reactFullMap = {}
    auxMetabMap = {}
    for x in metabFullMap:
        auxMetabMap[x] = x
    leftClasses = list(ReactEquiv.keys())
    rightClasses = list(ReactEquiv.values())
    cnt = 0
    reversible = []
    newReacts = []
    problematic = []
    for i in range(nEquiv):
        curReactLeft = leftClasses[i]
        curReactRight = rightClasses[i]
        (singleL, singleR) = (0,0)
        if type(curReactLeft) == type(0):
            singleL = 1
            curReactLeft = [curReactLeft]
        else:
            curReactLeft = list(curReactLeft)
        if type(curReactRight) == type(0):
            singleR = 1
            curReactRight = [curReactRight]
        else:
            curReactRight = list(curReactRight)
        if (singleL, singleR) == (0,0):
            print(("Error: many-to-many mapping !" + str(i)))
        reactsLeft = [mapReact(Reacts0[x], metabMapNew0, exempt0) for x in curReactLeft]
        reactsRight = [mapReact(Reacts1[x], metabMapNew1, exempt1) for x in curReactRight]
        # eliminate possible duplicates
        reactsLeft = list(set([tuple([tuple(y) for y in x]) for x in reactsLeft]))
        reactsRight = list(set([tuple([tuple(y) for y in x]) for x in reactsRight]))
        # add the reactions together
        reactLeftTotal = addReacts(reactsLeft)
        reactRightTotal = addReacts(reactsRight)
        # compare the total reactions
        result = compareReact(reactLeftTotal, reactRightTotal, auxMetabMap, [])
        if type(result) == type(0):
            if singleL:
                newReacts.append(reactLeftTotal)
            else:
                newReacts.append(reactRightTotal)
            reactFullMap[cnt] = (curReactLeft, curReactRight)
            # determine reversibility
            revL = len([x for x in curReactLeft if x in rev0])
            revR = len([x for x in curReactRight if x in rev1])
            if revL or revR:
                reversible.append(cnt)
            cnt += 1
        else:
            problematic.append((i, result))
    for i in range(len(reactsUnmapped0)):
        curInd = reactsUnmapped0[i]
        curReact = mapReact(Reacts0[curInd], metabMapNew0, exempt0)
        newReacts.append(curReact)
        reactFullMap[cnt] = ([curInd], [])
        if curInd in rev0:
            reversible.append(cnt)
        cnt += 1
    for i in range(len(reactsUnmapped1)):
        curInd = reactsUnmapped1[i]
        curReact = mapReact(Reacts1[curInd], metabMapNew1, exempt1)
        newReacts.append(curReact)
        reactFullMap[cnt] = ([], [curInd])
        if curInd in rev1:
            reversible.append(cnt)
        cnt += 1
    return (newReacts, reactFullMap, reversible, problematic, metabFullMap, externals)

def prepareNumbering(MetabEquiv, nMetab0, nMetab1, externals0, externals1, exempt0 = [], exempt1 = []):
    # This function prepares a global numbering of metabolites
    metabs0 = set(range(nMetab0))
    metabs1 = set(range(nMetab1))
    metabsMapped0 = [x for x in list(MetabEquiv.keys()) if type(x)==type(0)]
    metabsMapped0 += sum([[y for y in x] for x in list(MetabEquiv.keys()) if type(x)!=type(0)],[])
    metabsMapped1 = [x for x in list(MetabEquiv.values()) if type(x)==type(0)]
    metabsMapped1 += sum([[y for y in x] for x in list(MetabEquiv.values()) if type(x)!=type(0)],[])
    metabsUnmapped0 = sorted(list(metabs0.difference(metabsMapped0 + exempt0)))
    metabsUnmapped1 = sorted(list(metabs1.difference(metabsMapped1 + exempt1)))
    metabMapNew0 = {}
    metabMapNew1 = {}
    metabFullMap = {}
    externals = []
    nEquiv = len(MetabEquiv)
    leftClasses = list(MetabEquiv.keys())
    rightClasses = list(MetabEquiv.values())
    cnt = 0
    for i in range(nEquiv):
        (extL, extR) = (0,0)
        # true if at least one element is external
        curLeft = leftClasses[i]
        curRight = rightClasses[i]
        if type(curLeft) == type(0):
            metabMapNew0[curLeft] = cnt
            if curLeft in externals0:
                extL = 1
            curLeft = [curLeft]
        else:
            for x in curLeft:
                metabMapNew0[x] = cnt
                if x in externals0:
                    extL = 1
            curLeft = list(curLeft)
        if type(curRight) == type(0):
            metabMapNew1[curRight] = cnt
            if curRight in externals0:
                extR = 1
            curRight = [curRight]
        else:
            for x in curRight:
                metabMapNew1[x] = cnt
                if x in externals1:
                    extR = 1
            curRight = list(curRight)
        metabFullMap[cnt] = (curLeft, curRight)
        if extL and extR:
            externals.append(cnt)
        cnt += 1
    for i in range(len(metabsUnmapped0)):
        curLeft = metabsUnmapped0[i]
        metabMapNew0[curLeft] = cnt
        metabFullMap[cnt] = ([curLeft], [])
        if curLeft in externals0:
            externals.append(cnt)
        cnt += 1
    for i in range(len(metabsUnmapped1)):
        curRight = metabsUnmapped1[i]
        metabMapNew1[curRight] = cnt
        metabFullMap[cnt] = ([],[curRight])
        if curRight in externals1:
            externals.append(cnt)
        cnt += 1
    return (metabMapNew0, metabMapNew1, metabFullMap, externals)

def MetabNumbering(MetabMatch, externals0, externals1, del0 = [], del1 = []):
    # Prepares a global numbering for metabolites in the merged network
    # externals0 and externals1 are the lists of external metabolites
    # del0 and del1 are sets of metabolites to delete from the network
    # NOTE: The deleted metabolites must NOT be part of the matching!
    # NOTE: The MetabMatch matrix needs to have all the -1 changed to 0
    # ALSO: Assumes that the cover is a transitive cover, and that all
    # external compounds are matched only to each other. Fails if not!
    # ALSO: Assumes that external compounds are numbered consecutively
    # and placed at the end of their metabolite lists. May fail if not!
    nMetab0 = len(MetabMatch)
    nMetab1 = len(MetabMatch[0])
    internals0 = [x for x in range(nMetab0) if x not in externals0]
    internals1 = [x for x in range(nMetab1) if x not in externals1]
    MetabMatchInt = [[MetabMatch[x][y] for y in internals1] for x in internals0]
    MetabMatchExt = [[MetabMatch[x][y] for y in externals1] for x in externals0]
    (Groups0I, Groups1I) = RectangleCover(MetabMatchInt)
    (Groups0E, Groups1E) = RectangleCover(MetabMatchExt)
    # Translate into the original indices!
    Groups0E = [[x + len(internals0) for x in group] for group in Groups0E]
    Groups1E = [[x + len(internals1) for x in group] for group in Groups1E]
    unmatched0 = [x for x in range(nMetab0) if not [_f for _f in MetabMatch[x] if _f]]
    unmatchedInt0 = [x for x in unmatched0 if x in internals0]
    unmatchedExt0 = [x for x in unmatched0 if x in externals0]
    unmatched1 = [y for y in range(nMetab1) if not [_f for _f in [MetabMatch[x][y] for x in range(nMetab0)] if _f]]
    unmatchedInt1 = [x for x in unmatched1 if x in internals1]
    unmatchedExt1 = [x for x in unmatched1 if x in externals1]
    LI, LE = len(Groups0I), len(Groups0E)
    UI0, UI1 = len(unmatchedInt0), len(unmatchedInt1)
    UE0, UE1 = len(unmatchedExt0), len(unmatchedExt1)
    nInternal = LI + UI0 + UI1
    nExternal = LE + UE0 + UE1
    nTotal = nInternal + nExternal - (len(del0) + len(del1))
    Numbering0 = [-1 for i in range(nMetab0)]
    Numbering1 = [-1 for i in range(nMetab1)]
    Numbering = {}
    Total = 0
    for i in range(LI):
        group0, group1 = Groups0I[i], Groups1I[i]
        for j in group0:
            Numbering0[j] = Total
        for j in group1:
            Numbering1[j] = Total
        Numbering[Total] = [(x,0) for x in group0] + [(x,1) for x in group1]
        Total += 1
    for i in range(UI0):
        cur = unmatchedInt0[i]
        if cur not in del0:
            Numbering0[cur] = Total
            Numbering[Total] = [(cur,0)]
            Total += 1
    for i in range(UI1):
        cur = unmatchedInt1[i]
        if cur not in del1:
            Numbering1[cur] = Total
            Numbering[Total] = [(cur,1)]
            Total += 1
    externals = list(range(Total, nTotal))
    for i in range(LE):
        group0, group1 = Groups0E[i], Groups1E[i]
        for j in group0:
            Numbering0[j] = Total
        for j in group1:
            Numbering1[j] = Total
        Numbering[Total] = [(x,0) for x in group0] + [(x,1) for x in group1]
        Total += 1
    for i in range(UE0):
        cur = unmatchedExt0[i]
        if cur not in del0:
            Numbering0[cur] = Total
            Numbering[Total] = [(cur,0)]
            Total += 1
    for i in range(UE1):
        cur = unmatchedExt1[i]
        if cur not in del1:
            Numbering1[cur] = Total
            Numbering[Total] = [(cur,1)]
            Total += 1
    if Total != nTotal:
        print("Error: this should never happen during MetabNumbering!")
    return (Numbering0, Numbering1, Numbering, externals)

def ReactNumbering(ReactMatch, irrev0, irrev1, del0 = [], del1 = [], Opposite = []):
    # Prepares a global numbering for reactions in the merged network
    # irrev0 and irrev1 are the sets of initially irreversible reactions 
    # del0 and del1 are sets of reactions to be deleted from the network
    # Opposite contains matched reactions that are oppositely directed.
    # NOTE: The deleted reactions must NOT also be part of the matching!
    # NOTE: The ReactMatch matrix needs to have all the -1 changed to 0
    # ALSO: Assumes that the cover is a transitive cover. Fails if not!
    # ALSO: Assumes that all the many-to-many matches have been resolved.
    nReact0 = len(ReactMatch)
    nReact1 = len(ReactMatch[0])
    (Groups0, Groups1) = RectangleCover(ReactMatch)
    Numbering0 = [[] for i in range(nReact0)]
    Numbering1 = [[] for i in range(nReact1)]
    Numbering = {}
    irrev = []
    # NOTE: [] in Numbering0 or Numbering1 means reaction is deleted!
    unmatched0 = [x for x in range(nReact0) if not [_f for _f in ReactMatch[x] if _f]]
    unmatched1 = [y for y in range(nReact1) if not [_f for _f in [ReactMatch[x][y] for x in range(nReact0)] if _f]]
    L, U0, U1, D0, D1 = len(Groups0), len(unmatched0), len(unmatched1), len(del0), len(del1)
    lens0, lens1 = [len(x) for x in Groups0], [len(x) for x in Groups1]
    maxes = [max(lens0[x], lens1[x]) for x in range(L)]
    nTotal = sum(maxes) + (U0 - D0) + (U1 - D1)
    newOpposite = []
    Total = 0
    for i in range(L):
        group0, group1 = Groups0[i], Groups1[i]
        rev = False
        if len(group0) == maxes[i]:
            j1 = group1[0]
            if j1 not in irrev1:
                rev = True
            Numbering1[j1] = list(range(Total, Total + maxes[i]))
            for j in group0:
                Numbering0[j] = [Total]
                Numbering[Total] = [(j,0), (j1,1)]
                if j not in irrev0:
                    rev = True
                if [j, j1] in Opposite:
                    rev = True
                    newOpposite.append(Total)
                Total += 1
        else: # len(group1) = maxes[i]
            j0 = group0[0]
            Numbering0[j0] = list(range(Total, Total + maxes[i]))
            if j0 not in irrev0:
                rev = True
            for j in group1:
                Numbering1[j] = [Total]
                Numbering[Total] = [(j0,0), (j,1)]
                if j not in irrev1:
                    rev = True
                if [j0, j] in Opposite:
                    rev = True
                    newOpposite.append(Total)
                Total += 1
        if not rev:
            irrev += list(range(Total - maxes[i], Total))
    for i in range(U0):
        j = unmatched0[i]
        if j not in del0:
            Numbering0[j] = [Total]
            Numbering[Total] = [(j,0)]
            if j in irrev0:
                irrev += [Total]
            Total += 1
    for i in range(U1):
        j = unmatched1[i]
        if j not in del1:
            Numbering1[j] = [Total]
            Numbering[Total] = [(j,1)]
            if j in irrev1:
                irrev += [Total]
            Total += 1
    if Total != nTotal:
        print("Error: this should never happen during ReactNumbering!")
    return (Numbering0, Numbering1, Numbering, irrev, newOpposite)

def RenumberReactions(Reactions, Numbering):
    # This function renumbers the reactions (in pair form) based on a map specified by Numbering
    return [[[Numbering[x[0]], x[1]] for x in react if Numbering[x[0]] != -1] for react in Reactions]

def NetworkMerge(MetabMatch, externals0, externals1, MetabDel0, MetabDel1, React0, React1, ReactMatch, irrev0, irrev1, del0, del1, Opposite):
    # Creates a merged network on the basis of the given information about each network and the reaction and metabolite matching matrices
    (MetNumbering0, MetNumbering1, MetNumbering, externals) = MetabNumbering(MetabMatch, externals0, externals1, MetabDel0, MetabDel1)
    (RxnNumbering0, RxnNumbering1, RxnNumbering, irrev, newOpposite) = ReactNumbering(ReactMatch, irrev0, irrev1, del0, del1, Opposite)
    NewReact0 = RenumberReactions(React0, MetNumbering0)
    NewReact1 = RenumberReactions(React1, MetNumbering1)
    NewReact = [[] for i in range(len(RxnNumbering))]
    for i in range(len(RxnNumbering)):
        reacts = RxnNumbering[i]
        if len(reacts) == 2:
            react0, react1 = reacts[0][0], reacts[1][0]
            # decide whether this is a 1-1 match or a 1-many/many-1 match
            c0, c1 = RxnNumbering0.count(i), RxnNumbering1.count(i)
            if c0 == 1 and c1 == 1: # 1-1 match
                NewReact[i] = CombineReacts(NewReact0[react0], NewReact1[react1], (i in newOpposite))
            elif c0 == 1: # 1-many match; just use the second reaction
                NewReact[i] = NewReact1[react1]
            else:
                NewReact[i] = NewReact0[react0]
        else: # a single reaction...
            if reacts[0][1] == 0: # which comes from network0
                NewReact[i] = NewReact0[reacts[0][0]]
            else: # which comes from network1
                NewReact[i] = NewReact1[reacts[0][0]]
    return (NewReact, RxnNumbering0, RxnNumbering1, RxnNumbering, irrev, MetNumbering0, MetNumbering1, MetNumbering, externals)

def MetabMerge(MetabMatch, externals0, externals1, MetabDel0, MetabDel1, React0, React1, irrev0, irrev1):
    # Creates a merged network based only on the information about matching metabolites (note: also collects together isozymes)
    (MetNumbering0, MetNumbering1, MetNumbering, externals) = MetabNumbering(MetabMatch, externals0, externals1, MetabDel0, MetabDel1)
    nReact0, nReact1 = len(React0), len(React1)
    RxnNumbering0, RxnNumbering1 = [-1 for x in range(nReact0)], [-1 for x in range(nReact1)]
    NewReact0 = RenumberReactions(React0, MetNumbering0)
    NewReact1 = RenumberReactions(React1, MetNumbering1)
    # now find the isozymes!
    NewReact = NewReact0 + NewReact1
    # sort the metabolites in the right order
    NewReact = [sorted(x) for x in NewReact]
    Isozymes = groupIdentical(NewReact)
    # keep one reaction per isozyme set!
    anchors = [x[0] for x in Isozymes]
    first = anchors[0]
    if len(NewReact[first]) == 0:
        anchors = anchors[1:]
    NewReact = [NewReact[x] for x in anchors]
    L = len(NewReact)
    RxnNumbering = [-1 for x in range(L)]
    irrev = list(range(L))
    for i in range(L):
        curIrrev = True
        curGroup = Isozymes[i]
        curL = len(curGroup)
        aux = [int(x >= nReact0) for x in curGroup]
        RxnNumbering[i] = [(curGroup[j] - aux[j]*nReact0, aux[j]) for j in range(curL)]
        for j in range(curL):
            if aux[j] == 0:
                RxnNumbering0[curGroup[j]] = i
                if curGroup[j] not in irrev0:
                    curIrrev = False
            else:
                RxnNumbering1[curGroup[j] - nReact0] = i
                if (curGroup[j] - nReact0) not in irrev1:
                    curIrrev = False
        if not curIrrev:
            irrev.remove(i)
    return (NewReact, RxnNumbering0, RxnNumbering1, RxnNumbering, irrev, MetNumbering0, MetNumbering1, MetNumbering, externals)

def MetabMergeExcept(MetabMatch, externals0, externals1, MetabExcept0, MetabExcept1, React0, React1, irrev0, irrev1):
    # Creates a merged network based only on the information about matching metabolites (note: also collects together isozymes)
    # The difference between this version and the previous version is that we keep ALL the metabolites, but allow to identify
    # as isozymes two reactions which differ only by the presence of "exceptional" metabolites, MetabExcept0/1, respectively.
    # In such a case the reaction with the most exceptional metabolites is taken to be the one representing the isozyme subset.
    (MetNumbering0, MetNumbering1, MetNumbering, externals) = MetabNumbering(MetabMatch, externals0, externals1)
    nReact0, nReact1 = len(React0), len(React1)
    RxnNumbering0, RxnNumbering1 = [-1 for x in range(nReact0)], [-1 for x in range(nReact1)]
    NewReact0 = RenumberReactions(React0, MetNumbering0)
    NewReact1 = RenumberReactions(React1, MetNumbering1)
    # now find the isozymes!
    NewReact = NewReact0 + NewReact1
    # sort the metabolites in the right order
    NewReact = [sorted(x) for x in NewReact]
    Except0 = [MetNumbering0[x] for x in MetabExcept0]
    Except1 = [MetNumbering1[x] for x in MetabExcept1]
    Except = Except0 + Except1
    NewReactNoExcept = [[y for y in x if y[0] not in Except] for x in NewReact]
    Isozymes = groupIdentical(NewReactNoExcept)
    # keep the "longest" reaction from each isozyme set!
    anchors = [max(x, key = lambda y:len(NewReact[y])) for x in Isozymes]
    # check if there are any reactions that only involve exceptional metabolites
    first = anchors[0]
    if len(NewReactNoExcept[first]) == 0:
        # keep all such reactions instead of just one
        anchors = Isozymes[0] + anchors[1:]
    NewReact = [NewReact[x] for x in anchors]
    L = len(NewReact)
    RxnNumbering = [-1 for x in range(L)]
    irrev = list(range(L))
    for i in range(L):
        curIrrev = True
        curGroup = Isozymes[i]
        curL = len(curGroup)
        aux = [int(x >= nReact0) for x in curGroup]
        RxnNumbering[i] = [(curGroup[j] - aux[j]*nReact0, aux[j]) for j in range(curL)]
        for j in range(curL):
            if aux[j] == 0:
                RxnNumbering0[curGroup[j]] = i
                if curGroup[j] not in irrev0:
                    curIrrev = False
            else:
                RxnNumbering1[curGroup[j] - nReact0] = i
                if (curGroup[j] - nReact0) not in irrev1:
                    curIrrev = False
        if not curIrrev:
            irrev.remove(i)
    return (NewReact, RxnNumbering0, RxnNumbering1, RxnNumbering, irrev, MetNumbering0, MetNumbering1, MetNumbering, externals)

def prepareMergedModel(MergedModel, growth, Genes, MetabNames, Externals, ReactFeatures = [], MetabFeatures = [], GeneFeatures = []):
    # This function prepares a merged metabolic model for conversion into an SBML file.
    # MergedModel is the 9-tuple (NewReact, RxnNumbering0, RxnNumbering1, RxnNumbering, irrev, MetNumbering0, MetNumbering1, MetNumbering, externals)
    # The other inputs are lists of size 2, containing the corresponding information for model 0 and model 1, respectively.
    # The exception is growth, which is a tuple whose first element is the index of a growth reaction and the second, that of the chosen model.
    MergedGrowth, MergedGenes, MergedMetabNames, MergedExternal, MergedReactFeatures, MergedMetabFeatures, MergedGeneFeatures = [], [], [], [], [], [], []
    # unpacking
    NewReact, RxnNumbering0, RxnNumbering1, RxnNumbering, irrev, MetNumbering0, MetNumbering1, MetNumbering, externals = MergedModel
    ExtNumbering = [MetNumbering[x] for x in externals]
    MetNumbering = [MetNumbering[x] for x in range(len(MetNumbering)) if x not in externals]
    # repacking
    RxnNumberingFull = [RxnNumbering0, RxnNumbering1]
    MetNumberingFull = [MetNumbering0, MetNumbering1]
    # computing
    MergedGrowth = RxnNumberingFull[0][growth[0]]
    print('Processing genes')
    MergedGenes = [sorted(sum([Genes[x[1]][x[0]] for x in RxnNumbering[y]],[])) for y in range(len(RxnNumbering))]
    print('Processing names')
    MergedMetabNames = [' and '.join(sorted([MetabNames[x[1]][x[0]] for x in MetNumbering[y]])) for y in range(len(MetNumbering))]
    # NOTE: for MergedExternal, we only choose ONE of the possible antecedents
    print('Processing externals')
    MergedExternal = [MetNumberingFull[x[1]][Externals[x[1]][x[0]-len(MetabNames[x[1]])]] for x in [z[0] for z in ExtNumbering]]
    if ReactFeatures:
        print('Processing reactions')
        ReactFeatures0, ReactFeatures1 = ReactFeatures[0], ReactFeatures[1]
        Attributes0, Attributes1 = ReactFeatures0[0], ReactFeatures1[0]
        commonAttributes = [x for x in Attributes0 if x in Attributes1]
        uniqueAttributes0, uniqueAttributes1 = [x for x in Attributes0 if x not in commonAttributes], [x for x in Attributes1 if x not in commonAttributes]
        Attributes =  uniqueAttributes0 + uniqueAttributes1 + commonAttributes
        inds0, inds1 = [Attributes0.index(x) for x in uniqueAttributes0], [Attributes1.index(x) for x in uniqueAttributes1]
        C0, C1 = [Attributes0.index(x) for x in commonAttributes], [Attributes1.index(x) for x in commonAttributes]
        MergedReactFeatures = [[]]*(len(RxnNumbering) + 1)
        MergedReactFeatures[0] = Attributes
        for y in range(len(RxnNumbering)):
            L = RxnNumbering[y]
            L0, L1 = [x[0] for x in L if x[1] == 0], [x[0] for x in L if x[1] == 1]
            MergedReactFeatures[y+1] = [' or '.join([ReactFeatures0[i][j] for i in L0]) for j in inds0] + [' or '.join([ReactFeatures1[i][j] for i in L1]) for j in inds1]
            MergedReactFeatures[y+1] += [' OR '.join(x) for x in zip([' or '.join([ReactFeatures0[i][j] for i in L0]) for j in C0],[' or '.join([ReactFeatures1[i][j] for i in L1]) for j in C1])]
    if MetabFeatures:
        print('Processing metabolites')
        MetabFeatures0, MetabFeatures1 = MetabFeatures[0], MetabFeatures[1]
        Attributes0, Attributes1 = MetabFeatures0[0], MetabFeatures1[0]
        commonAttributes = [x for x in Attributes0 if x in Attributes1]
        uniqueAttributes0, uniqueAttributes1 = [x for x in Attributes0 if x not in commonAttributes], [x for x in Attributes1 if x not in commonAttributes]
        Attributes =  uniqueAttributes0 + uniqueAttributes1 + commonAttributes
        inds0, inds1 = [Attributes0.index(x) for x in uniqueAttributes0], [Attributes1.index(x) for x in uniqueAttributes1]
        C0, C1 = [Attributes0.index(x) for x in commonAttributes], [Attributes1.index(x) for x in commonAttributes]
        MergedMetabFeatures = [[]]*(len(MetNumbering) + 1)
        MergedMetabFeatures[0] = Attributes
        for y in range(len(MetNumbering)):
            L = MetNumbering[y]
            L0, L1 = [x[0] for x in L if x[1] == 0], [x[0] for x in L if x[1] == 1]
            MergedMetabFeatures[y+1] = [' or '.join([MetabFeatures0[i][j] for i in L0]) for j in inds0] + [' or '.join([MetabFeatures1[i][j] for i in L1]) for j in inds1]
            MergedMetabFeatures[y+1] += [' OR '.join(x) for x in zip([' or '.join([MetabFeatures0[i][j] for i in L0]) for j in C0],[' or '.join([MetabFeatures1[i][j] for i in L1]) for j in C1])]
    if GeneFeatures:
        print('Processing gene features')
        GeneFeatures0,  GeneFeatures1  = GeneFeatures[0],  GeneFeatures[1]
        # TO BE COMPLETED LATER!
    return (NewReact, MergedGrowth, irrev, MergedGenes, MergedMetabNames, MergedExternal, MergedReactFeatures, MergedMetabFeatures, MergedGeneFeatures)
