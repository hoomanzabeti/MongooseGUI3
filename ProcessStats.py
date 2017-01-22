import os, itertools
from ModelProcessing import *

def processModelStats(model):
    source = model.description['Source File']
    m1, n1 = getSize(model.fullMatrix)
    m2, n2 = getSize(model.Matrix)
    if n1 != n2:
        print('Error: dimension mismatch; this should never happen!')
    m3, n3 = getSize(model.reducedMatrix)
    allMetStatuses = [x.reductionStatus for x in model.metabolites if not x.external]
    allRxnStatuses = [x.reductionStatus for x in model.reactions]
    TMet  = allMetStatuses.count(1)
    TRxn  = allRxnStatuses.count(1)
    ThRxn = allRxnStatuses.count(2)
    SRxn  = allRxnStatuses.count(3)
    initRev = [i for i,x in enumerate(model.reactions) if x.reversible]
    finalIrrev = sum([[y[0] for y in x.pairs if not x.reversible] for x in model.reactionSubsets],[])
    diffRev = [x for x in initRev if x in finalIrrev]
    negMult = sum([[y[0] for y in x.pairs if y[1] < 0] for x in model.reactionSubsets],[])
    EFor  = len([x for x in diffRev if x not in negMult])
    EBak  = len([x for x in diffRev if x in negMult])
    SubRxn= sum([len(x.pairs) for x in model.reactionSubsets if len(x.pairs) > 1])
    NSub  = len([i for i,x in enumerate(model.reactionSubsets) if len(x.pairs) > 1])
    RedMet= allMetStatuses.count(6)
    curStatus = -1
    curBio = model.findBiomassReaction()
    if curBio != -1:
        curStatus = allRxnStatuses[curBio]
    result = [source, m1, m2, n2, m3, n3, TMet, TRxn, ThRxn, SRxn, EFor, EBak, SubRxn, NSub, RedMet, curStatus]
    return result

def processShelveStats(shelfName):
    s = shelve.open(shelfName)
    Tab = {}
    for key in sorted(s.keys()):
        print(('Processing model ' + key))
        model = s[key]
        Tab[key] = processModelStats(model)
    s.close()
    return Tab

##q = shelve.open('TabulatedResults')
##TabR = q['TabR']
##q.close()
##initDir = os.getcwd()
##myDir = '../GSMNs/NewGSMNs'
##os.chdir(myDir)
##shelfList = ['CobraSBMLs','LastCobraSBMLs', '../../MetaMerge/CorrectedModelsNew', '../../MetaMerge/ReparsedSBMLs', '../../MetaMerge/BoundaryModels', 'NewExcels', 'NewSBMLs', 'LastSBMLs', '../../../../OldNetworks']
##fullTab = {}
##for shelfName in shelfList:
##    print('Processing shelf ' + shelfName)
##    fullTab[shelfName] = processShelveStats(shelfName)
##os.chdir(initDir)
##File = '../../../../SupplementaryData1.txt'
##f = open(File, 'r')
##g = f.readlines()
##f.close()
##Table = [x.strip().split('\t') for x in g][0]
##Inds=[i for i,x in enumerate(Table) if x.startswith('\r') and len(x) > 1]
##Inds=[0] + Inds + [len(Table)]
##TableS=[x.replace('\r','') for x in Table]
##TableR=[TableS[Inds[ind]:Inds[ind+1]] for ind in range(len(Inds)-1)]
##s=shelve.open('WebModels')
##allShelves=[]
##for x in shelfList:
##	allShelves.append(shelve.open(x))
##for key in TabR:
##	cur = TabR[key]
##	for item in cur:
##		relShelf=allShelves[shelfList.index(item[0])]
##		if item[1] in relShelf:
##			print(item[1])
##			relItem = relShelf[item[1]]
##			relKey = item[1]
##			if relKey.endswith('Reduced'):
##				relKey=relKey.replace('Reduced','')
##			if relKey in s:
##				print('Error: duplicating key ' + relKey)
##			else:
##				s[relKey]=relItem
##		else:
##			print('Error: missing key!')
##s.close()

def extractSigns(key):
    filename = key + 'EnergyReductionFull.txt'
    Dict = {}
    if filename in os.listdir('.'):
        print(('Processing model ' + str(key)))
        f = open(filename)
        g = f.readlines()
        f.close()
        G = g[-12].strip()
        if len(G):
            G = G.split('\t')
            L = len(G)
            n = len(G[0])
            if any([len(x) != n for x in G]):
                return Dict
            print(('Found ' + str(L) + ' sign patterns of length ' + str(n)))
            m = n
            while (L % 2 == 0):
                L /= 2
                m -= 1
            target = 2 ** m
            if (target == 1):
                return Dict
            M = max(m,5)
            print(('Processing combinations of size up to ' + str(M)))
            print(('The target is ' + str(target)))
            index = 0
            for k in range(1, M + 1):
                curSize = 2**k
                curOpts = [''.join(x) for x in itertools.product('+-', repeat = k)]
                for subset in generateSubsets(n, k):
                    index += 1
                    if index % 100 == 0:
                        print(('Processed ' + str(index) + ' combinations so far'))
                    if not any([all([x in subset for x in y]) for y in list(Dict.keys())]): # ignore supersets!
                        List = [''.join([x[i] for i in subset]) for x in G]
                        curEntry = dict([(opt, List.count(opt)) for opt in curOpts])
                        if 0 in list(curEntry.values()) and len(set(curEntry.values())) == 2:
                            Dict[tuple(subset)] = curEntry
                            L = (L * curSize) / len([_f for _f in list(curEntry.values()) if _f])
                            if L >= target:
                                return Dict
                            break
    return Dict

def processEBA(listKeys):
    numRE = re.compile('[0-9]+')
    Dict = {}
    for filename in os.listdir('.'):
        if filename.endswith('EnergyReduction.txt') and not filename.endswith('NewEnergyReduction.txt'):
            key = filename[:-len('EnergyReduction.txt')]
            f = open(filename)
            g = f.read()
            f.close()
            numbers = [int(x) for x in re.findall(numRE, g)]
            numbers = numbers[:4]
            Dict[key] = numbers
    return Dict

def processStats(statDict):
    fullStats = {}
    for key in sorted(statDict.keys()):
        print(('Processing model ' + key))
        (m1, n1, m2, n2, m3, n3, TMet, TRxn, ThRxn, SRxn, EFor, EBak, SubRxn, NSub, RedMet, Loops, curStatus, numCorr, minMedia, numEss, numLethal) = statDict[key]
        if (n1 != n2):
            print("Error: the number of reactions in the full and processed matrices should be identical!")
        sizeInit = n2
        sizeFinal = n3
        factRed = float(sizeInit)/sizeFinal
        nMet0 = m2
        nMet1 = nMet0 - TMet
        fracMet1 = 1 - float(nMet1)/nMet0
        nMet2 = nMet1 - RedMet
        fracMet2 = 1 - float(nMet2)/nMet1
        if (nMet2 != m3):
            print("Error: the number of metabolites should be identical!")
        nRxn0 = n2
        nRxn1 = nRxn0 - TRxn
        nRxn2 = nRxn1 - ThRxn
        nRxn3 = nRxn2 - SRxn
        fracRxn1 = 1 - float(nRxn3)/nRxn0
        fracFor  = float(EFor)/nRxn3
        fracBak  = float(EBak)/nRxn3
        fracSub  = float(SubRxn)/nRxn3
        nRxn4 = nRxn3 - (SubRxn - NSub)
        nRxn5 = nRxn4 - Loops
        fracLoops = float(Loops)/NSub
        nRev  = n3 - nRxn5
        if (nRev < 0):
            print("Error: negative number of reversible reactions!")
        nBlocked = nRxn4 - nRxn0
        check = min([n2, m2 + ThRxn]) - nBlocked
##        if (m3 < check):
##            print('The check fails for model ' + key)
        fracEss = float(numEss)/(nRxn3 - numCorr)
        nonEss = nRxn3 - numCorr - numEss
        fracLet = 2 * float(numLethal)/(nonEss * (nonEss - 1))
        fullStats[key] = [factRed, fracMet1, fracMet2, fracRxn1, fracFor, fracBak, fracSub, fracLoops, fracEss, fracLet, numCorr, minMedia, curStatus]
    return fullStats

def makeAverages(procStatDict, checkLast = True, goodVals = [0,5]): # Computes average statistics; if checkLast, only considers the entries where the last value is one of goodVals
    if checkLast:
        goodKeys = [key for key, value in procStatDict.items() if value[-1] in goodVals]
    else:
        goodKeys = list(procStatDict.keys())
    N = len(goodKeys)
    K = len(procStatDict[goodKeys[0]]) - int(checkLast)
    runningSums = [0] * K
    highest = [''] * K
    lowest  = [''] * K
    highestVal = [0]     * K
    lowestVal  = [10000] * K
    for key in goodKeys:
        cur = procStatDict[key]
        for ind in range(K):
            curElt = cur[ind]
            runningSums[ind] += curElt
            if curElt > highestVal[ind]:
                highestVal[ind] = curElt
                highest[ind] = key
            if curElt < lowestVal[ind]:
                lowestVal[ind]  = curElt
                lowest[ind]  = key
    Averages = [float(Sum)/N for Sum in runningSums]
    return (Averages, highest, lowest, highestVal, lowestVal)

def makeTable(statDict, Titles, filename = 'Table.txt', sep = "\t", skip = [1]):
    f = open(filename, 'w')
    f.write(sep.join([x for i,x in enumerate(Titles) if i not in skip]) + '\n')
    for key in sorted(statDict.keys()):
        value = statDict[key]
        f.write(sep.join([str(x) for i,x in enumerate(value) if i not in skip]) + '\n')
    f.close()
            
        
            
                
                
        
