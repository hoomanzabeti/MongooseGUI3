# This file contains functions for reducing metabolic models (remotely)
# Created by: Leonid Chindelevitch
# Last modified: December 20, 2012

import os, shelve, ModelProcessing, re 
from ClassDefinitions import *
from Utilities import *

numRE = re.compile('[0-9]+')
modelRE = re.compile('[A-Z][A-Z][0-9][c]?')

def findNumEssential(model, biomassIndex):
    for subset in model.reactionSubsets:
        if biomassIndex in zip(*subset.pairs)[0]:
            return len(subset.pairs)
    return 0

def findSmallestIrrev(model, biomassIndex):
    N = model.Matrix
    Irr = [i for i,x in enumerate(model.reactions) if not x.reversible]
    return minimalThermoUnblock(N, Irr, biomassIndex)

def findSmallestMetabs(model, biomassIndex):
    N = model.Matrix
    Irr = [i for i,x in enumerate(model.reactions) if not x.reversible]
    return minimalUnblock(N, Irr, biomassIndex)

def findMedia(model, biomassIndex):
    N = model.Matrix
    Irr = [i for i,x in enumerate(model.reactions) if not x.reversible]
    Ext = [i for i,x in enumerate(model.metabolites) if x.external]
    Exch = classifyExchange(model.fullMatrix, Ext, Irr, extra = True)
    Exch = sum([Exch[x] for x in range(6) if x != 2], []) # omit irreversible export reactions
    return findMinimalMedia(N, biomassIndex, Exch, rec = False, I = Irr)

def processAll(inputShelf = 'ProcessedNetworks', outputShelf = 'ExtraAnalyses', firstTime = False):
    v = shelve.open(inputShelf)
    u = shelve.open(outputShelf)
    if not firstTime:
        for key in sorted(v.keys()):
            if key.endswith('Reduced') and not key.startswith('All'):
                redKey = key[:-7]
                # filename = redKey + 'Reduction.txt'
                filenameL = redKey + 'ReductionFull.txt'
                ListNew = []
                if not filenameL in os.listdir('.'):
                    print(('Processing the ' + redKey + ' model'))
                    ListNew.append(key)
                    cur = v[key]
                    cur.reduceNetwork(redKey + 'Reduction.txt')
                    v[redKey + 'Reduced'] = cur
##                f = open(filenameL, 'r')
##                g = [x.strip().split('\t') for x in f.readlines()]
##                f.close()
##                emptyLines = [x for x, line in enumerate(g) if line == ['']]
##                firstEmptyLine = emptyLines[0]
##                allEmptyLines  = [firstEmptyLine + 2*x for x in range(8)]
##                if not all([(line in emptyLines) for line in allEmptyLines]):
##                    print('Error: this should never happen!')
##                Reduction = [[]] * 8
##                Reduction[0] = [[Fraction(y) for y in x] for x in g[:firstEmptyLine]]
##                for ind in range(1,8):
##                    curIndex = allEmptyLines[ind] - 1
##                    if curIndex not in emptyLines:
##                        Reduction[ind] = [eval(y) for y in g[curIndex]]
##                cur.applyReduction(Reduction)
##                v[redKey + 'Reduced'] = cur
    if not firstTime:
        # D = v['AllReduced']
        D = {}
        for key in sorted(v.keys()):
            shortKey = key[:-7]
            if key.endswith('Reduced') and not key.startswith('All'):
                Found = True
                print(('Processing the ' + shortKey + ' model'))
                cur = v[key]
                m1, n1 = getSize(cur.fullMatrix)
                m2, n2 = getSize(cur.Matrix)
                redMat = cur.reducedMatrix
                # if type(redMat) == type(None):
                # curD = D[shortKey]
                # if curD[0] == curD[2]: # no external metabolites!
                    # Found = False
                    # print("Didn't find the reduced matrix; reducing!")
                    # cur.reduceNetwork(shortKey + 'Reduction.txt')
                    # v[shortKey + 'Reduced'] = cur
                m3, n3 = getSize(cur.reducedMatrix)
                f = open(shortKey + 'Reduction.txt')
                g = f.read()
                f.close()
                numbers = [int(x) for x in re.findall(numRE, g)]
                curBio = cur.findBiomassReaction()
                curStatus, numCorr, minMedia = -1, -1, -1
                if not Found:
                    v[shortKey + 'Reduced'] = cur
                if curBio != -1:
                    curStatus = cur.reactions[curBio].reductionStatus
                    if curStatus in [0, 4, 5]:
                        numCorr = findNumEssential(cur, curBio)
                        curMedia = findMedia(cur, curBio)
                        v[shortKey + 'Media'] = curMedia
                        lenMedia = [len(x) for x in curMedia]
                        if lenMedia:
                            minMedia = min(lenMedia)
                    elif curStatus > 0:
                        S1 = findSmallestMetabs(cur, curBio)
                        u[shortKey + 'Relaxed'] = S1
                        if curStatus == 2:
                            S2 = findSmallestIrrev(cur, curBio)
                            if S2 is not None:
                                u[shortKey + 'RelaxedI'] = S2
                                if len(S2) <= len(S1):
                                    S1 = S2
                                    print("Might as well remove irreversibility constraints")
                                else:
                                    print("It's better to remove metabolite constraints")
                        numCorr = len(S1)
                D[shortKey] = [m1, n1, m2, n2, m3, n3] + numbers + [curStatus, numCorr, minMedia]            
        u['AllReduced'] = D
    v.close()
    u.close()
    return

def findCutsets(inputShelf = 'ProcessedNetworks', outputShelf = 'NewCutsets'):
    u = shelve.open(inputShelf)
    D = u['AllReduced']
    s = shelve.open(outputShelf)
    todo = sorted([x for x,val in D.items() if val[-5] in [0,4]])
    for shortKey in todo:
        key = shortKey + 'Reduced'
        Filename = shortKey + 'Cutsets.txt'
        if key.endswith('Reduced') and shortKey != 'All': # and Filename not in os.listdir('.'):
            print(('Processing the ' + shortKey + ' model'))
            cur = u[key]
            curBio = cur.findBiomassReaction()
            if curBio != -1:
                curStatus = cur.reactions[curBio].reductionStatus
                if curStatus in [0, 4, 5]:
                    redBio = [i for i, subset in enumerate(cur.reactionSubsets) if curBio in zip(*subset.pairs)[0]][0]
                else:
                    print('The biomass reaction is blocked')
                    redBio = -1
                if redBio != -1:
                    print('Processing subsets of size up to 2')
                    curMat = cur.reducedMatrix
                    print(('The dimensions are ' + str(getSize(curMat))))
                    Lens = [len(x.pairs) for x in cur.reactionSubsets]
                    Irr = [ind for ind, x in enumerate(cur.reactionSubsets) if not x.reversible]
                    Essential, Lethal = findEssentialLethal(curMat, redBio, rec = False, I = Irr)
                    f = open(Filename, 'w')
                    f.write(str(Essential) + '\n' + str(Lethal) + '\n')
                    f.close()
                    s[shortKey + 'Essential'] = Essential
                    s[shortKey + 'Lethal']    = Lethal
                    s[shortKey + 'EssentialNum'] = sum([Lens[x] for x in Essential])
                    s[shortKey + 'LethalNum']    = sum([Lens[x[0]] * Lens[x[1]] for x in Lethal])
            else:
                print('A biomass reaction was not found')
    u.close()
    s.close()
    return

def processLast(inputShelf = 'ProcessedNetworks'):
    u = shelve.open(inputShelf)
    D = u['AllReduced']
    todo = sorted([x for x,val in D.items() if val[-5] in [0,4,5] and val[-3] == 0])
    for key in todo:
        print(('Processing the ' + key + ' model'))
        cur = u[key + 'Reduced']
        curBio = cur.findBiomassReaction()
        if curBio != -1:
            curMedia = findMedia(cur, curBio)
            u[key + 'Media'] = curMedia
            lenMedia = [len(x) for x in curMedia]
            if lenMedia:
                minMedia = min(lenMedia)
                D[key][-3] = minMedia
    u['AllReduced'] = D
    u.close()

def processFinal(inputShelf = 'ProcessedNetworks'):
    u = shelve.open(inputShelf)
    D = u['AllReduced']
    todo = sorted([x for x,val in D.items() if val[-3] == 2])
    for shortKey in todo:
        print(('Processing the ' + shortKey + ' model'))
        cur = u[shortKey + 'Reduced']
        curBio = cur.findBiomassReaction()
        curStatus = D[shortKey][-3]
        S1 = findSmallestMetabs(cur, curBio)
        u[shortKey + 'Relaxed'] = S1
        if curStatus == 2:
            S2 = findSmallestIrrev(cur, curBio)
            if S2 is not None:
                u[shortKey + 'RelaxedI'] = S2
                if len(S2) <= len(S1):
                    print("Might as well remove irreversibility constraints")
                    D[shortKey][-2] = len(S2)
                else:
                    print("It's better to remove metabolite constraints")
                    D[shortKey][-2] = len(S1)
    u['AllReduced'] = D
    u.close()

def processEnergy(inputShelf = 'ProcessedNetworks', outputShelf = 'FinalReductions'):
    u = shelve.open(inputShelf)
    # v = shelve.open(outputShelf)
    D = u['AllReduced']
    todo = ['CT1', 'MB1', 'PP2', 'RF1', 'SO1', 'SP1', 'VV1']
    # todo = sorted([x for x in D.keys() if x > 'VV2' or (x > 'EC4' and x <= 'HP2')])
    # todo = sorted([x for x,val in D.iteritems() if val[-5] in [0,4,5]])
    for key in todo:
        # if key not in v:
        print(('Processing the ' + key + ' model'))
        cur = u[key + 'Reduced']
        cur.restoreZeroLoops(shortName = key)
        Irrev = [i for i,r in enumerate(cur.reactionSubsets) if not r.reversible]
        External = cur.findExchangeReactionsReduced()
        # v[key] = fullIterativeReduce(cur.reducedMatrix, Irrev, External, key + 'Reduction.txt')
        fullIterativeReduce(cur.reducedMatrix, Irrev, External, key + 'Reduction.txt')
        # energyBalanceReduce(cur.reducedMatrix, Irrev, External, key + 'NewEnergyReduction.txt')
    u.close()
    # v.close()

def processDistances(inputShelf = 'ProcessedNetworks', outputShelf = 'DistancesLinear'):
    u = shelve.open(inputShelf)
    v = shelve.open(outputShelf)
    D = u['AllReduced']
    todo = sorted([x for x in list(D.keys()) if x > 'PP1'])
    for key in todo:
        print(('Processing the ' + key + ' model'))
        cur = u[key + 'Reduced']
        curBio = cur.findBiomassReaction()
        if curBio != -1:
            curIrrev = [] # [i for i,r in enumerate(cur.reactions) if not r.reversible]
            v[key] = findDistance(cur.Matrix, curBio, curIrrev)
    u.close()
    v.close()

def checkCompleted(inputShelf = 'ProcessedNetworks'):
    u = shelve.open(inputShelf)
    D=u['AllReduced']
    allFiles = os.listdir('.')
    for key in sorted(D.keys()):
        goodFiles = [x for x in allFiles if x.startswith(key + 'ReductionFlux') and not x.endswith('Full.txt')]
        numFiles = [int(x[len(key + 'ReductionFlux'):len(key + 'ReductionFlux')+1]) for x in goodFiles if len(x) > len(key + 'ReductionFlux.txt')]
        if not numFiles:
            print(('No files found for ' + key))
        maxFile = key + 'ReductionFlux' + str(max(numFiles)) + '.txt'
        f = open(maxFile)
        g = f.read()
        f.close()
        numbers = [int(x) for x in re.findall(numRE, g)]
        if any(numbers):
            print((key + ' : ' + str(numbers)))
        maxFileFull = key + 'ReductionEnergy' + str(max(numFiles) + 1) + 'Full.txt'
        try:
            f = open(maxFileFull, 'r')
        except:
            maxFileFull = key + 'ReductionEnergy' + str(max(numFiles)) + 'Full.txt'
        f = open(maxFileFull, 'r')
        g = [x.strip().split('\t') for x in f.readlines()]
        f.close()
        emptyLines = [x for x, line in enumerate(g) if line == ['']]
        firstEmptyLine = emptyLines[0]
        if firstEmptyLine == 0:
            print((key + ' model got reduced to nothing!'))
            continue
        allEmptyLines  = [firstEmptyLine + 2*x for x in range(9)]
        if not all([(line in emptyLines) for line in allEmptyLines]):
            print('Error: this should never happen!')
        allReacts = list(range(len(g[0])))
        try:
            Irrev = [int(x) for x in g[firstEmptyLine + 1]]
        except:
            print((key + ' has no irreversible reactions!'))
            Irrev = []
        try:
            Extern = [int(x) for x in g[firstEmptyLine + 3]]
        except:
            print((key + ' has no external reactions!'))
            Extern = []
        remain = [x for x in allReacts if x not in Irrev + Extern]
        numR = len(remain)
        print((key + ' has ' + str(numR) + ' remaining reactions'))
        if numR < 20:
            signLine = firstEmptyLine + 5
            if signLine not in emptyLines:
                print('The signs have been processed')
            else:
                print('Warning: the signs have not been processed!')
    u.close()

def prepareEFMs(inputShelf = 'ExtraAnalyses'):
    s = shelve.open(inputShelf)
    Di = s['FinalBiomassReactions']
    todo = sorted([key for key,val in Di.items() if val[0] != -1])
    for key in todo:
        print(('Processing the ' + str(key) + ' model'))
        curFile = Di[key][2]
        f = open(curFile, 'r')
        g = [x.strip().split('\t') for x in f.readlines()]
        f.close()
        emptyLines = [x for x, line in enumerate(g) if line == ['']]
        firstEmptyLine = emptyLines[0]
        print(('The matrix has size ' + str(firstEmptyLine - 1) + ' by ' + str(len(g[0]))))
        # Matrix = [[Fraction(y) for y in x] for x in g[:firstEmptyLine]]
        Irrev = [int(x) for x in g[firstEmptyLine + 1]]
        print((str(len(Irrev)) + ' of the reactions are irreversible'))
        # allReacts = range(len(g[0]))
        # allRev = [False if x in Irrev else True for x in allReacts]
        # WriteASCIIMatrix(prepareForCplex(Matrix), key + 'EnergyMatrixNew.txt')
        # WriteASCIIMatrix([allRev], key + 'ReversibilitiesNew.txt')
    s.close()

def prepareCplexAll(inputShelf = 'ProcessedNetworks'):
    u = shelve.open(inputShelf)
    D=u['AllReduced']
    todo = sorted([x for x,val in D.items() if val[-5] not in [0,4,5]])
    for key in todo:
        print(('Processing the ' + str(key) + ' model'))
        cur = u[key + 'Reduced']
        N = cur.Matrix
        Irrev = [i for i,x in enumerate(cur.reactions) if not x.reversible]
        special = cur.findBiomassReaction()
        if special == -1:
            print('No biomass reaction found!')
            continue
        findFeasible(N, special, Irrev, pos = True, Filename = key + 'CPLEX.lp', Cplex = True)
    u.close()

def traceReaction(reactionIndex, cols, Enzymes, energy, External = None):
    # Returns the index of the reaction after reduction, or -1 if it has been deleted
    # Takes as input the column record and the enzyme/isozyme subset record structure
    # Additionally, takes the list of external reactions if energy is set to be True.        
    free = 0
    unidir = 3 if energy else 4
    enzyme = 4 if energy else 5
    if reactionIndex == -1:
        return -1
    if energy:
        survivors = [i for i,x in enumerate(cols) if x in [free, unidir]]
        Internal = [x for x in range(len(survivors)) if x not in External]
        ExternalInit = [survivors[x] for x in External]
        InternalInit = [x for x in range(len(cols)) if x not in ExternalInit]
        if reactionIndex in ExternalInit:
            reducedIndex = External[ExternalInit.index(reactionIndex)]
            return reducedIndex
        else: # reducing everything to only the internal reactions!
            cols = [cols[x] for x in InternalInit]
            reactionIndex = InternalInit.index(reactionIndex)
    status = cols[reactionIndex] 
    if not status in [free, unidir, enzyme]:
        print(('The reaction is blocked due to ' + ('energy' if energy else 'flux') + ' with status ' + str(status)))
        return -1
    anchor_to_subset = {}
    for i, subset in enumerate(Enzymes):
        if energy:
            anchor = subset[0]
            anchor_to_subset[anchor] = subset
        else:
            anchor = subset[0][0]
            anchor_to_subset[anchor] = [x[0] for x in subset]
    anchors = list(anchor_to_subset.keys())
    reducedIndex = -1
    index = 0
    for i, status in enumerate(cols):
        if status in [free, unidir]:
            if reactionIndex == i:
                reducedIndex = index
                break
            index += 1
        if i in anchors:
            subsetReacts = anchor_to_subset[i]
            if reactionIndex in subsetReacts:
                reducedIndex = index
                break
            if not energy:
                index += 1
    if energy:     # mapping the index to its actual value if energy is True
        if reducedIndex == -1:
            print('Error: this should never happen!')
        else:
            reducedIndex = Internal[reducedIndex]
    return reducedIndex

def fullTraceReactions(inputShelf = 'ProcessedNetworks'):
    u = shelve.open(inputShelf)
    D=u['AllReduced']
    todo = sorted([x for x,val in D.items()]) # if val[-5] in [0,4,5]])
    allFiles = os.listdir('.')
    currentIndices = {}
    for key in todo:
        print(('Processing the ' + key + ' model'))
        cur = u[key + 'Reduced']
        cur.restoreZeroLoops(shortName = key)
        curBio = cur.findBiomassReaction()
        if curBio != -1:
            curStatus = cur.reactions[curBio].reductionStatus
            if curStatus in [0, 4, 5]:
                redBio = [i for i, subset in enumerate(cur.reactionSubsets) if curBio in zip(*subset.pairs)[0]][0]
            else:
                print('The biomass reaction is blocked')
                redBio = -1
        base = key + 'Reduction'
        Iter = 0
        reductionFile = ''
        currentIndex = redBio
        fluxFiles = [x for x in allFiles if x.startswith(base + 'Flux') and not x.endswith('Full.txt')]
        fluxNumbers = [int(x[len(base + 'Flux'):len(base + 'Flux')+1]) for x in fluxFiles if len(x) > len(base + 'Flux.txt')]
        energyFiles = [x for x in allFiles if x.startswith(base + 'Energy') and not x.endswith('Full.txt')]
        energyNumbers = [int(x[len(base + 'Energy'):len(base + 'Energy')+1]) for x in energyFiles if len(x) > len(base + 'Energy.txt')]
##        while (True):
##            # print('The current index is ' + str(currentIndex))
##            Iter += 1
##            energy = Iter % 2
##            miniIter = (Iter + 1) / 2
##            extension = 'Energy' if energy else 'Flux'
##            check = (miniIter in energyNumbers) if energy else (miniIter in fluxNumbers)
##            if not check:
##                break
##            reductionFile = base + extension + str(miniIter) + 'Full.txt'
##            # print('Processing file ' + reductionFile)
##            f = open(reductionFile, 'r')
##            g = [x.strip().split('\t') for x in f.readlines()]
##            f.close()
##            emptyLines = [x for x, line in enumerate(g) if line == ['']]
##            firstEmptyLine = emptyLines[0]
##            if firstEmptyLine == 0:
##                print('The model got reduced to nothing!')
##                currentIndex = -1
##                break
##            try:
##                relevantLine = firstEmptyLine + 9 if energy else firstEmptyLine + 5
##                cols = [int(x) for x in g[relevantLine]]
##                External = [int(x) for x in g[firstEmptyLine + 3]] if energy else None
##                Enzymes = [eval(x) for x in g[relevantLine + 6]] if len(g[relevantLine + 6][0]) else []
##            except:
##                print('The file does not have a record of the reactions!')
##                currentIndex = -1
##                break
##            currentIndex = traceReaction(currentIndex, cols, Enzymes, energy, External)
##        # print('The relevant reaction is ' + str(currentIndex))
##        currentIndices[key] = [currentIndex, Iter, reductionFile] 
##    return currentIndices
        totalNumbers = [0] * 7
        for miniIter in energyNumbers:
            reductionFile = base + 'Energy' + str(miniIter) + '.txt'
            f = open(reductionFile, 'r')
            g = f.read()
            f.close()
            curNumbers = [int(x) for x in re.findall(numRE, g)]
            totalNumbers = [curNumbers[i] + y for i,y in enumerate(totalNumbers)]
        currentIndices[key] = totalNumbers
    u.close()
    return currentIndices

def getNumInternal(inputShelf = 'ProcessedNetworks'):
    u = shelve.open(inputShelf)
    D=u['AllReduced']
    todo = sorted([x for x,val in D.items()])
    NumInt = {}
    for key in todo:
        print(('Processing the ' + key + ' model'))
        cur = u[key + 'Reduced']
        cur.restoreZeroLoops(shortName = key)
        Exch = cur.findExchangeReactionsReduced()
        NumInt[key] = len(cur.reactionSubsets) - len(Exch)
    u.close()
    return NumInt

def getFinalSizes(inputShelf = 'ProcessedNetworks'):
    u = shelve.open(inputShelf)
    D=u['AllReduced']
    todo = sorted([x for x,val in D.items()])
    Sizes = {}
    for key in todo:
        print(('Processing the ' + key + ' model'))
        base = key + 'Reduction'
        fluxFiles = [x for x in os.listdir('.') if x.startswith(base + 'Flux') and not x.endswith('Full.txt')]
        fluxNumbers = [int(x[len(base + 'Flux'):len(base + 'Flux')+1]) for x in fluxFiles if len(x) > len(base + 'Flux.txt')]
        energyFiles = [x for x in os.listdir('.') if x.startswith(base + 'Energy') and not x.endswith('Full.txt')]
        energyNumbers = [int(x[len(base + 'Energy'):len(base + 'Energy')+1]) for x in energyFiles if len(x) > len(base + 'Energy.txt')]
        maxF, maxE = max(fluxNumbers), max(energyNumbers)
        (maxN, extension) = (maxF, 'Flux') if maxF >= maxE else (maxE, 'Energy')
        maxFile = base + extension + str(maxN) + 'Full.txt'  
        f = open(maxFile, 'r')
        g = [x.strip().split('\t') for x in f.readlines()]
        f.close()
        firstEmptyLine = g.index([''])
        Sizes[key] = firstEmptyLine, len(g[0])
    u.close()
    return Sizes

def processTest(inputShelf = 'TestStuff'):
    u = shelve.open(inputShelf)
    Mat = u['TestMatrixExt']
    Irr = u['TestIrrExt']
    External = u['TestExternalExt']
    Result = reduceMatrix(Mat, Irr, 'TestReductionExt.txt')
    u['TestOutputExt'] = Result
    New = Result[0]
    NewIrr = Result[1]
    NewExt = findExternal(External, Result[2:])
    Final = fullIterativeReduce(New, NewIrr, NewExt, Filename = 'TestIterativeReductionExt.txt')
    u['TestOutputFullReduceExt'] = Final
    u.close()
    return

##os.chdir('../')
##processTest()
##os.chdir('../GSMNs')
##findCutsets()
##Di = fullTraceReactions()
##prepareEFMs()
##Nums = getNumInternal()
##Sizes = getFinalSizes()

##s = shelve.open('YeastConsensus')
##model = s['YM']
##biomassIndex = model.findBiomassReaction()
##redBio = [i for i, subset in enumerate(model.reactionSubsets) if biomassIndex in zip(*subset.pairs)[0]][0]
##curMat = model.reducedMatrix
##Irr = [ind for ind, x in enumerate(model.reactionSubsets) if not x.reversible]
##s['YMEssential'], s['YMLethal'] = findEssentialLethal(curMat, redBio, rec = False, I = Irr)
##s.close()

##os.chdir('../GSMNs')
##curList = os.listdir('.')
##allProc = [x for x in curList if x.endswith('QS.lp')]
##print('There are a total of ' + str(len(allProc)) + ' files to process')
##for ind, filename in enumerate(allProc):
##    if ind % 10 == 0:
##        print('Processed ' + str(ind) + ' so far')
##    val = processFile(filename, destroyIn = False)
##    if not (type(val0) == type([]) and len(val0) == 0):
##        print('Error on ' + filename + ': this should never happen!')
