def strainOptimize(N, target, biomass, Lambda = 1,weight = [1]):
    # Finds a small knockout that results in a production of at least Lambda units of
    # flux through the target reaction for each unit of flux through the biomass one.
    # If the desired target flux level cannot be reached, the function returns False.
    # NOTE: For future processing, may need to implement a cutset minimality check (!)
    extraB = {target: (-1.0/Lambda)}
    (val, vec) = findMinAdded(N, biomass, Filename = 'OptStrain.lp', extra = extraB, weight = weight)
    return [int(x[1:]) for x in vec if x.startswith('T')]

def augmentWeight(w,listOfReactions,factor = 2):
    #increase the weight in w of all the reactions in <listOfReactions>
    #by <factor>.  This is used to punish reactions of non-target-minimal
    #cutsets 
    for rx in listOfReactions:
        w[rx] *= 2
    return w

def iterativeStrainOptimize(N,target,biomass, numTrials = 100):
    w = [1]*(len(N[0]))# should probably be len(N[0])
    #initialize

    for i in range(numTrials):
        listOfReactions = strainOptimize(N,target,biomass,1,w)
        #guaranteed to be a cut set
        try:
            listOfReactions.remove(target)
        except ValueError:
            #target is not in list of reactions
            pass

        if testCutSet(listOfReactions,N,biomass):
            augmentWeight(w,listOfReactions)
        else:
            print(max(w))
            return listOfReactions

    print('Failure')
    print(listOfReactions)
    print(max(w))
    
def findMinCutSetsByShuffling(cutset, N,biomass,numiter = 10):
    """Starting with a cutset"""
    l = []
    for i in range(numiter):
        random.shuffle(cutset)
        l.append(sorted(extractMinimal(cutset, testCutSet, (N,biomass))))
    return l
