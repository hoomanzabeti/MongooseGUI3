# This file contains functions in the user interface of MetaMerge
# Created by: Leonid Chindelevitch
# Last modified: January 13, 2012

import shelve
from MatchProcessing import *
from ReactionMatching import *
from GeneProcessing import findOverlap
from OutputProcessing import CreateSMatrix

def MatchingMaster(Reacts0, Reacts1, Metabs0, Metabs1, ReactFeats0, ReactFeats1, MetabFeats0, MetabFeats1, ReactMatch, MetabMatch, external0, external1, cutoff = 3, startData = []):
    # This function proposes a merger of two networks and records the history of the user's decisions
    # It returns the final match matrices which can then be post-processed to get the merged network
    # Also returns a record of the user's decisions; a list of pairs of oppositely directed reactions
    # Inputs: the sets of reactions (in list representation), the sets of metabolites, the features
    # of both reactions and metabolites (in table form), the matrix of reaction matching scores and
    # metabolite matching scores, the indices of metabolites having external versions of themselves,
    # and a cutoff value for the initial matching to be proposed (the default value we use is 30).
    # If an earlier run had started but did not complete, it can be loaded from the startData shelf!
    # NOTE: Side effects might include a modification of the ReactMatch and the MetabMatch matrices!
    # NOTE: There is an assumption that each reaction is sorted in increasing order of coefficients!
    Iter = 0

    # stores length of parameters into vars
    m0, n0 = len(Metabs0), len(Reacts0)
    m1, n1 = len(Metabs1), len(Reacts1)

    #if the length doesnt match matrix size, output error
    if len(ReactMatch) != n0 or len(ReactMatch[0]) != n1:
        print(('Incorrect matrix size! It should be ' + str(n0) + ' by ' + str(n1)))
        return
    if len(MetabMatch) != m0 or len(MetabMatch[0]) != m1:
        print(('Incorrect matrix size! It should be ' + str(m0) + ' by ' + str(m1)))
    ReactBin, MetabBin = [[0 for y in range(n1)] for x in range(n0)], [[0 for y in range(m1)] for x in range(m0)]
    Matches = sorted([(x,y, ReactMatch[x][y]) for x in range(n0) for y in range(n1) if ReactMatch[x][y] >= cutoff], key = lambda z:z[2], reverse=True)
    N = len(Matches)
    index = 0
    cmd = 'X'
    ReactHistory, MetabHistory = [], []
    Opposite, Deferred = [], []
    Extras0, Extras1 = [], []
    Codes = {0: 'not decided', 1: 'accepted', -1: 'rejected'}
    print("Welcome to the Matching Master, a smart program that assists you with the matching of metabolic networks!")
    print("The 1st network we are looking at today contains " + str(m0) + " metabolites and " + str(n0) + " reactions.")
    print("The 2nd network we are looking at today contains " + str(m1) + " metabolites and " + str(n1) + " reactions.\n\n")
    Unused = list(range(N))
    remain = N
    if startData:
        # substitute everything with the values that have been saved!
        s = shelve.open(startData)
        ReactBin, MetabBin = s['ReactBin'], s['MetabBin']
        Opposite, Deferred = s['Opposite'], s['Deferred']
        ReactHistory, MetabHistory = s['ReactHistory'], s['MetabHistory']
        Extras0, Extras1 = s['Extras0'], s['Extras1']
        if 'Matches' in list(s.keys()):
            Matches = s['Matches']
            N = len(Matches)
            remain = N
            Unused = list(range(N))
        else:
            processHistory(MetabMatch, MetabHistory)
            remain = 0
    while(cmd != 'Q'):
        if (remain):
            print(("Displaying the " + str(Unused[index]) + end(Unused[index]) + " matching pair of reactions out of " + str(N)))
            # Display the current pair of matching reactions
            (react0,react1,z) = Matches[Unused[index]]
            print(('(' + str(react0) + ') ' + ReconstructReaction(Reacts0[react0], Metabs0, 0)))
            print(('(' + str(react1) + ') ' + ReconstructReaction(Reacts1[react1], Metabs1, 1)))
            print(('The matching score for this pair of reactions is ' + str(z)))
            print(('Your previous decision for this reaction pair was: ' + Codes[ReactBin[react0][react1]]))
            cmd = input("Please enter A to accept, R to reject, or H for help: ")
            if cmd == 'Q':
                continue
            elif cmd == 'H':
                helpme()
                continue
            elif cmd == 'P':
                index -= 1
                if index == -1:
                    print("You are currently at the beginning; transfering you to the end")
                    index += remain
                continue
            elif cmd == 'F':
                index += 1
                if index == remain:
                    print("You are currently at the end; transfering you to the beginning")
                    index -= remain
                continue
            elif cmd == 'R':
                print("You have chosen to reject the current pair of reactions\n\n")
                ReactHistory.append([react0,react1,-1])
                ReactBin[react0][react1] = -1
                Unused.pop(index)
                remain -= 1
                if index == remain:
                    print("You are currently at the end; transfering you to the beginning")
                    index -= remain
                continue
            elif cmd == 'I':
                displayFeatures(react0, react1, ReactFeats0, ReactFeats1, 'react')
                continue
            elif cmd == 'M':
                displayMatches (react0, react1, (Reacts0, Metabs0), (Reacts1, Metabs1), ReactBin, 'react')
                continue
            elif cmd == 'A':
                print("You have chosen to accept the current pair of reactions\n")
                ReactHistory.append([react0,react1,1])
                ReactBin[react0][react1] = 1
                Unused.pop(index)
                remain -= 1
                if index == remain:
                    print("You are currently at the end; transfering you to the beginning")
                    index -= remain
                # Display the proposed matching of metabolites
                react0, react1 = Reacts0[react0], Reacts1[react1]
                Mets0P, Mets0N, Mets1P, Mets1N = [t[0] for t in react0 if t[1] > 0],[t[0] for t in react0 if t[1] < 0],[t[0] for t in react1 if t[1] > 0],[t[0] for t in react1 if t[1] < 0]
                L0P, L0N, L1P, L1N = len(Mets0P), len(Mets0N), len(Mets1P), len(Mets1N)
                (mMets0PP, mMets1PP, uMets0PP, uMets1PP, scorePP) = maxMatching(MetabMatch, Mets0P, Mets1P)
                (mMets0NN, mMets1NN, uMets0NN, uMets1NN, scoreNN) = maxMatching(MetabMatch, Mets0N, Mets1N)
                (mMets0NP, mMets1NP, uMets0NP, uMets1NP, scoreNP) = maxMatching(MetabMatch, Mets0N, Mets1P)
                (mMets0PN, mMets1PN, uMets0PN, uMets1PN, scorePN) = maxMatching(MetabMatch, Mets0P, Mets1N)
                if (scorePP + scoreNN) >= (scorePN + scoreNP):
                    mMets0, mMets1 = (mMets0NN + mMets0PP), (mMets1NN + mMets1PP)
                    uMets0, uMets1 = (uMets0NN + uMets0PP), (uMets1NN + uMets1PP)
                    Mets0, Mets1 = (Mets0N + Mets0P), (Mets1N + Mets1P)
                else:
                    Opposite.append([Reacts0.index(react0),Reacts1.index(react1)])
                    # reverse the correspondence between the reactions!
                    react1 = sorted([[x[0],-x[1]] for x in react1], key = lambda x:x[1])
                    mMets0, mMets1 = (mMets0NP + mMets0PN), (mMets1NP + mMets1PN)
                    uMets0, uMets1 = (uMets0NP + uMets0PN), (uMets1NP + uMets1PN)
                    Mets0, Mets1 = (Mets0N + Mets0P), (Mets1P + Mets1N)
                L, L0, L1 = len(mMets0), len(Mets0), len(Mets1)
                mMetabs0, uMetabs0 = [react0[Mets0.index(x)] for x in mMets0], [react0[Mets0.index(x)] for x in uMets0]
                mMetabs1, uMetabs1 = [react1[Mets1.index(x)] for x in mMets1], [react1[Mets1.index(x)] for x in uMets1]
                allMetabs0, allMetabs1 = mMetabs0 + uMetabs0, mMetabs1 + uMetabs1
                print('The reactions with the proposed matching look like this (unmatched metabolites, if any, are highlighted):')
                print((ReconstructReaction(allMetabs0, Metabs0, 0, list(range(L, L0)), True)))
                print((ReconstructReaction(allMetabs1, Metabs1, 1, list(range(L, L1)), True)))
                dec = 'X'
                while(dec != 'Y' and dec != 'C'):
                    dec = input('Please enter Y to accept all pairs of metabolites with matching numbers, or C to continue: ')
                    if dec == 'Y':
                        print('You have decided to match all pairs of metabolites with matching numbers!')
                        for pos in range(L):
                            u, v = mMets0[pos], mMets1[pos]
                            MetabHistory.append([u,v,1])
                            MetabBin[u][v] = 1
                            if MetabMatch[u][v] < 100: # give a large bonus to the newlymatcheds!
                                MetabMatch[u][v] += 100
                            if u in external0 and v in external1: # automatically defer the decision on external pairs!
                                uExt, vExt = m0 - len(external0) + external0.index(u), m1 - len(external1) + external1.index(v)
                                if [uExt, vExt] not in Deferred:
                                    Deferred.append([uExt, vExt])
                        pos = L
                        break
                    elif dec == 'C':
                        pos = 0
                        break
                    else:
                        print(('You have entered an unrecognized command ' + dec))
                        continue
                while(pos < L):
                    u, v = mMets0[pos], mMets1[pos]
                    print((ReconstructReaction(allMetabs0, Metabs0, 0, [pos], True)))
                    print((ReconstructReaction(allMetabs1, Metabs1, 1, [pos], True)))
                    print(('The matching score for this pair of metabolites is ' + str(MetabMatch[u][v])))
                    print(('Your previous decision for this metabolite pair was: ' + Codes[MetabBin[u][v]]))
                    cmd = input("Please enter A to accept, R to reject, or H for help: ")
                    if cmd == 'H':
                        helpme()
                        continue
                    elif cmd == 'I':
                        displayFeatures(u, v, MetabFeats0, MetabFeats1, 'metab')
                        continue
                    elif cmd == 'M':
                        displayMatches (u, v, Metabs0, Metabs1, MetabBin, 'metab')
                        continue
                    elif cmd == 'R':
                        print("You have chosen to reject the current pair of metabolites\n")
                        uMets0.append(u)
                        uMets1.append(v)
                        MetabHistory.append([u,v,-1])
                        MetabBin[u][v] = -1
                        MetabMatch[u][v] = 0 # they will never again be matched automatically
                        pos += 1
                        continue
                    elif cmd == 'A':
                        print("You have chosen to accept the current pair of metabolites\n")
                        MetabHistory.append([u,v,1])
                        MetabBin[u][v] = 1
                        if MetabMatch[u][v] < 100: # give a large bonus to the newlymatcheds!
                            MetabMatch[u][v] += 100
                        pos += 1
                        # see if there are corresponding external metabolites to be matched
                        if u in external0 and v in external1:
                            uExt, vExt = m0 - len(external0) + external0.index(u), m1 - len(external1) + external1.index(v)
                            if MetabBin[uExt][vExt] == 0:
                                ext0, ext1 = Metabs0[uExt], Metabs1[vExt]
                                print(("Would you also like to match the external metabolites " + ext0 + " and " + ext1 + " ?"))
                                print(('The matching score for this pair of metabolites is ' + str(MetabMatch[uExt][vExt])))
                                print(('Your previous decision for this metabolite pair was: ' + Codes[MetabBin[uExt][vExt]]))
                                dec = 'X'
                                while (dec != 'A' and dec != 'R' and dec != 'D'):
                                    dec = input("Please enter A to accept, R to reject, or D to defer the decision: ")
                                    if dec == 'A':
                                        print("You have chosen to accept the current pair of metabolites\n")
                                        MetabHistory.append([uExt,vExt,1])
                                        MetabBin[uExt][vExt] = 1
                                        if MetabMatch[uExt][vExt] < 100: # give a large bonus to the newlymatcheds!
                                            MetabMatch[uExt][vExt] += 100
                                        break
                                    elif dec == 'R':
                                        print("You have chosen to reject the current pair of metabolites\n")
                                        MetabHistory.append([uExt,vExt,-1])
                                        MetabBin[uExt][vExt] = -1
                                        MetabMatch[uExt][vExt] = 0 # they will never again be matched automatically
                                        break
                                    elif dec == 'D':
                                        if [uExt, vExt] not in Deferred:
                                            Deferred.append([uExt, vExt])
                                        break
                                    else:
                                        continue
                        else:
                            continue
                    else:
                        print(("You have entered an unrecognized command: " + cmd))
                    continue
                else: # this means pos = L
                    unmatched0, unmatched1 = [Metabs0[i] for i in uMets0], [Metabs1[i] for i in uMets1]
                    prevNull0, prevNull1 = [i for i in uMets0 if i in Extras0], [i for i in uMets1 if i in Extras1]
                    justMetabs0, justMetabs1 = [i[0] for i in allMetabs0], [i[0] for i in allMetabs1]
                    unmatchedInds0, unmatchedInds1 = [justMetabs0.index(j) for j in uMets0], [justMetabs1.index(j) for j in uMets1]
                    prevNullInds0, prevNullInds1 = [justMetabs0.index(j) for j in prevNull0], [justMetabs1.index(j) for j in prevNull1]
                    if (unmatched0 and unmatched1):
                        print('Please note that you will need to match the remaining metabolites yourself')
                        if (prevNullInds0 or prevNullInds1):
                            print('The reactions with the metabolites previously matched to nothing highlighted look like this:')
                            print((ReconstructReaction(allMetabs0, Metabs0, 0, prevNullInds0, True)))
                            print((ReconstructReaction(allMetabs1, Metabs1, 1, prevNullInds1, True)))
                            dec = 'X'
                            while (dec != 'Y' and dec != 'C'):
                                dec = input("Please enter Y to match all these metabolites to nothing, or C to continue: ")
                                if dec == 'Y':
                                    print('You have decided to match all these metabolites to nothing!')
                                    for u in prevNull0:
                                        MetabHistory.append([u,-1,1])
                                        # update the data structures
                                        unmatched0.remove(Metabs0[u])
                                        unmatchedInds0.pop(uMets0.index(u))
                                        uMets0.remove(u)
                                    for v in prevNull1:
                                        MetabHistory.append([-1,v,1])
                                        # update the data structures
                                        unmatched1.remove(Metabs1[v])
                                        unmatchedInds1.pop(uMets1.index(v))
                                        uMets1.remove(v)
                                    break
                                elif dec == 'C':
                                    break
                                else:
                                    print(('You have entered an unrecognized command ' + dec))
                                    continue
                    while(unmatched0 and unmatched1):
                        print('The reactions with the unmatched metabolites highlighted look like this:')
                        print((ReconstructReaction(allMetabs0, Metabs0, 0, unmatchedInds0, True)))
                        print((ReconstructReaction(allMetabs1, Metabs1, 1, unmatchedInds1, True)))
                        cmd1 = input("Please use N for numeric matching, T for text matching or H for help: ")
                        if cmd1 == 'H':
                            helpme()
                            continue
                        elif cmd1[0] == 'N':
                            args = cmd1.split()
                            if len(args) != 3:
                                print("Sorry, your command should have the form N m1 m2 (where m1 and m2 are integers).")
                                continue
                            else:
                                k0, k1 = args[1], args[2]
                                try:
                                    n0, n1 = int(k0), int(k1)
                                    if n0 not in unmatchedInds0 and n0 != -1:
                                        print(("Sorry, the metabolite's index should be -1 or an element of " + ', '.join([str(x) for x in unmatchedInds0])))
                                        continue
                                    elif n1 not in unmatchedInds1 and n1 != -1:
                                        print(("Sorry, the metabolite's index should be -1 or an element of " + ', '.join([str(x) for x in unmatchedInds1])))
                                    elif n0 == -1 and n1 == -1:
                                        print("Sorry, matching nothing to nothing is not allowed")
                                    else:
                                        if n1 == -1:
                                            u = justMetabs0[n0]
                                            if u not in Extras0:
                                                Extras0.append(u)
                                            print(("You have chosen to match the metabolite " + Metabs0[u] + " to nothing\n"))
                                            MetabHistory.append([u,-1,1])
                                            # update the data structures
                                            unmatched0.remove(Metabs0[u])
                                            unmatchedInds0.pop(uMets0.index(u))
                                            uMets0.remove(u)
                                            continue
                                        elif n0 == -1:
                                            v = justMetabs1[n1]
                                            if v not in Extras1:
                                                Extras1.append(v)
                                            print(("You have chosen to match the metabolite " + Metabs1[v] + " to nothing\n"))
                                            MetabHistory.append([-1,v,1])
                                            # update the data structures
                                            unmatched1.remove(Metabs1[v])
                                            unmatchedInds1.pop(uMets1.index(v))
                                            uMets1.remove(v)
                                            continue
                                        else:
                                            u, v = justMetabs0[n0], justMetabs1[n1]
                                            print(("You have chosen to match the metabolites " + Metabs0[u] + " and " + Metabs1[v] + "\n"))
                                            MetabHistory.append([u,v,1])
                                            MetabBin[u][v] = 1
                                            if MetabMatch[u][v] < 100: # give a large bonus to the newlymatcheds!
                                                MetabMatch[u][v] += 100
                                            # update the data structures
                                            unmatched0.remove(Metabs0[u])
                                            unmatched1.remove(Metabs1[v])
                                            unmatchedInds0.pop(uMets0.index(u))
                                            unmatchedInds1.pop(uMets1.index(v))
                                            uMets0.remove(u)
                                            uMets1.remove(v)
                                            continue
                                except:
                                    print(("One of the numbers " + k0 + ", " + k1 + " you have entered is not a valid numeric index."))
                                    continue
                        elif cmd1[0] == 'T':
                            args = cmd1.split()
                            if len(args) != 3:
                                print("Sorry, your command should have the form T m1 m2 (where m1 and m2 are metabolites).")
                                continue
                            else:
                                k0, k1 = args[1], args[2]
                                if k0 not in unmatched0:
                                    print(("Sorry, the metabolite " + k0 + " is not in the unmatched set " + ', '.join(unmatched0)))
                                    continue
                                elif k1 not in unmatched1:
                                    print(("Sorry, the metabolite " + k1 + " is not in the unmatched set " + ', '.join(unmatched1)))
                                    continue
                                else:
                                    u, v = Metabs0.index(k0), Metabs1.index(k1)
                                    print(("You have chosen to match the metabolites " + k0 + " and " + k1 + "\n"))
                                    MetabHistory.append([u,v,1])
                                    MetabBin[u][v] = 1
                                    if MetabMatch[u][v] < 100: # give a large bonus to the newlymatcheds!
                                        MetabMatch[u][v] += 100
                                    # update the data structures
                                    unmatched0.remove(k0)
                                    unmatched1.remove(k1)
                                    unmatchedInds0.pop(uMets0.index(u))
                                    unmatchedInds1.pop(uMets1.index(v))
                                    uMets0.remove(u)
                                    uMets1.remove(v)
                                    continue
                        else:
                            print(("You have entered an unrecognized command: " + cmd1))
                            continue
                    else: # we have completed the matching
                        if unmatched0:
                            print(("The metabolites remaining in the 1st reaction: " + ', '.join(unmatched0) + " will be matched to nothing"))
                            for t in unmatched0:
                                u = Metabs0.index(t)
                                MetabHistory.append([u,-1,1])
                                if u not in Extras0:
                                    Extras0.append(u)
                        elif unmatched1:
                            print(("The metabolites remaining in the 2nd reaction: " + ', '.join(unmatched1) + " will be matched to nothing"))
                            for t in unmatched1:
                                v = Metabs1.index(t)
                                MetabHistory.append([-1,v,1])
                                if v not in Extras1:
                                    Extras1.append(v)
                        print("The matching for the current pair of reactions is complete!\n\n")
                        continue
            else:
                print(("You have entered an unrecognized command: " + cmd))
                continue
        else: # need to generate new matches
            # start by saving the results so far!
            s = shelve.open('MatchingMaster')
            s['ReactBin'], s['MetabBin'] = ReactBin, MetabBin
            s['ReactHistory'], s['MetabHistory'] = ReactHistory, MetabHistory
            s['Opposite'], s['Deferred'] = Opposite, Deferred
            s['Extras0'], s['Extras1'] = Extras0, Extras1
            s.close()
            print("Please wait while I generate some new matches for you!")
            (perfectMatch, consider) = expandMatch(Reacts0, Reacts1, Metabs0, Metabs1, ReactBin, MetabBin, 'NewMatches' + str(Iter) + '.txt')
            Iter += 1
            if perfectMatch:
                Matches = sorted([(x[0],x[1],ReactMatch[x[0]][x[1]]) for x in perfectMatch], key = lambda z:z[2], reverse=True)
                N = len(Matches)
                Unused = list(range(N))
                remain = N
                index = 0
                continue
            elif consider:
                # have to go with imperfect matches...
                Matches = sorted([(x[0],x[1],ReactMatch[x[0]][x[1]]) for x in consider], key = lambda z:z[2], reverse=True)
                N = len(Matches)
                Unused = list(range(N))
                remain = N
                index = 0
                continue
            else:
                # no matches are possible...
                print("It looks like all the good matches have been exhausted!")
                if (Deferred):
                    print("Please make a decision on the following external metabolites")
                    for pair in Deferred:
                        uExt, vExt = pair[0], pair[1]
                        ext0, ext1 = Metabs0[uExt], Metabs1[vExt]
                        print(("Would you also like to match the external metabolites " + ext0 + " and " + ext1 + " ?"))
                        print(('The matching score for this pair of metabolites is ' + str(MetabMatch[uExt][vExt])))
                        print(('Your previous decision for this metabolite pair was: ' + Codes[MetabBin[uExt][vExt]]))
                        dec = 'X'
                        while (dec != 'A' and dec != 'R'):
                            dec = input("Please enter A to accept or R to reject: ")
                            if dec == 'A':
                                print("You have chosen to accept the current pair of metabolites\n")
                                MetabHistory.append([uExt,vExt,1])
                                MetabBin[uExt][vExt] = 1
                                if MetabMatch[uExt][vExt] < 100: # give a large bonus to the newlymatcheds!
                                    MetabMatch[uExt][vExt] += 100
                                break
                            elif dec == 'R':
                                print("You have chosen to reject the current pair of metabolites\n")
                                MetabHistory.append([uExt,vExt,-1])
                                MetabBin[uExt][vExt] = -1
                                MetabMatch[uExt][vExt] = 0 # they will never again be matched automatically
                                break
                            elif dec == 'I':
                                displayFeatures(uExt, vExt, MetabFeats0, MetabFeats1, 'metab')
                                continue
                            elif dec == 'M':
                                displayMatches(uExt, vExt, Metabs0, Metabs1, MetabBin, 'metab')
                                continue
                            else:
                                print(("You have entered an unrecognized command: " + dec))
                                continue
                    Deferred = []
                    continue
                else:
                    break
    print("Saving the information and closing the session. Thank you for using the Matching Master!")
    return (ReactBin, MetabBin, ReactHistory, MetabHistory, Opposite, Extras0, Extras1)

def helpme():
    # The help function for the Matching Master
    print("The following is a list of commands you can use in the matching generation mode (the default mode):")
    print("Entering P lets you move to the Previous matching reaction pair; F, to the Following matching reaction pair")
    print("When viewing a pair of matching reactions, you may enter A to Accept the proposed pair or R to Reject it.")
    print("You can also request more information about the currently proposed pair of reactions by using the command I.")
    print("You can also see any reactions that are Matched to the currently proposed match by using the command M.\n")
    print("After you accept a pair, you will be asked to specify a course of action for the unmatched metabolites.")
    print("The metabolites will be presented to you in a way that matches them well according to the program's knowledge.")
    print("The metabolites that matched will be displayed with matching numbers and the unmatched ones will be highlighted.")
    print("You have the option to accept all the matched metabolites immediately by entering Y, or continue by entering C.")
    print("Possible metabolite matches will then be highlighted in their respective reactions, e.g. << atp >> << ATP >>.")
    print("When viewing a pair of matching metabolites, you may enter A to Accept the proposed pair or R to Reject it.")
    print("You can also request more Information about the currently highlighted metabolites by using the command I.")
    print("You can also see any metabolites that are Matched to the currently highlighted ones by using the command M.\n")
    print("If there are unmatched metabolites in an accepted pair of matching reactions, you will be asked to match them.")
    print("You have the option Y to match to nothing those that were previously matched to nothing, or press C to continue.")
    print("The command N is for Numerical match; it lets you match a pair of metabolites specified by their numbers.")
    print("e.g. N 5 4 will match the 5-th metabolite in the first reaction to the 4-th metabolite in the second one.")
    print("You can also match a metabolite to nothing by using -1, the largest index not used in the numbering scheme.")
    print("e.g. N 5 -1 will mean that the 5-th metabolite in the first reaction is not present in the second reaction.\n")
    print("The command T is for Textual match; it lets you add a pair of metabolites to the matching in text form.")
    print("For example, T atp ATP will match the metabolite atp in the first reaction to ATP in the second one.")
    print("Lastly, the command Q (Quit program) allows you to save what you have done so far and quit the program.")
    print("To see this message again at any time, type H (Help) followed, just like for any other command, by Enter.\n")
    return

def displayFeatures(Element0, Element1, Features0, Features1, option = 'metab'):
    # This function compares the features of two given elements; the valid options are "metab" and "react"
    # Supply the indices of the two elements as well as the corresponding lists of all the features
    Titles = getTitles(option)
    L = len(Titles)
    rel0, rel1 = Features0[Element0], Features1[Element1]
    # print the feature name with the value if it's known, N/A otherwise
    print(('\n'.join([Titles[i]+': '+str(rel0[i])*int(bool(rel0[i]))+'N/A'*(1-int(bool(rel0[i])))+' vs. '+str(rel1[i])*int(bool(rel1[i]))+'N/A'*(1-int(bool(rel1[i]))) for i in range(L)])))
    return

def getTitles(option):
    # This function returns the titles of the features corresponding to the option (either 'metab' or 'react')
    if option == 'react':
        Titles = ['Gene names', 'Protein names', 'Reaction name', 'Pathway/Subsystem', 'Enzyme names/Reaction name']
    elif option == 'metab':
        Titles = ['name', 'Species name', 'IUPAC name', 'CAS number', 'Formula', 'Biocyc', 'kegg-id']
    else:
        print(('Error: unrecognized option! ' + option))
        Titles = []
    return Titles

def displayMatches(Element0, Element1, List0, List1, Matches, option = 'metab'):
    # This function displays the matches of two given elements; the valid options are "metab" and "react"
    # For metabolites, supply the indices of the two metabolites as well as the full lists of metabolites
    # For reactions, supply the indices and the full lists of reactions and metabolites for each network!
    if option == 'metab':
        Item0, Match0 = List0[Element0], [(x,List1[x]) for x in range(len(List1)) if Matches[Element0][x] > 0]
        Item1, Match1 = List1[Element1], [(x,List0[x]) for x in range(len(List0)) if Matches[x][Element1] > 0]
    elif option == 'react':
        Reacts0, Metabs0, Reacts1, Metabs1 = List0[0], List0[1], List1[0], List1[1]
        item0, match0 = Reacts0[Element0], [x for x in range(len(Reacts1)) if Matches[Element0][x] > 0]
        item1, match1 = Reacts1[Element1], [x for x in range(len(Reacts0)) if Matches[x][Element1] > 0]
        Item0, Match0 = ReconstructReaction(item0, Metabs0, 0), [(x, ReconstructReaction(Reacts1[x], Metabs1, 1)) for x in match0]
        Item1, Match1 = ReconstructReaction(item1, Metabs1, 1), [(x, ReconstructReaction(Reacts0[x], Metabs0, 0)) for x in match1]
    else:
        print(('Error: unrecognized option! ' + option))
        return
    if Match0:
        print(('(' + str(Element0) + ') ' + Item0 + ' has already been matched to ' + ', '.join(['(' + str(x[0]) + ') ' + x[1] for x in Match0])))
    else:
        print(('(' + str(Element0) + ') ' + Item0 + ' has not been matched to anything yet!'))
    if Match1:
        print(('(' + str(Element1) + ') ' + Item1 + ' has already been matched to ' + ', '.join(['(' + str(x[0]) + ') ' + x[1] for x in Match1])))
    else:
        print(('(' + str(Element1) + ') ' + Item1 + ' has not been matched to anything yet!'))
    if Match0 and Match1:
        Between = [str(x[0]) + ' : ' + str(y[0]) for x in Match1 for y in Match0 if Matches[x[0]][y[0]] > 0]
        if Between:
            print(('Out of these, ' + ', '.join(Between) + ' have also been matched to each other'))
        else:
            print('However, none of these have been matched to each other')
    return

def processHistory(Matches, History, bonus = 10):
    # This function modifies the Matches matrix in place(!) based on a recorded sequence of decisions
    for item in History:
        (x,y,z) = item
        if z == -1:
            Matches[x][y] = 0
        elif z == 1 and Matches[x][y] < bonus:
            Matches[x][y] += bonus
    return
