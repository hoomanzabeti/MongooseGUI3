# This file contains functions for preparing the features for the network merging process
# Created by: Leonid Chindelevitch
# Last modified: November 30, 2016

import time, re
from urllib.request import urlopen
from Utilities import *
from functools import reduce
funcWords = ['al.', 'and', 'are', 'but', 'can', 'et.', 'for', 'has', 'is', 'not', 'of', 'or', 'out', 'see', 'the', 'to', 'via']
CASRE = re.compile("[0-9]*-[0-9]*-[0-9]")
ProteinRE = re.compile('[A-Z]{1}[a-z]{2,3}[A-Z]{0,3}[0-9]{0,2}')
atomRE = re.compile('[A-Z]{1}[a-z]*[0-9]*')

def flattenFeature(reactions, featureName):
    # For each specified reaction, transforms a list corresponding to featureName into a string 
    for reaction in reactions:
        if reaction.description and featureName in reaction.description:
            reaction.description[featureName] = ' '.join(reaction.description[featureName])
    return

def collectSpeciesNames(species):
    # For each specified species, adds the species name into the description of the species
    for spec in species:
        if spec.description:
            spec.description['Species name'] = spec.name
    return

def collectReactionNames(reactions):
    # For each specified reaction, adds the reaction name into the description of the reaction
    for reaction in reactions:
        if reaction.description:
            reaction.description['Reaction name'] = reaction.name
    return

def fixGeneCombinations(reactions):
    # For each specified reaction, changes a missing gene combination into an empty one
    for reaction in reactions:
        if reaction.geneCombination is None:
            reaction.geneCombination = CNF([[]])
    return

def collectAllGeneNames(reactions):
    # For each specified reaction, collects the gene names into a list
    for reaction in reactions:
        if reaction.geneCombination and reaction.geneCombination.clauses and reaction.description:
            reaction.description['Gene names'] = sorted(list(reduce(lambda x,y:set(x).union(y), reaction.geneCombination.clauses, [])))
    return

def splitAllProteinNames(reactions, andSymbol = '+', orSymbol = ','):
    # For each specified reaction, splits the protein combination into individual protein names
    for reaction in reactions:
        if reaction.description and 'Protein' in reaction.description:
            curProteins = reaction.description['Protein']
            splitProteins = [x.strip() for x in curProteins.split(orSymbol)]
            splitProteins = sum([[x.strip() for x in protein.split(andSymbol)] for protein in splitProteins], [])
            reaction.description['Protein names'] = splitProteins
    return

def findAllProteinNames(reactions):
    # For each specified reaction, finds the likely protein names
    for reaction in reactions:
         if reaction.name:
             if reaction.description:
                 reaction.description['Protein names'] = findProteinNames(reaction.name)
    return

def findProteinNames(Line):
    # Returns the substrings that are likely to be protein names
    if not Line:
        return []
    else:
        line = Line.replace(',','')
        line = line.replace('(','')
        line = line.replace(')','')
        line = line.replace('/',' ')
        line = line.split(' ')
        line = [x.title() for x in line]
        found = sorted([x for x in line if re.match(ProteinRE, x) and len(re.match(ProteinRE, x).group(0)) == len(x)])
        return found

def findAllEnzymeNames(reactions, Map = {}, delay = 0.5):
    # For each specified reaction, finds the names of enzymes based on their EC identifiers via a web search
    # Saves the results into an optionally specified dictionary; delay is the time between consecutive queries
    allECIDs = sum([reaction.description['EC numbers'] for reaction in reactions if 'EC numbers' in reaction.description], [])
    allECIDs = sorted(list(set(allECIDs)))
    print(('There are ' + str(len(allECIDs)) + ' enzymes to process'))
    for ind, ECID in enumerate(allECIDs):
        if ind % 10 == 0:
            print(('Processed ' + str(ind) + ' enzymes so far'))
        Map[ECID] = findEnzymeName(ECID)
        time.sleep(delay)
    for reaction in reactions:
        if 'EC numbers' in reaction.description:
            curNumbers = reaction.description['EC numbers']
            names = [_f for _f in [Map[x] for x in curNumbers] if _f]
            reaction.description['Enzyme names'] = names
    return

def findEnzymeName(ECID):
    # This function returns the accepted name of an enzyme given by its E.C. number
    if not ECID:
        return ''
    else:
        url = 'http://www.expasy.org/enzyme/' + ECID
        page = urlopen(url).read()
        index = page.find("<strong>Accepted Name</strong>")
        if index == -1:
            return ''
        else:
            index1 = page.find("<strong>", index + len("<strong>"))
            index2 = page.find("</strong>", index1)
            return page[index1 + len("<strong>") : index2].strip().strip('.')

def getAllSpeciesInfo(species):
    # Finds the formula and CAS number for each given species via a web search
    # If information is already available for a given species, it has priority
    for ind, spec in enumerate(species):
        if ind % 10 == 0:
            print(('Processed ' + str(ind) + ' species so far'))
        curDescription = spec.description
        if curDescription is not None:
            descriptors = list(curDescription.keys())
            foundCAS, foundFormula = True, True
            if not ('CAS number' in descriptors and curDescription['CAS number']):
                foundCAS = False
            if not ('Formula' in descriptors and curDescription['Formula']):
                foundFormula = False
            if not (foundCAS and foundFormula):
                if 'Biocyc' in descriptors:
                    (formula, CAS) = getMoreInfo(curDescription['Biocyc'], 'Bio')
                    if formula:
                        if 'Formula' in descriptors and curDescription['Formula'] != formula:
                            print(('Error: discrepancy between the formula of ' + spec.name + ' and ' + formula + ' using KEGG ID'))
                        else:
                            curDescription['Formula'] = formula
                        foundFormula = True
                    if CAS:
                        if 'CAS number'  in descriptors and curDescription['CAS number'] != CAS:
                           print(('Error: discrepancy between the CAS number of ' + spec.name + ' and ' + CAS  + ' using KEGG ID'))
                        else:
                            curDescription['CAS number'] = CAS
                        foundCAS = True
                if not (foundCAS and foundFormula) and 'kegg-id' in descriptors:
                    (formula, CAS) = getMoreInfo(curDescription['kegg-id'], 'KEGG')
                    if formula:
                        if 'Formula' in descriptors and curDescription['Formula'] != formula:
                            print(('Error: discrepancy between the formula of ' + spec.name + ' and ' + formula + ' using Biocyc ID'))
                        else:
                            curDescription['Formula'] = formula
                    if CAS:
                        if 'CAS number'  in descriptors and curDescription['CAS number'] != CAS:
                           print(('Error: discrepancy between the CAS number of ' + spec.name + ' and ' + CAS + ' using Biocyc ID'))
                        else:
                            curDescription['CAS number'] = CAS
    return

def getMoreInfo(CompoundID, option = 'KEGG'):
    # This function returns the empirical formula and CAS number of a given compound from its ID.
    # The possible options are 'KEGG' if this is the KEGG ID and 'Bio' if it is the BioCyc ID.
    Formula, CASNum = '', ''
    if not CompoundID:
        return ('', '')
    elif option == 'KEGG':
        url = 'http://www.genome.jp/dbget-bin/www_bget?cpd:' + CompoundID
        page = urlopen(url).read()
        cleanPage = cleanupTags(page)
        if 'Formula' in cleanPage:
            index0 = cleanPage.index('Formula')
            Formula = cleanPage[index0 + 1]
        if 'CAS:' in cleanPage:
            index1 = cleanPage.index('CAS:')
            CASNum = re.match(CASRE, cleanPage[index1 + 1])
            if CASNum:
                CASNum = CASNum.group(0)
            else:
                CASNum = ''.group(0)
    elif option == 'Bio':
        url = 'http://biocyc.org/META/NEW-IMAGE?type=COMPOUND&object=' + CompoundID
        page = urlopen(url).read()
        index = page.find("Formula:")
        if index != -1:
            index1 = page.find("<p", index)
            Formula = page[index + len("Formula:") : index1].strip()
            Formula = Formula.replace('<SUB>','').replace('</SUB>','')
        index2 = page.find("CAS:")
        if index2 != -1:
            CASNum = re.match(CASRE, page[index2 + len("CAS:"):])
            if CASNum:
                CASNum = CASNum.group(0)
            else:
                CASNum = ''
    else:
        print(('Error: unrecognized option ' + option))    
    return (Formula, CASNum)

def convertAllFormulas(species):
    # For each given species, converts its formula into a dictionary format for further comparisons
    for ind, spec in enumerate(species):
        if ind % 10 == 0:
            print(('Processed ' + str(ind) + ' species so far'))
        if spec.description is not None and 'Formula' in spec.description and spec.description['Formula']:
            spec.description['Formula'] = convertFormula(spec.description['Formula'])
    return

def convertFormula(Formula):
    # Converts a string encoding a chemical formula into a dictionary with key = atom name, value = multiplicity
    # Note that strings such as (C2H5)n cannot be dealt with even if n is a specified integer, not just a letter
    Dico = {}
    if type(Formula) == type(Dico):
        return Formula
    else:
        # start by removing groups!
        other = [x for x in range(len(Formula)) if not (Formula[x].isalpha() or Formula[x].isdigit())]
        if other:
            Formula = Formula[:other[0]]
        for item in re.findall(atomRE, Formula):
            curLast = max([i for i,x in enumerate(item) if x.isalpha()]) + 1
            myIncrement(Dico, item[:curLast], int(item[curLast:]) if curLast < len(item) else 1)
        return Dico

def countFeatures(speciesList):
    # Tabulates the features present in a given list of species
    Features = {}
    for species in speciesList:
        if species.description is not None:
            curDescription = species.description
            for feature in curDescription:
                if feature not in Features:
                    Features[feature] = 0
                if curDescription[feature] or curDescription[feature] == 0:
                    Features[feature] += 1
    return Features

def convertFeaturesToMatrix(featureDicts, featureNames):
    # Converts a list of feature dictionaries into a matrix (using featureNames as header row)
    # Note: assumes that each feature in featureNames is present in each feature dictionary!
    Matrix = [[]]*len(featureDicts)
    for row, featureDict in enumerate(featureDicts):
        Matrix[row] = [0]*len(featureNames)
        for col, featureName in enumerate(featureNames):
            Matrix[row][col] = featureDict[featureName]
    return Matrix
