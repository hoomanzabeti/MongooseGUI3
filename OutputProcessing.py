# This file contains functions for processing the output of the merging process
# Created by: Leonid Chindelevitch
# Last modified: January 13, 2012

import math, codecs
from Utilities import *

def CreateSMatrix(ReactList, m):
    # creates a stoichiometric matrix
    # from a given list of reactions,
    # assuming there are m metabolites.
    n = len(ReactList)
    # create a matrix of zeros
    M = []
    for i in range(m):
        M.append([0]*n)
    # fill it out with numerical values
    for k in range(len(ReactList)):
        Re = ReactList[k]
        for pair in Re:
            M[pair[0]][k] += pair[1]
    return M

def extractReactions(Mat):
    # Transforms a matrix into a list of metabolite-coefficient pairs (the reverse of CreateSMatrix)
    return [sorted([(x, Mat[x][y]) for x in range(len(Mat)) if Mat[x][y]], key = lambda x:x[1]) for y in range(len(Mat[0]))]

def WriteMatrix(M, Filename):
    # writes a matrix into a file
    # for subsequent use in MATLAB
    m = len(M)
    string  = '['
    string += ';'.join([','.join([str(x) for x in M[i]]) for i in range(m)])
    string += '];'
    f = open(Filename,'w')
    f.write(string)
    f.close()
    return

def ConvertToSBML(Reacts, growth, Irrev, Genes, MetabNames, External, ReactFeatures = [], MetabFeatures = [], GeneFeatures = [], sep = "\t", br = '\n', ModelName = 'Model', FileName = 'Model.xml'):
    # This function converts a given metabolic model into an SBML file in the level 2, version 4 format. See www.sbml.org.
    # The model contains a list of reactions in the usual format, the growth reaction, the list of all irreversible reactions,
    # the genes that catalyze every one of them in OR of ANDs format, the names of metabolites and a mapping(!) of external ones.
    # The current version assumes that there are only two compartments, intracellular and extracellular; this can be refined.
    # Optional inputs are React-, Metab- and Gene-Features, which are a table whose first row lists the attributes and then each
    # subsequent row contains the values of each attribute (empty or "" if not available) for each reaction, metabolite, or gene.
    # The last arguments are the separator and the linebreak to use (default: \t and \n), the model name and the file name to use.
    # begin the model
    Time = ConvertTime()
    allGenes = sorted(list(set(sum([sum(x,[]) for x in Genes],[]))))
    k = len(allGenes)
    l = len(External)
    m = len(MetabNames)
    n = len(Reacts)
    start  = ('<?xml version="1.0" encoding="UTF-8"?>' + br)
    start += ('<sbml xmlns="http://www.sbml.org/sbml/level2/version4" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/" level="2" version="4">' + br)
    start += (sep + '<model metaid="metaid_' + ModelName + '" name="' + ModelName + '">' + br)
    start += (sep*2 + '<annotation>' + br)
    start += (sep*3 + '<rdf:RDF>' + br)
    start += (sep*4 + '<rdf:Description rdf:about="#metaid_' + ModelName + '">' + br)
    start += (sep*5 + '<dc:creator rdf:parseType="Resource">' + br)
    start += (sep*6 + '<rdf:Bag>' + br)
    start += (sep*7 + '<rdf:li rdf:parseType="Resource">' + br)
    start += (sep*8 + '<vCard:N rdf:parseType="Resource">' + br)
    start += (sep*9 + '<vCard:Family>Chindelevitch</vCard:Family>' + br)
    start += (sep*9 + '<vCard:Given>Leonid</vCard:Given>' + br)
    start += (sep*8 + '</vCard:N>' + br)
    start += (sep*8 + '<vCard:EMAIL>leonidus@mit.edu</vCard:EMAIL>' + br)
    start += (sep*8 + '<vCard:ORG>' + br)
    start += (sep*9 + '<vCard:Orgname>Massachusetts Institute of Technology</vCard:Orgname>' + br)
    start += (sep*8 + '</vCard:ORG>' + br)
    start += (sep*7 + '</rdf:li>' + br)
    start += (sep*6 + '</rdf:Bag>' + br)
    start += (sep*5 + '</dc:creator>' + br)
    start += (sep*5 + '<dcterms:created rdf:parseType="Resource">' + br)
    start += (sep*6 + '<dcterms:W3CDTF>' + Time + '</dcterms:W3CDTF>' + br)
    start += (sep*5 + '</dcterms:created>' + br)
    start += (sep*5 + '<dcterms:modified rdf:parseType="Resource">' + br)
    start += (sep*6 + '<dcterms:W3CDTF>' + Time + '</dcterms:W3CDTF>' + br)
    start += (sep*5 + '</dcterms:modified>' + br)
    start += (sep*4 + '</rdf:Description>' + br)
    start += (sep*3 + '</rdf:RDF>' + br)
    start += (sep*2 + '</annotation>' + br)
    start += (sep*2 + '<listOfUnitDefinitions>' + br)
    start += (sep*3 + '<unitDefinition id="u" name="millimoles per gram (dry weight) per hour">' + br)
    start += (sep*4 + '<listOfUnits>' + br)
    start += (sep*5 + '<unit kind="mole" scale="-3"/>' + br)
    start += (sep*5 + '<unit kind="gram" exponent="-1"/>' + br)
    start += (sep*5 + '<unit kind="second" exponent="-1" multiplier="0.00027777"/>' + br)
    start += (sep*4 + '</listOfUnits>' + br)
    start += (sep*3 + '</unitDefinition>' + br)
    start += (sep*2 + '</listOfUnitDefinitions>' + br)
    print('Prepared the starting portion')
    # list the compartments
    comps  = (sep*2 + '<listOfCompartments>' + br)
    comps += (sep*3 + '<compartment metaid="metaid_c_01" id="c_01" name="cytoplasm">' + br)
    comps += (sep*3 + '</compartment>' + br)
    comps += (sep*3 + '<compartment metaid="metaid_c_02" id="c_02" name="extracellular">' + br)
    comps += (sep*3 + '</compartment>' + br)
    comps += (sep*2 + '</listOfCompartments>' + br)
    print('Prepared the compartments')
    # list the species types (internal metabolites and genes)
    L = int(math.ceil(math.log(k + m, 10)))
    # number of digits to use for the encoding
    if MetabFeatures:
        Attributes = MetabFeatures[0]
    specs  = (sep*2 + '<listOfSpeciesTypes>' + br)
    for i in range(m):
        specs += (sep*3 + '<speciesType metaid="metaid_t_' + Pad(i, L) + '" id="t_' + Pad(i, L) + '" name="' + MetabNames[i] + '">' + br)
        if MetabFeatures:
            specs += (sep*4 + '<notes>' + br)
            specs += (sep*5 + '<body xmlns="http://www.w3.org/1999/xhtml">' + br)
            for j in range(len(Attributes)):
                if Attributes[j].lower() != "name" and MetabFeatures[i+1][j]:
                    specs += (sep*6 + '<p>' + Attributes[j] + ':' + MetabFeatures[i+1][j] + '</p>' + br)
            specs += (sep*5 + '</body>' + br)
            specs += (sep*4 + '</notes>' + br)
        specs += (sep*3 + '</speciesType>' + br)
    print('Prepared the metabolites')
    if GeneFeatures:
        Attributes = GeneFeatures[0]
        GeneNames = [x[0] for x in GeneFeatures[1:]]
    for i in range(k):
        specs += (sep*3 + '<speciesType metaid="metaid_t_' + Pad(m + i, L) + '" id="t_' + Pad(m + i, L) + '" name="' + allGenes[i] + '">' + br)
        if GeneFeatures:
            ind = GeneNames.find(allGenes[i])
            if ind != -1:
                specs += (sep*4 + '<notes>' + br)
                specs += (sep*5 + '<body xmlns="http://www.w3.org/1999/xhtml">' + br)
                for j in range(1, len(Attributes)):
                    if GeneFeatures[i+1][j]:
                        specs += (sep*6 + '<p>' + Attributes[j] + ':' + GeneFeatures[i+1][j] + '</p>' + br)
                specs += (sep*5 + '</body>' + br)
                specs += (sep*4 + '</notes>' + br)
        specs += (sep*3 + '</speciesType>' + br)
    specs += (sep*2 + '</listOfSpeciesTypes>' + br)
    print('Prepared the genes')
    # list the species (all metabolites and genes)
    M = int(math.ceil(math.log(k + l + m, 10)))
    metab  = (sep*2 + '<listOfSpecies>' + br)
    for i in range(m):
        metab += (sep*3 + '<species id="s_' + Pad(i, M) + '" name="' + MetabNames[i] + '[cytoplasm]" speciesType="t_' + Pad(i, L) + '" compartment="c_01">' + br)
        metab += (sep*3 + '</species>' + br)
    for i in range(l):
        metab += (sep*3 + '<species id="s_' + Pad(m + i, M) + '" name="' + MetabNames[External[i]] + '[external]" speciesType="t_' + Pad(External[i], L) + '" compartment="c_02">' + br)
        metab += (sep*3 + '</species>' + br)
    # metab += (sep*2 + '</listOfSpecies>' + br) -> REMOVED (thanks to Johnjoe McFadden for pointing out the duplication)
    for i in range(k):
        metab += (sep*3 + '<species id="s_' + Pad(m + l + i, M) + '" name="' + allGenes[i] + '[cytoplasm]" speciesType="t_' + Pad(m + i, L) + '" compartment="c_01">' + br)
        metab += (sep*3 + '</species>' + br)
        # Mapping[m + l + i] = m + i
    metab += (sep*2 + '</listOfSpecies>' + br)
    print('Prepared the list of species')
    # list the reactions
    N = int(math.ceil(math.log(n, 10)))
    react  = (sep*2 + '<listOfReactions>' + br)
    for i in range(n):
        if ReactFeatures:
            Attributes = ReactFeatures[0]
        curReact = Reacts[i]
        curReags = [x for x in curReact if x[1] < 0]
        curProds = [x for x in curReact if x[1] > 0]
        curGenes = Genes[i]
        allCurGenes = set(sum(curGenes, []))
        react += (sep*3 + '<reaction id="r_' + Pad(i, N) + '"' + ' reversible="false"'*int(i in Irrev) + '>' + br)
        react += (sep*4 + '<notes>' + br)
        react += (sep*5 + '<body xmlns="http://www.w3.org/1999/xhtml">' + br)
        react += (sep*6 + '<p>GENE_ASSOCIATION:' + ' and '.join(['('*int(len(y) > 1) + ' or '.join(y) + ')'*int(len(y) > 1) for y in curGenes]) + '</p>' + br)
        if ReactFeatures:
            for j in range(len(Attributes)):
                if Attributes[j].lower() != "name" and ReactFeatures[i+1][j]:
                    react += (sep*6 + '<p>' + Attributes[j] + ':' + ReactFeatures[i+1][j] + '</p>' + br)
        react += (sep*5 + '</body>' + br)
        react += (sep*4 + '</notes>' + br)
        react += (sep*4 + '<listOfReactants>' + br)
        for pair in curReags:
            react += (sep*5 + '<speciesReference species="s_' + Pad(pair[0], M) + '"' + (' stoichiometry="' + str(abs(pair[1])) + '"')*int(pair[1] != -1.0) + '/>' + br)
        react += (sep*4 + '</listOfReactants>' + br)
        react += (sep*4 + '<listOfProducts>' + br)
        for pair in curProds:
            react += (sep*5 + '<speciesReference species="s_' + Pad(pair[0], M) + '"' + (' stoichiometry="' + str(abs(pair[1])) + '"')*int(pair[1] != 1.0) + '/>' + br)
        react += (sep*4 + '</listOfProducts>' + br)
        react += (sep*4 + '<listOfModifiers>' + br)
        for gene in allCurGenes:
            react += (sep*5 + '<modifierSpeciesReference species="s_' + Pad(m + l + allGenes.index(gene), M) + '"/>' + br)
        react += (sep*4 + '</listOfModifiers>' + br)
        react += (sep*4 + '<kineticLaw>' + br)    
        react += (sep*5 + '<listOfParameters>' + br)
        react += (sep*6 + '<parameter id="LOWER_BOUND" value="' + '-INF'*int(i not in Irrev) + '0'*int(i in Irrev) + '" units="u"/>' + br)
        react += (sep*6 + '<parameter id="UPPER_BOUND" value="INF" units="u"/>' + br)
        react += (sep*6 + '<parameter id="OBJECTIVE_COEFFICIENT" value="' + str(int(i == growth)) + '" units="dimensionless"/>' + br)
        react += (sep*5 + '</listOfParameters>' + br)
        react += (sep*4 + '</kineticLaw>' + br)
        react += (sep*3 + '</reaction>' + br)
    react += (sep*2 + '</listOfReactions>' + br)
    print('Prepared the reactions')
    # end the model
    final  = (sep + '</model>' + br)
    final += ('</sbml>' + br)
    print('Finished the model')
    f = codecs.open(FileName, encoding='utf-8', mode='w')
    f.write(start)
    f.write(comps)
    f.write(specs)
    f.write(metab)
    f.write(react)
    f.write(final)
    f.close()
    return allGenes
