from ModelParsing import *

model = parseSBML('Acinetobacter Baumannii.xml')
model.adjustCompartments('xt', startPos = -2)
model.biomassCoefficients[760] = 1
s = shelve.open('ParsedModel')
s['AB1'] = model
model.reduceNetwork()
# model.checkReduced()
# model.createMatrices()
# model.reducedMatrix
# model.applyReduction(model.reducedMatrix)
# model.checkElementalBalance()
# model.findBiomassReaction()
# model.findEssentialReactions()
# model.findIrrevBlockedReactions()

# print(get_used())

