from TasksThreeFunctions import SIRV

# We want to compare two proportons
S, I, R, V = 0, 1, 2, 3

# We want to see the effects of improving the vaccine:
newVaccinationRates = [.005, .01, .05, .25, .5]
newVaccineImunityLoss = [.005, 0.025, 0.01, .1, .2]

SIRV(0.95, 0.05, 0, 0, 0.2, .1, .05, .0001, .005, 200, 
     mainPlot = False,
     plotChoice1 = S, 
     plotChoice2 = V, 
     newNu = newVaccinationRates, 
     newLoss = newVaccineImunityLoss
     )

