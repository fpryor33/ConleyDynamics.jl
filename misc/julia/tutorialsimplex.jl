using ConleyDynamics

labels = ["A","B","C","D","E","F"]
simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
mvf = [["A","AC"],["B","AB"],["C","BC"],["D","BD"],["E","DE"]]

sc = create_simplicial_complex(labels,simplices)
fieldnames(typeof(sc))

# Create proper coordinates

coordA = [[0,0],[100,0],[50,100],[200,0],[300,0],[250,100]]
labeld = [3, 3, 1, 3, 3, 1]

fnameA = "/Users/wanner/Desktop/tutorialsimplex.svg"
fnameB = "/Users/wanner/Desktop/tutorialforman.svg"

plot_planar_simplicial(sc, coordA, fnameA,
                       labeldir=labeld, hfac=1.5, vfac=1.4)
plot_planar_simplicial(sc, coordA, fnameB, mvf=mvf,
                       labeldir=labeld, hfac=1.5, vfac=1.4)

