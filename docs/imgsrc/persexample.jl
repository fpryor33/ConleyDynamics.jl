using ConleyDynamics

# Creation of the filtration figures

labels1    = ["A","B","C","D","E","G"]
simplices1 = [["A","B"],["B","C","G"],["C","D"],["D","E"]]
coords1    = [[0,0],[2,0],[4,0],[6,0],[8,0],[4,2]]
ldir1      = [3,3,3,3,3,1]
fname1     = "/Users/wanner/Desktop/persistence1.pdf"
sc1 = create_simplicial_complex(labels1,simplices1)
plot_planar_simplicial(sc1,coords1,fname1,labeldir=ldir1,labeldis=10,hfac=2,vfac=1.5,sfac=50)

labels2    = ["A","B","C","D","E","F","G","H"]
simplices2 = [["A","B"],["A","F"],["B","F"],["B","C","G"],["H"],["C","D"],["D","E"]]
coords2    = [[0,0],[2,0],[4,0],[6,0],[8,0],[1,2],[4,2],[6,2]]
ldir2      = [3,3,3,3,3,1,1,1]
fname2     = "/Users/wanner/Desktop/persistence2.pdf"
sc2 = create_simplicial_complex(labels2,simplices2)
plot_planar_simplicial(sc2,coords2,fname2,labeldir=ldir2,labeldis=10,hfac=2,vfac=1.5,sfac=50)

labels3    = ["A","B","C","D","E","F","G","H"]
simplices3 = [["A","B"],["A","F"],["B","F"],["B","C","G"],["C","D"],["D","E"],["D","H"],["E","H"]]
coords3    = [[0,0],[2,0],[4,0],[6,0],[8,0],[1,2],[4,2],[6,2]]
ldir3      = [3,3,3,3,3,1,1,1]
fname3     = "/Users/wanner/Desktop/persistence3.pdf"
sc3 = create_simplicial_complex(labels3,simplices3)
plot_planar_simplicial(sc3,coords3,fname3,labeldir=ldir3,labeldis=10,hfac=2,vfac=1.5,sfac=50)

labels4    = ["A","B","C","D","E","F","G","H"]
simplices4 = [["A","B"],["A","F"],["B","F"],["B","C","G"],["D","E","H"],["C","D"],["G","H"]]
coords4    = [[0,0],[2,0],[4,0],[6,0],[8,0],[1,2],[4,2],[6,2]]
ldir4      = [3,3,3,3,3,1,1,1]
fname4     = "/Users/wanner/Desktop/persistence4.pdf"
sc4 = create_simplicial_complex(labels4,simplices4)
plot_planar_simplicial(sc4,coords4,fname4,labeldir=ldir4,labeldis=10,hfac=2,vfac=1.5,sfac=50)

# Computation of the persistence intervals

labels    = ["A","B","C","D","E","F","G","H"]
simplices = [["A","B"],["A","F"],["B","F"],["B","C","G"],["D","E","H"],["C","D"],["G","H"]]
sc = create_simplicial_complex(labels,simplices)

filtration = [1,1,1,1,1,2,1,2,
              1,2,1,2,1,1,1,1,3,3,4,
              1,4]
ph = persistent_homology(sc, filtration)

length.(ph[1])







# Example for a subcomplex

tmpfil = fill(Int(0),sc.ncells)
tmpfil[sc.indices["BC"]]  = 1
tmpfil[sc.indices["BG"]]  = 1
tmpfil[sc.indices["CG"]]  = 1
tmpfil[sc.indices["DE"]]  = 2
tmpfil[sc.indices["DH"]]  = 2
tmpfil[sc.indices["EH"]]  = 2
tmpfil[sc.indices["CD"]]  = 3
tmpfil[sc.indices["DEH"]] = 4
tmpfil[sc.indices["BCG"]] = 5
scsub, filtrationsub = lefschetz_filtration(sc, tmpfil)
phsub = persistent_homology(scsub, filtrationsub)

length.(phsub[1])

