import liboofem
dr=liboofem.OOFEMTXTDataReader("patch100.in")
pb=liboofem.InstanciateProblem(dr,liboofem.problemMode._processor,0)
pb.checkProblemConsistency()
pb.setRenumberFlag()
pb.solveYourself()
pb.terminateAnalysis()

fm=pb.giveContext().giveFieldManager()
print fm
# why is this one None?
print fm.giveField(liboofem.FieldType.FT_Displacements)
