import oofempy

# engngModel
problem = oofempy.linearStatic(nSteps=3, outFile="test2.out")

# domain (if no engngModel specified to domain, it is asigned to the last one created)
domain = oofempy.domain(1, 1, problem, oofempy.domainType._2dBeamMode, tstep_all=True, dofman_all=True, element_all=True)
problem.setDomain(1, domain, True)

ltf1 = oofempy.peakFunction(1, domain, t=1, f_t=1)
ltf2 = oofempy.peakFunction(2, domain, t=2, f_t=1)
ltf3 = oofempy.peakFunction(3, domain, t=3, f_t=1)
ltfs = (ltf1, ltf2, ltf3)

# boundary conditions
# loadTimeFunction parameter can be specified as int value or as LoadTimeFunction itself (valid for all objects with giveNumber() method)
bc1   = oofempy.boundaryCondition(    1, domain, loadTimeFunction=1,    prescribedValue=0.0)
bc2   = oofempy.boundaryCondition(    2, domain, loadTimeFunction=2,    prescribedValue=-.006e-3)
eLoad = oofempy.constantEdgeLoad(     3, domain, loadTimeFunction=1, components=(0.,10.,0.), loadType=3, ndofs=3)
nLoad = oofempy.nodalLoad(            4, domain, loadTimeFunction=1,    components=(-18.,24.,0.))
tLoad = oofempy.structTemperatureLoad(5, domain, loadTimeFunction=3, components=(30.,-20.))
bcs = (bc1, bc2, eLoad, nLoad, tLoad)

# nodes
# if one value is passed as parameter where oofem expects array of values, it must be passed as tuple or list (see load in n4)
n1 = oofempy.node(1, domain, coords=(0.,  0., 0. ), bc=(0,1,0))
n2 = oofempy.node(2, domain, coords=(2.4, 0., 0. ), bc=(0,0,0))
n3 = oofempy.node(3, domain, coords=(3.8, 0., 0. ), bc=(0,0,bc1))
n4 = oofempy.node(4, domain, coords=(5.8, 0., 1.5), bc=(0,0,0), load=(4,))
n5 = oofempy.node(5, domain, coords=(7.8, 0., 3.0), bc=(0,1,0))
n6 = oofempy.node(6, domain, coords=(2.4, 0., 3.0), bc=(bc1,1,bc2))
nodes = (n1, n2, n3, n4, n5, n6)

# material and cross section
mat = oofempy.isoLE(1, domain, d=1., E=30.e6, n=0.2, tAlpha=1.2e-5)
cs  = oofempy.simpleCS(1, domain, area=0.162, Iy=0.0039366, beamShearCoeff=1.e18, thick=0.54)

# elements
e1 = oofempy.beam2d(1, domain, nodes=(1,n2),  mat=1,   crossSect=1,  boundaryLoads=(3,1), bodyLoads=(5,))
e2 = oofempy.beam2d(2, domain, nodes=(2,3),   mat=mat, crossSect=1,  DofsToCondense=(6,), bodyLoads=[tLoad])
e3 = oofempy.beam2d(3, domain, nodes=(n3,4),  mat=1,   crossSect=cs, dofstocondense=[3])
e4 = oofempy.beam2d(4, domain, nodes=(n4,n5), mat=mat, crossSect=cs)
e5 = oofempy.beam2d(5, domain, nodes=(n6,2),  mat=1,   crossSect=1,  DofsToCondense=(6,))
elems = (e1, e2, e3, e4, e5)

# add eveything to domain (resize container first to save some time, but it is not necessary 0 see ltfs)
domain.resizeDofManagers(len(nodes))
for n in nodes:
   domain.setDofManager(n.number, n)
domain.resizeElements(len(elems))
for e in elems:
   domain.setElement(e.number, e)
domain.resizeMaterials(1)
domain.setMaterial(1, mat)
domain.resizeCrossSectionModels(1)
domain.setCrossSection(1, cs)
domain.resizeBoundaryConditions(len(bcs))
for bc in bcs:
   domain.setBoundaryCondition(bc.number, bc)
domain.resizeFunctions(len(ltfs))
for ltf in ltfs:
   domain.setFunction(ltf.number, ltf)


print("\nSolving problem")
problem.checkProblemConsistency()
problem.init()
problem.postInitialize()
problem.setRenumberFlag()
problem.solveYourself()
problem.terminateAnalysis()
print("\nProblem solved")
