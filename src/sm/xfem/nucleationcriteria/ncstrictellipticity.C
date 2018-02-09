/*
 * ncstrictellipticity.C
 *
 *  Created on: Mar 1, 2017
 *      Author: svennine
 */

#include "ncstrictellipticity.h"

#include "error.h"
#include "xfem/enrichmentitem.h"
#include "domain.h"
#include "element.h"
#include "gausspoint.h"

#include "Materials/structuralms.h"
#include "Materials/structuralmaterial.h"
#include "Materials/structuralfe2material.h"

#include "xfem/enrichmentitems/crack.h"
#include "xfem/xfemmanager.h"
#include "xfem/enrichmentfunction.h"
#include "xfem/enrichmentfronts/enrichmentfrontlinbranchfunconeel.h"
#include "xfem/enrichmentfronts/enrichmentfrontcohesivebranchfunconeel.h"
#include "xfem/enrichmentfronts/enrichmentfrontreducefront.h"
#include "xfem/propagationlaws/plhoopstresscirc.h"
#include "xfem/propagationlaws/plprincipalstrain.h"
#include "xfem/propagationlaws/plstrictellipticity.h"
#include "dynamicdatareader.h"
#include "dynamicinputrecord.h"
#include "geometry.h"
#include "classfactory.h"
#include "spatiallocalizer.h"
#include "crosssection.h"
#include "engngm.h"

#include <memory>
#include <cmath>

namespace oofem {
REGISTER_NucleationCriterion(NCStrictEllipticity)

NCStrictEllipticity::NCStrictEllipticity(Domain *ipDomain):
NucleationCriterion(ipDomain),
mInitialCrackLength(0.0),
mIncrementLength(1.0),
mPropStrainThreshold(0.0),
mCutOneEl(false),
mAllowInsInEnrEl(false)
{

}

NCStrictEllipticity::~NCStrictEllipticity() {

}

std::vector<std::unique_ptr<EnrichmentItem>> NCStrictEllipticity::nucleateEnrichmentItems() {

	printf("Entering NCStrictEllipticity::nucleateEnrichmentItems()\n");

	SpatialLocalizer *octree = this->mpDomain->giveSpatialLocalizer();
	XfemManager *xMan = mpDomain->giveXfemManager();

	std::vector<std::unique_ptr<EnrichmentItem>> eiList;

	// Center coordinates of newly inserted cracks
	std::vector<FloatArray> center_coord_inserted_cracks;


	// Loop over all elements and all bulk GP.
	for(auto &el : mpDomain->giveElements() ) {

		int numIR = el->giveNumberOfIntegrationRules();

		//int csNum = el->giveCrossSection()->giveNumber();

		if(true) { // check cross section index

			for(int irInd = 0; irInd < numIR; irInd++) {
				IntegrationRule *ir = el->giveIntegrationRule(irInd);

				int numGP = ir->giveNumberOfIntegrationPoints();

				for(int gpInd = 0; gpInd < numGP; gpInd++) {
					GaussPoint *gp = ir->getIntegrationPoint(gpInd);


//						StructuralMaterialStatus *ms = dynamic_cast<StructuralMaterialStatus*>(gp->giveMaterialStatus());
				        StructuralFE2MaterialStatus *ms = dynamic_cast< StructuralFE2MaterialStatus * >( gp->giveMaterialStatus() );

						if(ms != NULL) {



//							const FloatArray &strain = ms->giveTempStrainVector();


							FloatMatrix D9 = ms->giveTangent();

//					        printf("D9: "); D9.printYourself();

					        if( D9.giveNumberOfColumns() == 0 || D9.giveNumberOfRows() == 0 ) {
					        	TimeStep *tStep = mpDomain->giveEngngModel()->giveCurrentStep();
					        	ms->computeTangent(tStep);
					        	D9 = ms->giveTangent();
					        }

//					        printf("D9 recomputed: "); D9.printYourself();

							FloatMatrix D, Dsym;
					        StructuralMaterial::giveReducedSymMatrixForm(D, D9, _PlaneStress);
//					        printf("D: "); D.printYourself();

					        Dsym.beTranspositionOf(D);
					        Dsym.add(D);
					        Dsym.times(0.5);

							FloatArray eig_vals;
							FloatMatrix eig_vecs;
							int num_ev = Dsym.giveNumberOfColumns();
							Dsym.jaco_(eig_vals, eig_vecs, num_ev);


							double min_eig_val = eig_vals(0);
							int min_eig_val_index = 0;
							for(int i = 0; i < eig_vals.giveSize();i++) {
								double e = eig_vals(i);
								if(e < min_eig_val) {
									min_eig_val = e;
									min_eig_val_index = i;
								}
							}

							const double eig_val_tol = 0.0;

							FloatArray eig_vec;
							eig_vec.beColumnOf(eig_vecs, min_eig_val_index+1);
//							printf("eig_vec: "); eig_vec.printYourself();


							////////////////////////////////////////////////////////////////////////////
							// Check if there is more than one eigenvalue below the tolerance
							int num_ev_below_tol = 0;
							for(int i = 0; i < eig_vals.giveSize();i++) {
								double e = eig_vals(i);

								if( e < eig_val_tol ) {
									num_ev_below_tol++;
								}

							}

							if( num_ev_below_tol > 0 ) {
								printf("eig_vals: ");
								eig_vals.printYourself();

								printf("eig_vecs: ");
								eig_vecs.printYourself();
							}

							if( num_ev_below_tol > 1 ) {
								printf("num_ev_below_tol: %d\n", num_ev_below_tol );

								// If there is more than one negative eigenvalue,
								// we choose the one that is most "diagonally dominant",
								// i.e. an eigenvector dominated by mode I opening.
								// For now, measure with infinity norm.

								double max_mode_I = 0.0;
								int max_mode_I_ind = 0;

								for(int i = 0; i < eig_vals.giveSize();i++) {
									double e = eig_vals(i);

									if( e < eig_val_tol ) {

										FloatArray eig_vec_temp;
										eig_vec_temp.beColumnOf(eig_vecs, i+1);

										if( fabs(eig_vec_temp(0)) > max_mode_I ) {
											max_mode_I = fabs(eig_vec_temp(0));
											max_mode_I_ind = i;
										}

										if( fabs(eig_vec_temp(1)) > max_mode_I ) {
											max_mode_I = fabs(eig_vec_temp(1));
											max_mode_I_ind = i;
										}


										printf("e: %e\n", e);
										printf("eig_vec_temp: "); eig_vec_temp.printYourself();
									}
								}


								min_eig_val = eig_vals(max_mode_I_ind);
								eig_vec.beColumnOf(eig_vecs, max_mode_I_ind+1);
							}

							////////////////////////////////////////////////////////////////////////////


							if(min_eig_val < eig_val_tol) {

								// Find direction of softening

//								printf("min_eig_val: %e i: %d\n", min_eig_val, min_eig_val_index );
//								printf("\n\n//////////////////////////////////////////////////////////////////////////////\n");
//								printf("A discontinuity should be injected.\n\n\n\n");


								FloatArray crackNormal;
								if( findLocalizationDirectionFromStrain(crackNormal, eig_vec) ) {

									printf("Found localization direction.\n");


									FloatArray crackTangent = {-crackNormal(1), crackNormal(0)};
									crackTangent.normalize();
	//								printf("crackTangent: "); crackTangent.printYourself();

									// Create geometry
									FloatArray pc = {gp->giveGlobalCoordinates()(0), gp->giveGlobalCoordinates()(1)};
	//								printf("Global coord: "); pc.printYourself();


									FloatArray ps = pc;
									ps.add(-0.5*mInitialCrackLength, crackTangent);

									FloatArray pe = pc;
									pe.add(0.5*mInitialCrackLength, crackTangent);


									if(mCutOneEl) {
										// If desired, ensure that the crack cuts exactly one element.
										Line line(ps, pe);
										std::vector<FloatArray> intersecPoints;
										line.computeIntersectionPoints(el.get(), intersecPoints);


										if(intersecPoints.size() == 2) {
											ps = std::move(intersecPoints[0]);
											pe = std::move(intersecPoints[1]);

											FloatArray t;
											t.beDifferenceOf(pe, ps);

											pe.add( 1.0e-3, t);
											ps.add(-1.0e-3, t);

										}
										else {
	//										printf("intersecPoints.size(): %d\n", intersecPoints.size() );
											OOFEM_ERROR("intersecPoints.size() != 2")
										}

	//									printf("ps: "); ps.printYourself();
	//									printf("pc: "); pc.printYourself();
	//									printf("pe: "); pe.printYourself();
									}


		//							FloatArray principalVals;
		//							FloatMatrix principalDirs;
		//							StructuralMaterial::computePrincipalValDir(principalVals, principalDirs, strain, principal_strain);

									FloatArray points = {ps(0), ps(1), pc(0), pc(1), pe(0), pe(1)};


									// Check if nucleation is allowed, by checking for already existing cracks close to the GP.
									// Idea: Nucleation is not allowed if we are within an enriched element. In this way, branching is not
									// completely prohibited, but we avoid initiating multiple similar cracks.
									bool insertionAllowed = true;

#if 1
									Element *el_s = octree->giveElementContainingPoint(ps);
									if(el_s && !mCutOneEl) {
	//									if( xMan->isAllElementNodesEnriched(el_s) && !mAllowInsInEnrEl ) {
										if( xMan->isElementEnriched(el_s) && !mAllowInsInEnrEl ) {
											insertionAllowed = false;
										}
									}

									Element *el_c = octree->giveElementContainingPoint(pc);
									if(el_c) {
	//									if( xMan->isAllElementNodesEnriched(el_c) && !mAllowInsInEnrEl ) {
										if( xMan->isElementEnriched(el_c) && !mAllowInsInEnrEl ) {
											insertionAllowed = false;
										}
									}

									Element *el_e = octree->giveElementContainingPoint(pe);
									if(el_e && !mCutOneEl) {
	//									if( xMan->isAllElementNodesEnriched(el_e) && !mAllowInsInEnrEl ) {
										if( xMan->isElementEnriched(el_e) && !mAllowInsInEnrEl ) {
											insertionAllowed = false;
										}
									}
#endif

									for(const auto &x: center_coord_inserted_cracks) {
										if( x.distance(pc) <  2.0*mInitialCrackLength) {
											insertionAllowed = false;
											break;
											printf("Preventing insertion.\n");
										}
									}

									if(insertionAllowed) {
										int n = xMan->giveNumberOfEnrichmentItems() + 1;
										std::unique_ptr<Crack> crack(new Crack(n, xMan, mpDomain));


										// Geometry
										std::unique_ptr<BasicGeometry> geom = std::unique_ptr<BasicGeometry>(new PolygonLine());
										geom->insertVertexBack(ps);
										geom->insertVertexBack(pc);
										geom->insertVertexBack(pe);
										crack->setGeometry(std::move(geom));

										// Enrichment function
										EnrichmentFunction *ef = new HeavisideFunction(1, mpDomain);
										crack->setEnrichmentFunction(ef);

										// Enrichment fronts
		//									EnrichmentFront *efStart = new EnrFrontLinearBranchFuncOneEl();
										EnrichmentFront *efStart = new EnrFrontCohesiveBranchFuncOneEl();
//										EnrichmentFront *efStart = new EnrFrontReduceFront();
										crack->setEnrichmentFrontStart(efStart);

		//									EnrichmentFront *efEnd = new EnrFrontLinearBranchFuncOneEl();
										EnrichmentFront *efEnd = new EnrFrontCohesiveBranchFuncOneEl();
//										EnrichmentFront *efEnd = new EnrFrontReduceFront();
										crack->setEnrichmentFrontEnd(efEnd);




										///////////////////////////////////////
										// Propagation law

										// Options
	//									PLDoNothing *pl = new PLDoNothing();
	//									PLPrincipalStrain *pl = new PLPrincipalStrain();
	//									pl->setRadius(0.1*mIncrementLength);
	//									pl->setIncrementLength(mIncrementLength);
	//									pl->setStrainThreshold(mPropStrainThreshold);

										PLStrictEllipticity *pl = new PLStrictEllipticity();
										pl->setRadius(0.1*mIncrementLength);
										pl->setIncrementLength(mIncrementLength);


										crack->setPropagationLaw(pl);

										crack->updateDofIdPool();

										center_coord_inserted_cracks.push_back(pc);
										eiList.push_back( std::unique_ptr<EnrichmentItem>(std::move(crack)) );

										printf("NCStrictEllipticity: Nucleating a crack. min_eig_val: %e\n", min_eig_val );
										printf("crackNormal: "); crackNormal.printYourself();

										// We only introduce one crack per element in a single time step.
										break;
									}
								}
								else {
									printf("Could not find localization direction.\n");
								}

	#if 0

								if(principalVals[0] > 1.0) { // Check strain threshold



								}
	#endif

							}

					}
				}
			}

		} // If correct csNum
	}

	return std::move( eiList );

}


IRResultType NCStrictEllipticity::initializeFrom(InputRecord *ir) {

    IRResultType result; // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, mInitialCrackLength, _IFT_NCStrictEllipticity_InitialCrackLength);
    printf("mInitialCrackLength: %e\n", mInitialCrackLength);

    IR_GIVE_FIELD(ir, mIncrementLength, _IFT_NCStrictEllipticity_IncrementLength);
    printf("mIncrementLength: %e\n", mIncrementLength);

    IR_GIVE_FIELD(ir, mPropStrainThreshold, _IFT_NCStrictEllipticity_PropStrainThreshold);
    printf("mPropStrainThreshold: %e\n", mPropStrainThreshold);

    mCutOneEl = ir->hasField(_IFT_NCStrictEllipticity_CutOneEl);
    if(mCutOneEl) {
    	printf("Cutting one element.\n");
    }


    mAllowInsInEnrEl = ir->hasField(_IFT_NCStrictEllipticity_AllowInsertionInEnrEl);
    if(mAllowInsInEnrEl) {
    	printf("Allowing insertion in already enriched elements.\n");
    }

    return NucleationCriterion::initializeFrom(ir);
}

void NCStrictEllipticity :: appendInputRecords(DynamicDataReader &oDR)
{
    DynamicInputRecord *ir = new DynamicInputRecord();

    ir->setRecordKeywordField( this->giveInputRecordName(), 1 );

    ir->setField(mInitialCrackLength, _IFT_NCStrictEllipticity_InitialCrackLength);
    ir->setField(mIncrementLength, _IFT_NCStrictEllipticity_IncrementLength);
    ir->setField(mPropStrainThreshold, _IFT_NCStrictEllipticity_PropStrainThreshold);

    if(mCutOneEl) {
    	ir->setField(_IFT_NCStrictEllipticity_CutOneEl);
    }

    if(mAllowInsInEnrEl) {
    	ir->setField(_IFT_NCStrictEllipticity_AllowInsertionInEnrEl);
    }

    oDR.insertInputRecord(DataReader :: IR_crackNucleationRec, ir);

    // Enrichment function
    DynamicInputRecord *efRec = new DynamicInputRecord();
    mpEnrichmentFunc->giveInputRecord(* efRec);
    oDR.insertInputRecord(DataReader :: IR_enrichFuncRec, efRec);
}

bool NCStrictEllipticity :: findLocalizationDirectionFromStrain(FloatArray &oN, const FloatArray &iEps)
{
	return PLStrictEllipticity :: findLocalizationDirectionFromStrain(oN, iEps);
}



} /* namespace oofem */

