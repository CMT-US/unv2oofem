/*
 * plstrictellipticity.C
 *
 *  Created on: Apr 7, 2017
 *      Author: svennine
 */

#include "plstrictellipticity.h"
#include "xfem/propagationlaw.h"
#include "xfem/tipinfo.h"
#include "classfactory.h"
#include "mathfem.h"
#include "dynamicinputrecord.h"
#include "spatiallocalizer.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "xfem/enrichmentitem.h"
#include "feinterpol.h"
#include "xfem/xfemmanager.h"

#include "Materials/structuralms.h"
#include "Materials/structuralmaterial.h"
#include "Materials/structuralfe2material.h"

#include "engngm.h"

#include "xfem/XFEMDebugTools.h"

namespace oofem {
REGISTER_PropagationLaw(PLStrictEllipticity)

PLStrictEllipticity::PLStrictEllipticity():
mRadius(0.0),
mIncrementLength(0.0)
{

}

PLStrictEllipticity::~PLStrictEllipticity() {

}

IRResultType PLStrictEllipticity :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, mRadius,                          _IFT_PLStrictEllipticity_Radius);
    IR_GIVE_FIELD(ir, mIncrementLength,         _IFT_PLStrictEllipticity_IncLength);

    return IRRT_OK;
}

void PLStrictEllipticity :: giveInputRecord(DynamicInputRecord &input)
{
    int number = 1;
    input.setRecordKeywordField(this->giveInputRecordName(), number);

    input.setField(mRadius,                             _IFT_PLStrictEllipticity_Radius);
    input.setField(mIncrementLength,            _IFT_PLStrictEllipticity_IncLength);
}

bool PLStrictEllipticity :: propagateInterface(Domain &iDomain, EnrichmentFront &iEnrFront, TipPropagation &oTipProp)
{
	printf("Entering PLStrictEllipticity :: propagateInterface().\n");

    if ( !iEnrFront.propagationIsAllowed() ) {
        return false;
    }

    // Fetch crack tip data
    const TipInfo &tipInfo = iEnrFront.giveTipInfo();

    SpatialLocalizer *localizer = iDomain.giveSpatialLocalizer();


    const FloatArray &xT    = tipInfo.mGlobalCoord;
    const FloatArray &t     = tipInfo.mTangDir;
//    const FloatArray &n     = tipInfo.mNormalDir;

    // It is meaningless to propagate a tip that is not inside any element
    Element *el = localizer->giveElementContainingPoint(tipInfo.mGlobalCoord);
    if ( el != NULL ) {

    	FloatArray x(xT);
    	x.add(mRadius, t);


        std :: vector< double >sigTTArray, sigRTArray;

        FloatArray strainVec;

        // Take stress from closest Gauss point
        int region = 1;
        bool useCZGP = false;
        GaussPoint &gp = * ( localizer->giveClosestIP(x, region, useCZGP) );


        // Compute stresses
        StructuralFE2MaterialStatus *ms = dynamic_cast< StructuralFE2MaterialStatus * >( gp.giveMaterialStatus() );
        if ( ms == NULL ) {
        	OOFEM_ERROR("failed to fetch MaterialStatus.");
        }



        const FloatMatrix &D9 = ms->giveTangent();
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

							printf("eig_vals: ");
							eig_vals.printYourself();

//							printf("eig_vecs: ");
//							eig_vecs.printYourself();

		double min_eig_val = eig_vals(0);
		int min_eig_val_index = 0;
		for(int i = 0; i < eig_vals.giveSize();i++) {
			double e = eig_vals(i);
			if(e < min_eig_val) {
				min_eig_val = e;
				min_eig_val_index = i;
			}
		}

		const double eig_val_tol = 1.0e3;

		FloatArray eig_vec;
		eig_vec.beColumnOf(eig_vecs, min_eig_val_index+1);
//							printf("eig_vec: "); eig_vec.printYourself();


		if(min_eig_val < eig_val_tol) {

			// Find direction of softening

//								printf("min_eig_val: %e i: %d\n", min_eig_val, min_eig_val_index );
//								printf("\n\n//////////////////////////////////////////////////////////////////////////////\n");
//								printf("A discontinuity should be injected.\n\n\n\n");


			FloatArray crackNormal;
			findLocalizationDirectionFromStrain(crackNormal, eig_vec);



			FloatArray propTangent = {-crackNormal(1), crackNormal(0)};
			propTangent.normalize();
			printf("propTangent: "); propTangent.printYourself();


			if( propTangent.dotProduct(t) < 0.0 ) {
				propTangent.times(-1.0);
			}

			printf("Propagating crack.\n");
			// Fill up struct
			oTipProp.mTipIndex = tipInfo.mTipIndex;
			oTipProp.mPropagationDir = propTangent;
			oTipProp.mPropagationLength = mIncrementLength;

			return true;
		}

//        // Compare with threshold
//        printf("Max principal strain: %e\n", principalVals[0]);
//        if ( principalVals[0] > mStrainThreshold ) {
//
//			FloatArray propNormal;
//			propNormal.beColumnOf(principalDirs, 1);
//
//			FloatArray propTangent = {-propNormal(1), propNormal(0)};
//
//			if( propTangent.dotProduct(t) < 0.0 ) {
//				propTangent.times(-1.0);
//			}
//
//            // Fill up struct
//            oTipProp.mTipIndex = tipInfo.mTipIndex;
//            oTipProp.mPropagationDir = propTangent;
//            oTipProp.mPropagationLength = mIncrementLength;
//
//            return true;
//        }



    } // Tip is inside an element.

    return false;
}

bool PLStrictEllipticity :: findLocalizationDirectionFromStrain(FloatArray &oN, const FloatArray &iEps)
{

	// For given strain iEps, compute the corresponding localization direction n from the the minimization
	// problem of finding d,n such that
	// norm( iEps - (d \otimes n)^sym ) -> min

	// Define several initial guesses to solve for.
	// Then evaluate the objective function to find the best solution.
	FloatArray theta_initial_guesses;
	int N = 20;
	double d_theta = M_PI/double(N);
	double t = -0.5*M_PI;
	for(int i = 0; i < N; i++) {
//		printf("t: %e\n", t);
		theta_initial_guesses.append(t);
		t += d_theta;
	}

	double best_theta = 0.;
	double best_f = 1.e20;

	for( double theta0 : theta_initial_guesses ) {

		const int maxIter = 20;
		const double absTol = 1.0e-6;

		double d1 = cos(theta0), d2 = sin(theta0), theta = theta0;
		FloatArray x = {d1,d2,theta}, res(3), dx(3);

		double eps11 = iEps(0);
		double eps22 = iEps(1);
		double eps12 = 0.5*iEps(2);

		FloatMatrix K(3,3);

		double ct, st;

		// Do not allow larger increments than 5 degrees
		double maxThetaInc = 5.*M_PI/180.;

		for(unsigned int iter = 0; iter < maxIter; iter++) {

			d1 		= x(0);
			d2 		= x(1);
			theta 	= x(2);

			ct = cos(theta);
			st = sin(theta);

			res(0) = -2.*eps11*ct + 2.*d1 - 4.*eps12*st + 2.*d2*ct*st;
			res(1) = -4.*eps12*ct + 2.*d2 + 2.*d1*st*ct - 2.*eps22*st;
			res(2) = 2.*eps11*st*d1 + 4.*eps12*(st*d2 - ct*d1) - 2.*eps22*ct*d2 + 2.*d1*d2*(ct*ct - st*st);

			res.times(-1.);

			double resNorm = res.computeNorm();
//			printf("iter: %d res: %e\n", iter, resNorm );

			if(resNorm < absTol) {
//				printf("Solution converged.\n");
//				printf("theta (deg): %e, d1: %e, d2: %e\n", theta*180./M_PI, d1, d2);
				double n1 = cos(theta);
				double n2 = sin(theta);
//				printf("n1d1: %e 0.5*(n1d2+n2d1): %e n2d2: %e\n", n1*d1, 0.5*(n1*d2+n2*d1), n2*d2 );

				// Compute objective function
				FloatArray b = {n1*d1-iEps(0), n2*d2-iEps(1), n1*d2+n2*d1-iEps(2)};
//				printf("f: %e\n", b.computeNorm() );
				if( b.computeNorm() < best_f ) {
//					printf("Best solution so far.\n\n");
					best_f = b.computeNorm();
					best_theta = theta;
				}


				break;
			}

			K(0,0) = 2.0;
			K(0,1) = 2.*ct*st;
			K(0,2) = 2.*eps11*st - 4.*eps12*ct + 2.*d2*ct*ct - 2.*d2*st*st;

			K(1,0) = 2.*ct*st;
			K(1,1) = 2.;
			K(1,2) = 4.*eps12*st + 2.*d1*ct*ct - 2.*d1*st*st - 2.*eps22*ct;

			K(2,0) = 2.*eps11*st - 4.*eps12*ct + 2.*d2*ct*ct - 2.*d2*st*st;
			K(2,1) = 4.*eps12*st + 2.*d1*ct*ct - 2.*d1*st*st - 2.*eps22*ct;
			K(2,2) = 2.*eps11*d1*ct + 4.*eps12*d2*ct + 4.*eps12*d1*st + 2.*eps22*d2*st + 2.*d1*d2*( -4.*ct*st );

			// Solve and update

			K.solveForRhs(res, dx);
//			if(!K.solveForRhs(res, dx)) {
//				printf("Warning in NCStrictEllipticity::nucleateEnrichmentItems():\nFailed to solve for n.\n");
//			}


			if( fabs(dx(2)) > maxThetaInc ) {
				dx.times( maxThetaInc/fabs(dx(2)) );
			}

//			printf("dtheta (deg): %e\n", dx(2)*180./M_PI);

			x.add(dx);

		}

	}

//								double n1 = a(2);
//								double n2 = a(3);
//
//								printf("n1: %e n2: %e\n", n1, n2);


	oN = {cos(best_theta), sin(best_theta)};
//	printf("oN: "); oN.printYourself();

	// Check if the solution is sufficiently good
//	double rel_tol = 5.0e-2;
	double rel_tol = 1.0;
	printf("best_f: %e rel_tol*iEps.computeNorm(): %e\n", best_f, rel_tol*iEps.computeNorm() );
	if( best_f < rel_tol*iEps.computeNorm() ) {
		return true;
	}
	else {
		return false;
	}


}

} /* namespace oofem */
