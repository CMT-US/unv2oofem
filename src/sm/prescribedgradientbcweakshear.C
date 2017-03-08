/*
 * prescribedgradientbcweakshear.cpp
 *
 *  Created on: Mar 6, 2017
 *      Author: svennine
 */

#include "prescribedgradientbcweakshear.h"

#include "classfactory.h"
#include "domain.h"
#include "spatiallocalizer.h"
#include "gausspoint.h"
#include "timestep.h"
#include "function.h"

namespace oofem {
REGISTER_BoundaryCondition(PrescribedGradientBCWeakShear);

PrescribedGradientBCWeakShear::PrescribedGradientBCWeakShear(int n, Domain *d) :
PrescribedGradientBCWeakPeriodic(n, d)
{
	printf("Entering PrescribedGradientBCWeakShear::PrescribedGradientBCWeakShear(int n, Domain *d)\n");

    mTractionDofIDs = {Trac_v};
    mDispLockDofIDs = {LMP_v};
    mRegularDispDofIDs = {D_v};

}

PrescribedGradientBCWeakShear::~PrescribedGradientBCWeakShear() {

}

IRResultType PrescribedGradientBCWeakShear :: initializeFrom(InputRecord *ir)
{
    mMeshIsPeriodic = true;

    return PrescribedGradientBCWeak :: initializeFrom(ir);
}

void PrescribedGradientBCWeakShear :: postInitialize()
{
    bool enforceCornerPeriodicity = true;
    int numSides = 2;
    clear();
    createTractionMesh(enforceCornerPeriodicity, numSides);
}


void PrescribedGradientBCWeakShear :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                                CharType type, ValueModeType mode,
                                                const UnknownNumberingScheme &s, FloatArray *eNorm)
{
    int dim = domain->giveNumberOfSpatialDimensions();

    if ( type == ExternalForcesVector ) {
        // The external force vector is given by
        // f_ext = int N^trac H . (x - x_c)


        for ( TracSegArray* el : mpTracElNew ) {

            FloatArray contrib;
            computeExtForceElContrib(contrib, *el, dim, tStep);

            IntArray rows;
            el->giveTractionLocationArray(rows, type, s, {Trac_v});
            answer.assemble(contrib, rows);

        }

    } else if ( type == InternalForcesVector ) {

        for ( TracSegArray* el : mpTracElNew ) {

            for ( GaussPoint *gp: *(el->mIntRule.get()) ) {

                // Contribution on gamma_plus
                FloatArray contrib_disp, contrib_trac;
                IntArray disp_loc_array, trac_loc_array;
                computeIntForceGPContrib(contrib_disp, disp_loc_array, contrib_trac, trac_loc_array, *el, *gp, dim, tStep, gp->giveGlobalCoordinates(), 1.0, mode, type, s);
                answer.assemble(contrib_disp, disp_loc_array);
                answer.assemble(contrib_trac, trac_loc_array);


                // Contribution on gamma_minus
                contrib_disp.clear(); contrib_trac.clear();
                disp_loc_array.clear(); trac_loc_array.clear();
                FloatArray xMinus;
                this->giveMirroredPointOnGammaMinus(xMinus, gp->giveGlobalCoordinates());
                computeIntForceGPContrib(contrib_disp, disp_loc_array, contrib_trac, trac_loc_array, *el, *gp, dim, tStep, xMinus, -1.0, mode, type, s);
                answer.assemble(contrib_disp, disp_loc_array);
                answer.assemble(contrib_trac, trac_loc_array);

            }
    	}



        if ( mpDisplacementLock != NULL ) {
            IntArray dispLockRows;
            mpDisplacementLock->giveLocationArray(giveDispLockDofIDs(), dispLockRows, s);

            FloatArray fe_dispLock;

            int lockNodePlaceInArray = domain->giveDofManPlaceInArray(mLockNodeInd);
            FloatArray nodeUnknowns;
            domain->giveDofManager(lockNodePlaceInArray)->giveUnknownVector(nodeUnknowns,this->giveRegularDispDofIDs(), mode, tStep);

            for ( int i = 0; i < domain->giveNumberOfSpatialDimensions(); i++ ) {
                fe_dispLock.push_back(nodeUnknowns [ i ]);
            }

            fe_dispLock.times(mDispLockScaling);

            answer.assemble(fe_dispLock, dispLockRows);
        }

    }
}

void PrescribedGradientBCWeakShear :: computeExtForceElContrib(FloatArray &oContrib, TracSegArray &iEl, int iDim, TimeStep *tStep)
{
	printf("Entering PrescribedGradientBCWeakShear :: computeExtForceElContrib.\n");

    oContrib.clear();
    FloatArray contrib_gp;

    for ( GaussPoint *gp: *(iEl.mIntRule.get()) ) {

        // Fetch global coordinate x
        const FloatArray &x = gp->giveGlobalCoordinates();

        // Compute H.[x]
        FloatArray temp;
        giveBoundaryCoordVector(temp, x);

        FloatArray Hx;
        FloatMatrix grad2D = mGradient;
        grad2D.resizeWithData(2,2);
        Hx.beProductOf(grad2D, temp);

        // For now, assume piecewise constant approx
        FloatArray Ntrac = FloatArray { 1.0*mTracDofScaling };

        // N-matrix
        FloatMatrix Nmat;
        Nmat.beNMatrixOf(Ntrac, iDim);


        // Assemble contribution to the global vector directly
        contrib_gp.beTProductOf(Nmat, Hx);
        double detJ = 0.5 * iEl.giveLength();
        double loadLevel = this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
        contrib_gp.times(-detJ * gp->giveWeight()*loadLevel);

        oContrib.add(contrib_gp);
    }

}

void PrescribedGradientBCWeakShear :: giveBoundaryCoordVector(FloatArray &oX, const FloatArray &iPos) const
{
	FloatArray xC = mLC;
	xC.add(mUC);
	xC.times(0.5);

    FloatArray xMinus;
    giveMirroredPointOnGammaMinus(xMinus, iPos);

    oX = {
        iPos [ 0 ] - xMinus [ 0 ] - xC[0], iPos [ 1 ] - xMinus [ 1 ] - xC[1]
    };
}

bool PrescribedGradientBCWeakShear :: pointIsOnGammaPlus(const FloatArray &iPos) const
{
    const double distTol = 1.0e-12;

    if ( iPos [ 0 ] > mUC [ 0 ] - distTol ) {
        return true;
    }

    return false;
}

bool PrescribedGradientBCWeakShear ::pointIsOnGammaR(const FloatArray &iPos) const
{
	printf("mUC: "); mUC.printYourself();

	const double distTol = 1.0e-12;

	if ( iPos [ 0 ] > mUC [ 0 ] - distTol ) {
		return true;
	}
	else {
		return false;
	}

}


void PrescribedGradientBCWeakShear :: removeSupressedSegments(TracSegArray &ioTSeg, const double &iAbsTol)
{
	// Idea:	Loop over segments and check if the currect segment is located on
	//			the horizontal or vertical part of the boundary. Only keep segments
	//			on the vertical part of the boundary.


//	printf("Entering PrescribedGradientBCWeakShear :: removeSupressedSegments()\n");

    std :: vector< Line > tmp;

    SpatialLocalizer *localizer = domain->giveSpatialLocalizer();
    const double tol2 = iAbsTol*iAbsTol;

    for ( auto &l : ioTSeg.mInteriorSegments ) {
        const FloatArray &xS = l.giveVertex(1);
        const FloatArray &xE = l.giveVertex(2);
        FloatArray xPlus = {0.5*(xS[0]+xE[0]), 0.5*(xS[1]+xE[1])};
        printf("xPlus: "); xPlus.printYourself();

        if ( pointIsOnGammaR(xPlus) ) {
            tmp.push_back(l);
        }
    }

    ioTSeg.mInteriorSegments = std::move(tmp);

//    printf("done.\n");


}


} /* namespace oofem */
