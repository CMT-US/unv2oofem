/*
 * trustregionsolver3.h
 *
 *  Created on: Mar 20, 2017
 *      Author: svennine
 */

#ifndef TRUSTREGIONSOLVER3_H_
#define TRUSTREGIONSOLVER3_H_

#if 1
#define _IFT_TrustRegionSolver3_Name "trustregionsolver3"
#define _IFT_TrustRegionSolver3_InitialSize "initialsize"
#define _IFT_TrustRegionSolver3_Beta "beta"
#define _IFT_TrustRegionSolver3_EigVecRecompute "eig_vec_recalc"


#include "nrsolver.h"

#include <slepceps.h>


namespace oofem {

class PetscSparseMtrx;

/**
 *
 * Trust-region algorithm.
 * Special version for weakly periodic BCs.
 * Implementation based on Conn et al. (2000).
 *
 */
class TrustRegionSolver3 : public NRSolver {
public:
	TrustRegionSolver3(Domain * d, EngngModel * m);
	virtual ~TrustRegionSolver3();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual NM_Status solve(SparseMtrx &k, FloatArray &R, FloatArray *R0,
                            FloatArray &X, FloatArray &dX, FloatArray &F,
                            const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                            int &nite, TimeStep *) override;

    bool checkConvergence(FloatArray &RT, FloatArray &F, FloatArray &rhs, FloatArray &ddX, FloatArray &X,
                          double RRT, const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange, bool printToScreen = true);

    void checkPetscError(PetscErrorCode iErrorCode) const;
    void calcSmallestEigVal(double &oEigVal, FloatArray &oEigVec, PetscSparseMtrx &K);

    virtual const char *giveClassName() const { return "TrustRegionSolver3"; }
    virtual const char *giveInputRecordName() const { return _IFT_TrustRegionSolver3_Name; }

protected:

    /// Trust region parameters
    double mEta1, mEta2, mGamma1, mGamma2;

    /// Trust region size
    double mTrustRegionSize;

    // How large the eigenvector perturbation should be
    double mBeta;

    // How often the eigenvector should be recomputed
    int mEigVecRecalc;

    /// Eigenvalue solver context.
    EPS eps;
    /// Flag if context initialized.
    bool epsInit;

};

} /* namespace oofem */

#endif
#endif /* TRUSTREGIONSOLVER3_H_ */
