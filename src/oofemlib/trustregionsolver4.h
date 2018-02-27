/*
 * trustregionsolver4.h
 *
 *  Created on: Jun 20, 2017
 *      Author: svennine
 */

#ifndef TRUSTREGIONSOLVER4_H_
#define TRUSTREGIONSOLVER4_H_

#define _IFT_TrustRegionSolver4_Name "trustregionsolver4"

/// Initial size of trust region. The increment is restricted in L_inf norm.
#define _IFT_TrustRegionSolver4_InitialSize "initialsize"

/// Tolerance on smallest eigenvalue for adding perturbations
#define _IFT_TrustRegionSolver4_evtolpert "evtolpert"

#include "nrsolver.h"

#include <slepceps.h>

namespace oofem {

class PetscSparseMtrx;

class TrustRegionSolver4 : public NRSolver {
public:
	TrustRegionSolver4(Domain * d, EngngModel * m);
	virtual ~TrustRegionSolver4();

    virtual IRResultType initializeFrom(InputRecord *ir);

    NM_Status solve(SparseMtrx &k, FloatArray &R, FloatArray *R0,
                            FloatArray &X, FloatArray &dX, FloatArray &F,
                            const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                            int &nite, TimeStep *) override;

    void updateTrustRegionSize(const double &iOldRes, const double &iNewtonTrialRes, const double &iEigTrialRes);

    double giveMaxAbs(const FloatArray &iVec) const;
    void clipToLimit(FloatArray &ioVec, const double &iLimit, double &oIncrementRatio);

    void calcTrialRes(double &oTrialRes, FloatArray &iX, FloatArray &idX, FloatArray &iddX_trial, TimeStep *tStep, int nite, FloatArray *iR, FloatArray &R, FloatArray &F, FloatArray &RT);

    bool checkConvergence(FloatArray &RT, FloatArray &F, FloatArray &rhs, FloatArray &ddX, FloatArray &X,
                          double RRT, const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange, bool printToScreen = true);

    void checkPetscError(PetscErrorCode iErrorCode) const;
    void calcSmallestEigVal(double &oEigVal, FloatArray &oEigVec, PetscSparseMtrx &K);


protected:

    /// Trust region size
    double mTrustRegionSize;

    /// Perturbation tolerance: perturbations will be added if the smallest eignvalue is below this value.
    double mPertTol;

    /// Eigenvalue solver context.
    EPS eps;
    /// Flag if context initialized.
    bool epsInit;

};

} /* namespace oofem */

#endif /* TRUSTREGIONSOLVER4_H_ */
