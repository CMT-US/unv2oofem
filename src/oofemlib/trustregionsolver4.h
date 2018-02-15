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

#include "nrsolver.h"

#include <slepceps.h>

namespace oofem {

class PetscSparseMtrx;

class TrustRegionSolver4 : public NRSolver {
public:
	TrustRegionSolver4(Domain * d, EngngModel * m);
	virtual ~TrustRegionSolver4();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual NM_Status solve(SparseMtrx &k, FloatArray &R, FloatArray *R0,
                            FloatArray &X, FloatArray &dX, FloatArray &F,
                            const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                            int &nite, TimeStep *);

    void updateTrustRegionSize(const double &iOldRes, const double &iNewtonTrialRes, const double &iEigTrialRes);

    double giveMaxAbs(const FloatArray &iVec) const;
    void clipToLimit(FloatArray &ioVec, const double &iLimit, double &oIncrementRatio);

    void calcTrialRes(double &oTrialRes, FloatArray &iX, FloatArray &idX, FloatArray &iddX_trial, TimeStep *tStep, int nite, FloatArray *iR, FloatArray &R, FloatArray &F, FloatArray &RT);

    bool checkConvergence(FloatArray &RT, FloatArray &F, FloatArray &rhs, FloatArray &ddX, FloatArray &X,
                          double RRT, const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange, bool printToScreen = true);

    void calcSmallestEigVal(double &oEigVal, FloatArray &oEigVec, PetscSparseMtrx &K);

    void addOnDiagonal(const double &iVal, PetscSparseMtrx &K);


protected:

    /// Trust region size
    double mTrustRegionSize;

    // Variables for eigenvalue analysis
    PetscSparseMtrx *A;
    PetscSparseMtrx *B;

    /// Eigenvalue solver context.
    EPS eps;
    /// Flag if context initialized.
    bool epsInit;

};

} /* namespace oofem */

#endif /* TRUSTREGIONSOLVER4_H_ */
