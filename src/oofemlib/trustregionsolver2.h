/*
 * trustregionsolver2.h
 *
 *  Created on: Mar 2, 2017
 *      Author: svennine
 */

#ifndef TRUSTREGIONSOLVER2_H_
#define TRUSTREGIONSOLVER2_H_

#if 1
#define _IFT_TrustRegionSolver2_Name "trustregionsolver2"

/// Initial size of trust region. The increment is restricted in L_inf norm.
#define _IFT_TrustRegionSolver2_InitialSize "initialsize"

/**
 * Parameter controlling the size of the eigenvector perturbation when a negative eigenvalue is detected.
 * In this case, the solution is perturbed by beta*trust_region_size_du.
 * Default value: beta = 0.05
 */
#define _IFT_TrustRegionSolver2_Beta "beta"

/**
 * If the eigenvector used for perturbation should be updated each iteration or kept for a number
 * of iterations. For example, eig_vec_recalc = 5 implies that the eigenvector is recomputed every fifth
 * iteration.
 */
#define _IFT_TrustRegionSolver2_EigVecRecompute "eig_vec_recalc"

/**
 * If the trust region size should be fixed or updated based on how well the solution converges.
 * For example, fix_size_iter = 10 implies that the trust region size is kept fixed for the first 10 iterations.
 */
#define _IFT_TrustRegionSolver2_FixSizeIter "fix_size_iter"

/**
 * Parameters for determining if an update was successful.
 */
#define _IFT_TrustRegionSolver2_eta1 "eta1"
#define _IFT_TrustRegionSolver2_eta2 "eta2"


#include "nrsolver.h"

#include <slepceps.h>


namespace oofem {

class PetscSparseMtrx;

/**
 *
 * Trust-region algorithm.
 * Implementation based on Conn et al. (2000).
 *
 */
class TrustRegionSolver2 : public NRSolver {
public:
	TrustRegionSolver2(Domain * d, EngngModel * m);
	virtual ~TrustRegionSolver2();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual NM_Status solve(SparseMtrx &k, FloatArray &R, FloatArray *R0,
                            FloatArray &X, FloatArray &dX, FloatArray &F,
                            const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                            int &nite, TimeStep *);

    bool checkConvergence(FloatArray &RT, FloatArray &F, FloatArray &rhs, FloatArray &ddX, FloatArray &X,
                          double RRT, const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange, bool printToScreen = true);

    void checkPetscError(PetscErrorCode iErrorCode) const;

    void calcSmallestEigVal(double &oEigVal, FloatArray &oEigVec, PetscSparseMtrx &K);

    void addOnDiagonal(const double &iVal, PetscSparseMtrx &K);

protected:

    /// Trust region parameters
    double mEta1, mEta2, mGamma1, mGamma2;

    /// Trust region size
    double mTrustRegionSize;

    /// How large the eigenvector perturbation should be
    double mBeta;

    /// How often the eigenvector should be recomputed
    int mEigVecRecalc;

    /// How many iterations the trust region size should be kept fixed.
    int mFixIterSize;

    // Variables for eigenvalue analysis
    PetscSparseMtrx *A;
    PetscSparseMtrx *B;

    /// Eigenvalue solver context.
    EPS eps;
    /// Flag if context initialized.
    bool epsInit;

};

} /* namespace oofem */

#endif
#endif /* TRUSTREGIONSOLVER2_H_ */
