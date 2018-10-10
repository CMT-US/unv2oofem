/*
 * linesearchsolver.C
 *
 *  Created on: Feb 24, 2017
 *      Author: svennine
 */

#include "trustregionsolver.h"

#include "nrsolver.h"
#include "classfactory.h"
#include "dof.h"
#include "unknownnumberingscheme.h"
#include "dofmanager.h"
#include "node.h"
#include "element.h"
#include "generalboundarycondition.h"

#include "timestep.h"
#include "engngm.h"
#include "domain.h"
#include "exportmodulemanager.h"

#if 0
#include "petscsparsemtrx.h"
#include <petscvec.h>
#endif

namespace oofem {
#define nrsolver_ERROR_NORM_SMALL_NUM 1.e-6
#define NRSOLVER_MAX_REL_ERROR_BOUND 1.e20
#define NRSOLVER_MAX_RESTARTS 4
#define NRSOLVER_RESET_STEP_REDUCE 0.25
#define NRSOLVER_DEFAULT_NRM_TICKS 10


REGISTER_SparseNonLinearSystemNM(TrustRegionSolver)


TrustRegionSolver::TrustRegionSolver(Domain * d, EngngModel * m) :
NRSolver(d,m),
mEta1(-100.0),
mEta2(0.9),
mGamma1(0.5),
mGamma2(0.5),
mTrustRegionSize(1.0e-3)
{

//	mEta1(0.01),
//	mEta2(0.9),

}

TrustRegionSolver::~TrustRegionSolver() {
	// TODO Auto-generated destructor stub
}



#if 1
NM_Status
TrustRegionSolver :: solve(SparseMtrx &k, FloatArray &R, FloatArray *R0,
                  FloatArray &X, FloatArray &dX, FloatArray &F,
                  const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                  int &nite, TimeStep *tStep)
//
// this function solve the problem of the unbalanced equilibrium
// using NR scheme
//
//
{
#if 0

    // residual, iteration increment of solution, total external force
    FloatArray rhs, ddX, RT;
    double RRT;
    int neq = X.giveSize();
    bool converged, errorOutOfRangeFlag;
    ParallelContext *parallel_context = engngModel->giveParallelContext( this->domain->giveNumber() );

    if ( engngModel->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO("NRSolver: Iteration");
        if ( rtolf.at(1) > 0.0 ) {
            OOFEM_LOG_INFO(" ForceError");
        }
        if ( rtold.at(1) > 0.0 ) {
            OOFEM_LOG_INFO(" DisplError");
        }
        OOFEM_LOG_INFO("\n----------------------------------------------------------------------------\n");
    }

    l = 1.0;

    NM_Status status = NM_None;
    this->giveLinearSolver();

    // compute total load R = R+R0
    RT = R;
    if ( R0 ) {
        RT.add(* R0);
    }

    RRT = parallel_context->localNorm(RT);
    RRT *= RRT;

    ddX.resize(neq);
    ddX.zero();

    // Fetch the matrix before evaluating internal forces.
    // This is intentional, since its a simple way to drastically increase convergence for nonlinear problems.
    // (This old tangent is just used)
    // This improves convergence for many nonlinear problems, but not all. It may actually
    // cause divergence for some nonlinear problems. Therefore a flag is used to determine if
    // the stiffness should be evaluated before the residual (default yes). /ES

    engngModel->updateComponent(tStep, NonLinearLhs, domain);
    if ( this->prescribedDofsFlag ) {
        if ( !prescribedEqsInitFlag ) {
            this->initPrescribedEqs();
        }
        applyConstraintsToStiffness(k);
    }

    double coeff = 0.0e3;

    nite = 0;
    for ( nite = 0; ; ++nite ) {
        // Compute the residual
        engngModel->updateComponent(tStep, InternalRhs, domain);
        if (nite || iR == NULL) {
            rhs.beDifferenceOf(RT, F);
        } else {
            rhs = R;
            if (iR) {
                rhs.add(*iR); // add initial guess
            }
        }
        if ( this->prescribedDofsFlag ) {
            this->applyConstraintsToLoadIncrement(nite, k, rhs, rlm, tStep);
        }

        // convergence check
        converged = this->checkConvergence(RT, F, rhs, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag);
        double oldForceErr = rhs.computeNorm();
        if ( engngModel->giveProblemScale() == macroScale ) {
        	printf("oldForceErr: %e\n", oldForceErr );
        }

        if ( errorOutOfRangeFlag ) {
            status = NM_NoSuccess;
            OOFEM_WARNING("Divergence reached after %d iterations", nite);
            break;
        } else if ( converged && ( nite >= minIterations ) ) {
            status |= NM_Success;
            break;
        } else if ( nite >= nsmax ) {
            OOFEM_LOG_DEBUG("Maximum number of iterations reached\n");
            break;
        }


        //#ifdef __PETSC_MODULE
        #if 0
            PetscSparseMtrx *lhs = dynamic_cast< PetscSparseMtrx * >(&k);
            if ( lhs ) {
            	printf("Fetched Petsc matrix.\n");

            	int N = k.giveNumberOfRows();


            	// Monitor eigenvalues
            	//PetscErrorCode  KSPComputeEigenvalues(KSP ksp,PetscInt n,PetscReal r[],PetscReal c[],PetscInt *neig)
            	int n = 10;
            	int neig = 6;
//            	Vec eig_real;
//                VecCreate(PETSC_COMM_SELF, & eig_real);
//                VecSetType(eig_real, VECSEQ);
//                VecSetSizes(eig_real, PETSC_DECIDE, n);
//
//            	Vec eig_imag;
//                VecCreate(PETSC_COMM_SELF, & eig_imag);
//                VecSetType(eig_imag, VECSEQ);
//                VecSetSizes(eig_imag, PETSC_DECIDE, n);

#if 0
                if ( !lhs->kspInit ) {
                	printf("!lhs->kspInit.\n");
                }
                else {
                	printf("Computing eigenvalues.\n");

					PetscReal eig_real[n];
					PetscReal eig_imag[n];

					PetscErrorCode err = KSPComputeEigenvalues(lhs->ksp, n, eig_real, eig_imag, &neig);
                }
#endif

            	Vec petsc_mat_diag;
                VecCreate(PETSC_COMM_SELF, & petsc_mat_diag);
                VecSetType(petsc_mat_diag, VECSEQ);
                VecSetSizes(petsc_mat_diag, PETSC_DECIDE, N);

            	MatGetDiagonal(lhs->mtrx, petsc_mat_diag);


        //    	VecScale(petsc_mat_diag, coeff);
//            	VecScale(petsc_mat_diag, 1.0 + coeff);
            	for(int i = 0; i < N; i++) {
            		VecSetValue(petsc_mat_diag, i, coeff, ADD_VALUES);
            	}
            	printf("coeff: %e\n", coeff );
            	coeff *= 0.9;

            	MatDiagonalSet(lhs->mtrx, petsc_mat_diag, INSERT_VALUES);

                VecDestroy(& petsc_mat_diag);

        //    	FloatArray diag_scale_vec(N);
        //    	for(int i = 0; i < N; i++) {
        //////    		k.at(i,i) *= coeff;
        //    		petsc_mat_diag[i] *= coeff;
        //    	}


        //    	MatDiagonalSet(lhs->mtrx, petsc_mat_diag, INSERT_VALUES);

        //        Vec petsc_diag_scale_vec;
        //        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, diag_scale_vec.giveSize(), diag_scale_vec.givePointer(), & petsc_diag_scale_vec);
        //
        //        MatDiagonalScaleLocal(lhs->mtrx, petsc_diag_scale_vec);

            }
        #endif


        if ( nite > 0 || !mCalcStiffBeforeRes ) {
            if ( ( NR_Mode == nrsolverFullNRM ) || ( ( NR_Mode == nrsolverAccelNRM ) && ( nite % MANRMSteps == 0 ) ) ) {
                engngModel->updateComponent(tStep, NonLinearLhs, domain);
                applyConstraintsToStiffness(k);
            }
        }

        if ( ( nite == 0 ) && ( deltaL < 1.0 ) ) { // deltaL < 1 means no increment applied, only equilibrate current state
            rhs.zero();
            R.zero();
            ddX = rhs;
        } else {

//            if ( engngModel->giveProblemScale() == macroScale ) {
//            	k.writeToFile("k.txt");
//            }

            linSolver->solve(k, rhs, ddX);
        }

        //
        // update solution
        //
        if ( this->lsFlag && ( nite > 0 ) ) { // Why not nite == 0 ?
            // line search
            LineSearchNM :: LS_status LSstatus;
            double eta;
            this->giveLineSearchSolver()->solve(X, ddX, F, R, R0, prescribedEqs, 1.0, eta, LSstatus, tStep);
        } else if ( this->constrainedNRFlag && ( nite > this->constrainedNRminiter ) ) {
            ///@todo This doesn't check units, it is nonsense and must be corrected / Mikael
            if ( this->forceErrVec.computeSquaredNorm() > this->forceErrVecOld.computeSquaredNorm() ) {
                OOFEM_LOG_INFO("Constraining increment to be %e times full increment...\n", this->constrainedNRalpha);
                ddX.times(this->constrainedNRalpha);
            }
            //this->giveConstrainedNRSolver()->solve(X, & ddX, this->forceErrVec, this->forceErrVecOld, status, tStep);
        }


        /////////////////////////////////////////

        double maxInc = 0.0;
        for ( double inc : ddX ) {
            if(fabs(inc) > maxInc) {
                maxInc = fabs(inc);
            }
        }

        if(maxInc > maxIncAllowed) {
            if ( engngModel->giveProblemScale() == macroScale ) {
            	printf("Restricting increment. maxInc: %e\n", maxInc);
            }
        	ddX.times(maxIncAllowed/maxInc);
        }

        /////////////////////////////////////////


        FloatArray steepest_descent = rhs;
        double max_rhs = 0.0;
        for ( double a : steepest_descent ) {
            if(fabs(a) > max_rhs) {
            	max_rhs = fabs(a);
            }
        }
        if(max_rhs > 1.0e-12) {
        	double inc = maxInc;
            if(inc > maxIncAllowed) {
            	inc = maxIncAllowed;
            }
        	steepest_descent.times(0.0001*inc/max_rhs);
        }

        FloatArray X_ref = X;
        FloatArray dX_ref = dX;
        FloatArray ddX_ref = ddX;

//        std::vector<double> scaleFactors = {0.1,0.3,0.7, 1.0, 1.5, 2.0};
        std::vector<double> scaleFactors = {0.5, 1.0, 1.5};
//        std::vector<double> residuals;

        double lowest_res = 1.0e20;
        double lowest_res_sf = 1.0;


        for( double sf : scaleFactors ) {

        	// Compute residual for the trial update
        	ddX = ddX_ref;
        	ddX.times(sf);

        	X = X_ref;
        	X.add(ddX);

        	dX = dX_ref;
        	dX.add(ddX);


            // Compute the residual
            engngModel->updateComponent(tStep, InternalRhs, domain);
            if (nite || iR == NULL) {
                rhs.beDifferenceOf(RT, F);
            } else {
                rhs = R;
                if (iR) {
                    rhs.add(*iR); // add initial guess
                }
            }
            if ( this->prescribedDofsFlag ) {
                this->applyConstraintsToLoadIncrement(nite, k, rhs, rlm, tStep);
            }

        	converged = this->checkConvergence(RT, F, rhs, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag, false);
        	double res = rhs.computeNorm();
            if ( engngModel->giveProblemScale() == macroScale ) {
            	printf("sf: %e res: %e\n", sf, res );
            }

        	if( res < lowest_res ) {
        		lowest_res_sf = sf;
        		lowest_res = res;
        	}
        }

        // Update with the best increment length found
    	ddX = ddX_ref;
    	ddX.times(lowest_res_sf);

    	X = X_ref;
    	X.add(ddX);

    	dX = dX_ref;
    	dX.add(ddX);


        if(max_rhs > 1.0e-12 && false) {

            FloatArray lowest_X = X;
            FloatArray lowest_dX = dX;
            FloatArray lowest_ddX = ddX;

            for( double sf : scaleFactors ) {

				/////////////////////////////////
				// Try also with steepest descent direction
				ddX = steepest_descent;
				ddX.times(sf);

				X = X_ref;
				X.add(steepest_descent);

				dX = dX_ref;
				dX.add(steepest_descent);

				// Compute the residual
				engngModel->updateComponent(tStep, InternalRhs, domain);
				if (nite || iR == NULL) {
					rhs.beDifferenceOf(RT, F);
				} else {
					rhs = R;
					if (iR) {
						rhs.add(*iR); // add initial guess
					}
				}
				if ( this->prescribedDofsFlag ) {
					this->applyConstraintsToLoadIncrement(nite, k, rhs, rlm, tStep);
				}

				converged = this->checkConvergence(RT, F, rhs, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag, false);
				double res = rhs.computeNorm();
	            if ( engngModel->giveProblemScale() == macroScale ) {
	            	printf("steepest descent res: %e\n", res );
	            }

				if( res > lowest_res ) {
					ddX = ddX_ref;
					ddX.times(lowest_res_sf);

					X = X_ref;
					X.add(ddX);

					dX = dX_ref;
					dX.add(ddX);


				}
				else {
					lowest_res = res;

					lowest_X = X;
		            lowest_dX = dX;
		            lowest_ddX = ddX;

				}
				//////////////////////////////


            }

            X = lowest_X;
            dX = lowest_dX;
            ddX = lowest_ddX;
        }


        tStep->incrementStateCounter(); // update solution state counter
        tStep->incrementSubStepNumber();

        engngModel->giveExportModuleManager()->doOutput(tStep, true);
    }

    // Modify Load vector to include "quasi reaction"
    if ( R0 ) {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R.at( prescribedEqs.at(i) ) = F.at( prescribedEqs.at(i) ) - R0->at( prescribedEqs.at(i) ) - R.at( prescribedEqs.at(i) );
        }
    } else {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R.at( prescribedEqs.at(i) ) = F.at( prescribedEqs.at(i) ) - R.at( prescribedEqs.at(i) );
        }
    }

    this->lastReactions.resize(numberOfPrescribedDofs);

#ifdef VERBOSE
    if ( numberOfPrescribedDofs ) {
        // print quasi reactions if direct displacement control used
        OOFEM_LOG_INFO("\n");
        OOFEM_LOG_INFO("NRSolver:     Quasi reaction table                                 \n");
        OOFEM_LOG_INFO("NRSolver:     Node            Dof             Displacement    Force\n");
        double reaction;
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            reaction = R.at( prescribedEqs.at(i) );
            if ( R0 ) {
                reaction += R0->at( prescribedEqs.at(i) );
            }
            lastReactions.at(i) = reaction;
            OOFEM_LOG_INFO("NRSolver:     %-15d %-15d %-+15.5e %-+15.5e\n", prescribedDofs.at(2 * i - 1), prescribedDofs.at(2 * i),
                           X.at( prescribedEqs.at(i) ), reaction);
        }
        OOFEM_LOG_INFO("\n");
    }
#endif

    return status;

#endif
}
#endif




#if 0
NM_Status
TrustRegionSolver :: solve(SparseMtrx &k, FloatArray &R, FloatArray *R0, FloatArray *iR,
                  FloatArray &X, FloatArray &dX, FloatArray &F,
                  const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                  int &nite, TimeStep *tStep)
{

    // residual, iteration increment of solution, total external force
    FloatArray rhs, ddX, RT;
    double RRT;
    int neq = X.giveSize();
    bool converged, errorOutOfRangeFlag;
    ParallelContext *parallel_context = engngModel->giveParallelContext( this->domain->giveNumber() );

    if ( engngModel->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO("NRSolver: Iteration");
        if ( rtolf.at(1) > 0.0 ) {
            OOFEM_LOG_INFO(" ForceError");
        }
        if ( rtold.at(1) > 0.0 ) {
            OOFEM_LOG_INFO(" DisplError");
        }
        OOFEM_LOG_INFO("\n----------------------------------------------------------------------------\n");
    }

    l = 1.0;

    NM_Status status = NM_None;
    this->giveLinearSolver();

    // compute total load R = R+R0
    RT = R;
    if ( R0 ) {
        RT.add(* R0);
    }

    RRT = parallel_context->localNorm(RT);
    RRT *= RRT;

    ddX.resize(neq);
    ddX.zero();

    // Fetch the matrix before evaluating internal forces.
    // This is intentional, since its a simple way to drastically increase convergence for nonlinear problems.
    // (This old tangent is just used)
    // This improves convergence for many nonlinear problems, but not all. It may actually
    // cause divergence for some nonlinear problems. Therefore a flag is used to determine if
    // the stiffness should be evaluated before the residual (default yes). /ES

    engngModel->updateComponent(tStep, NonLinearLhs, domain);
    if ( this->prescribedDofsFlag ) {
        if ( !prescribedEqsInitFlag ) {
            this->initPrescribedEqs();
        }
        applyConstraintsToStiffness(k);
    }

    // Old and new residual norm
    double forceErr = 0.0;
    double oldForceErr = 0.0;
    double refErr = 0.0;

    double coeff = 5.0e-1;
    nite = 0;
    for ( nite = 0; ; ++nite ) {
        // Compute the residual
        engngModel->updateComponent(tStep, InternalRhs, domain);
        if (nite || iR == NULL) {
            rhs.beDifferenceOf(RT, F);
        } else {
            rhs = R;
            if (iR) {
                rhs.add(*iR); // add initial guess
            }
        }
        if ( this->prescribedDofsFlag ) {
            this->applyConstraintsToLoadIncrement(nite, k, rhs, rlm, tStep);
        }

        // convergence check
        converged = this->checkConvergence(RT, F, rhs, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag);

//        forceErr = this->forceErrVec.computeNorm();
        forceErr = rhs.computeNorm();
//        printf("forceErr: %e\n" , forceErr);
        if( nite == 0 ) {
        	oldForceErr = forceErr;
        	refErr = forceErr;
        }

        if ( errorOutOfRangeFlag ) {
            status = NM_NoSuccess;
            OOFEM_WARNING("Divergence reached after %d iterations", nite);
            break;
        } else if ( converged && ( nite >= minIterations ) ) {
            status |= NM_Success;
            break;
        } else if ( nite >= nsmax ) {
            OOFEM_LOG_DEBUG("Maximum number of iterations reached\n");
            break;
        }

//#ifdef __PETSC_MODULE
#if 1
    PetscSparseMtrx *lhs = dynamic_cast< PetscSparseMtrx * >(&k);
    if ( lhs && true ) {
//    	printf("Fetched Petsc matrix.\n");

    	int N = k.giveNumberOfRows();

    	Vec petsc_mat_diag;
        VecCreate(PETSC_COMM_SELF, & petsc_mat_diag);
        VecSetType(petsc_mat_diag, VECSEQ);
        VecSetSizes(petsc_mat_diag, PETSC_DECIDE, N);

    	MatGetDiagonal(lhs->mtrx, petsc_mat_diag);


    	if ( engngModel->giveProblemScale() == macroScale ) {
    		printf("(1.0+coeff): %e\n", (1.0+coeff) );
    	}

    	VecScale(petsc_mat_diag, (1.0+coeff) );
//    	VecScale(petsc_mat_diag, 1.5);

    	MatDiagonalSet(lhs->mtrx, petsc_mat_diag, INSERT_VALUES);

        VecDestroy(& petsc_mat_diag);

//    	FloatArray diag_scale_vec(N);
//    	for(int i = 0; i < N; i++) {
//////    		k.at(i,i) *= coeff;
//    		petsc_mat_diag[i] *= coeff;
//    	}


//    	MatDiagonalSet(lhs->mtrx, petsc_mat_diag, INSERT_VALUES);

//        Vec petsc_diag_scale_vec;
//        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, diag_scale_vec.giveSize(), diag_scale_vec.givePointer(), & petsc_diag_scale_vec);
//
//        MatDiagonalScaleLocal(lhs->mtrx, petsc_diag_scale_vec);

    }
#endif


#if 0
    PetscSparseMtrx *lhs = dynamic_cast< PetscSparseMtrx * >(&k);
    if ( lhs ) {
    	printf("Fetched Petsc matrix.\n");

    	int N = k.giveNumberOfRows();

    	Vec petsc_mat_diag;
        VecCreate(PETSC_COMM_SELF, & petsc_mat_diag);
        VecSetType(petsc_mat_diag, VECSEQ);
        VecSetSizes(petsc_mat_diag, PETSC_DECIDE, N);

    	MatGetDiagonal(lhs->mtrx, petsc_mat_diag);


//    	VecScale(petsc_mat_diag, coeff);
//            	VecScale(petsc_mat_diag, 1.0 + coeff);
    	for(int i = 0; i < N; i++) {
    		VecSetValue(petsc_mat_diag, i, coeff, ADD_VALUES);
    	}
    	printf("coeff: %e\n", coeff );
    	coeff *= 0.9;

    	MatDiagonalSet(lhs->mtrx, petsc_mat_diag, INSERT_VALUES);

        VecDestroy(& petsc_mat_diag);

//    	FloatArray diag_scale_vec(N);
//    	for(int i = 0; i < N; i++) {
//////    		k.at(i,i) *= coeff;
//    		petsc_mat_diag[i] *= coeff;
//    	}


//    	MatDiagonalSet(lhs->mtrx, petsc_mat_diag, INSERT_VALUES);

//        Vec petsc_diag_scale_vec;
//        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, diag_scale_vec.giveSize(), diag_scale_vec.givePointer(), & petsc_diag_scale_vec);
//
//        MatDiagonalScaleLocal(lhs->mtrx, petsc_diag_scale_vec);

    }
#endif


        if ( nite > 0 || !mCalcStiffBeforeRes ) {
            if ( ( NR_Mode == nrsolverFullNRM ) || ( ( NR_Mode == nrsolverAccelNRM ) && ( nite % MANRMSteps == 0 ) ) ) {
                engngModel->updateComponent(tStep, NonLinearLhs, domain);
                applyConstraintsToStiffness(k);
            }
        }

        if ( ( nite == 0 ) && ( deltaL < 1.0 ) ) { // deltaL < 1 means no increment applied, only equilibrate current state
            rhs.zero();
            R.zero();
            ddX = rhs;
        } else {


//            if ( engngModel->giveProblemScale() == macroScale ) {
//            	k.writeToFile("k.txt");
//            }

            linSolver->solve(k, rhs, ddX);
        }

        //
        // update solution
        //
        if ( this->lsFlag && ( nite > 0 ) ) { // Why not nite == 0 ?
            // line search
            LineSearchNM :: LS_status LSstatus;
            double eta;
            this->giveLineSearchSolver()->solve(X, ddX, F, R, R0, prescribedEqs, 1.0, eta, LSstatus, tStep);
        } else if ( this->constrainedNRFlag && ( nite > this->constrainedNRminiter ) ) {
            ///@todo This doesn't check units, it is nonsense and must be corrected / Mikael
            if ( this->forceErrVec.computeSquaredNorm() > this->forceErrVecOld.computeSquaredNorm() ) {
                OOFEM_LOG_INFO("Constraining increment to be %e times full increment...\n", this->constrainedNRalpha);
                ddX.times(this->constrainedNRalpha);
            }
            //this->giveConstrainedNRSolver()->solve(X, & ddX, this->forceErrVec, this->forceErrVecOld, status, tStep);
        }


        /////////////////////////////////////////

        double maxInc = 0.0;
        for ( double inc : ddX ) {
            if(fabs(inc) > maxInc) {
                maxInc = fabs(inc);
            }
        }

        bool restricted = false;
        double reductionFactor = 1.0;
        if(maxInc > mTrustRegionSize) {
            if ( engngModel->giveProblemScale() == macroScale ) {
            	printf("Restricting increment to stay within trust region. mTrustRegionSize: %e\n", mTrustRegionSize);
            }

            reductionFactor = mTrustRegionSize/maxInc;
            maxInc = mTrustRegionSize;
        	ddX.times(reductionFactor);

        	restricted = true;
        }

        /////////////////////////////////////////


        X.add(ddX);
        dX.add(ddX);


        ////////////////////////////////////////////////////////////////
        // Compute residual of trial point ( f(xk + sk) in Conn's notation ).
        engngModel->updateComponent(tStep, InternalRhs, domain);
//        engngModel->updateComponent(tStep, NonLinearLhs, domain);

        if (nite || iR == NULL || true) {
            rhs.beDifferenceOf(RT, F);
        } else {
            rhs = R;
            if (iR) {
                rhs.add(*iR); // add initial guess
            }
        }
        if ( this->prescribedDofsFlag ) {
            this->applyConstraintsToLoadIncrement(nite, k, rhs, rlm, tStep);
        }
        // convergence check
        bool printTosScreen = false;
        converged = this->checkConvergence(RT, F, rhs, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag, printTosScreen);

//        forceErr = this->forceErrVec.computeNorm();
        forceErr = rhs.computeNorm();
    	if ( engngModel->giveProblemScale() == macroScale ) {
    		printf("forceErr: %e oldForceErr: %e refErr: %e\n", forceErr, oldForceErr, refErr);
    	}
        ////////////////////////////////////////////////////////////////


//        double rho_k = ( oldForceErr - forceErr )/( (.5)*fabs( ddX.dotProduct(rhs) ) );
//        double rho_k = ( oldForceErr - forceErr )/( (.5)*maxInc );
//        double rho_k = ( oldForceErr - forceErr )/( oldForceErr*reductionFactor );
//        double rho_k = ( oldForceErr - forceErr )/( oldForceErr*maxInc );
//        double rho_k = ( oldForceErr - forceErr )/( oldForceErr*reductionFactor );

        double rho_c = 0.5;
        double denom1 = oldForceErr;
        if( denom1 > 1.0e-12 ) {
        	rho_c = ( oldForceErr - forceErr )/( denom1 );
        }

    	if ( engngModel->giveProblemScale() == macroScale ) {
    		printf("oldForceErr - forceErr: %e denom1: %e\n" , oldForceErr - forceErr, denom1 );
    	}

//        double rho_k = ( oldForceErr - forceErr )/( maxInc );

        double rho_k = rho_c;

        double rho_h = 0.5;
        double denom2 = refErr;
        if( denom2 > 1.0e-12 ) {
        	rho_h = ( refErr - forceErr )/( denom2 );
        }
        if( rho_h > rho_k  ) {
        	rho_k = rho_h;
        }

    	if ( engngModel->giveProblemScale() == macroScale ) {
    		printf("rho_c: %e rho_h: %e rho_k: %e\n", rho_c, rho_h, rho_k );
    	}

#if 0
        double c1 = 1.0;
        bool expand_region = false;
        double forceErrAlt = 0.0;

        if( rho_k < 0.0 || true ) {
//        	printf("rho_k < 0.0\n");

        	// Try another increment
            ddX.times(c1);
            X.add(ddX);
            dX.add(ddX);

            engngModel->updateComponent(tStep, InternalRhs, domain);
    //        engngModel->updateComponent(tStep, NonLinearLhs, domain);

            if (nite || iR == NULL || true) {
                rhs.beDifferenceOf(RT, F);
            } else {
                rhs = R;
                if (iR) {
                    rhs.add(*iR); // add initial guess
                }
            }
            if ( this->prescribedDofsFlag ) {
                this->applyConstraintsToLoadIncrement(nite, k, rhs, rlm, tStep);
            }
            // convergence check
            bool printTosScreen = false;
            converged = this->checkConvergence(RT, F, rhs, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag, printTosScreen);

    //        forceErr = this->forceErrVec.computeNorm();
            forceErrAlt = rhs.computeNorm();

            printf("forceErrAlt: %e\n\n\n", forceErrAlt );

            if( forceErrAlt < forceErr && forceErrAlt < oldForceErr ){
            	printf("\nTrust-region should be expanded.\n");
            	expand_region  = true;
            }


            // Reset
            X.subtract(ddX);
            dX.subtract(ddX);
            ddX.times(1.0/c1);

        }
#endif

        if( rho_k >= mEta1  ) {

        	if ( engngModel->giveProblemScale() == macroScale ) {
        		printf("Accepting trial solution.\n");
        	}

            oldForceErr = forceErr;

            refErr *= 0.9;
            coeff *= 0.9;
        }
        else {

        	if ( engngModel->giveProblemScale() == macroScale ) {
        		printf("Keeping old solution.\n");
        	}

            X.subtract(ddX);
            dX.subtract(ddX);

        }


        ////////////////////////////////////////////////////////////////
        // Update trust-region size

#if 0
        if( expand_region && restricted ) {
			if ( (c1+1.0)*maxInc > mTrustRegionSize ) {
	        	mTrustRegionSize = (c1+1.0)*mTrustRegionSize;
			}

			printf("mTrustRegionSize: %e\n", mTrustRegionSize );

            oldForceErr = forceErrAlt;
            ddX.times(c1);
            X.add(c1, ddX);
            dX.add(c1, ddX);

        }
        else {
#endif
			if( rho_k >= mEta2  ) {

		    	if ( engngModel->giveProblemScale() == macroScale ) {
					printf("rho_k >= mEta2.\n");
					printf("Very successful update.\n");
		    	}

				// Parameter on p.782 in Conn et al.
				double alpha1 = 2.5;
	//        	double alpha1 = 1.5;

				if ( alpha1*maxInc > mTrustRegionSize ) {
	//        		mTrustRegionSize = alpha1*maxInc;
					mTrustRegionSize = alpha1*mTrustRegionSize;
				}

//				coeff *= 0.5;

		    	if ( engngModel->giveProblemScale() == macroScale ) {
		    		printf("mTrustRegionSize: %e\n", mTrustRegionSize );
		    	}

			}
			else {

				if( rho_k >= mEta1 && rho_k < mEta2  ) {

			    	if ( engngModel->giveProblemScale() == macroScale ) {
						printf("rho_k >= mEta1 && rho_k < mEta2.\n");

						printf("Successful update.\n");
						printf("Keeping trust-region size.\n");

						printf("mTrustRegionSize: %e\n", mTrustRegionSize );
			    	}
				}
				else {

			    	if ( engngModel->giveProblemScale() == macroScale ) {
						printf("rho_k < mEta1.\n");

						printf("Unsuccessful update.\n");
						printf("Contracting trust-region.\n");
			    	}

					// Parameter on p.782 in Conn et al.
					double alpha2 = 0.25;
	//            	double alpha2 = 0.75;
	//        		mTrustRegionSize = alpha2*maxInc;

					mTrustRegionSize = alpha2*mTrustRegionSize;

	//        		if( mTrustRegionSize < 1.0e-5 ) {
	//        			mTrustRegionSize = 1.0e-5;
	//        		}

//					coeff *= 2.0;


			    	if ( engngModel->giveProblemScale() == macroScale ) {
			    		printf("mTrustRegionSize: %e\n", mTrustRegionSize );
			    	}
				}

			}

#if 0
        }
#endif

        tStep->incrementStateCounter(); // update solution state counter
        tStep->incrementSubStepNumber();

        engngModel->giveExportModuleManager()->doOutput(tStep, true);

    }




    // Modify Load vector to include "quasi reaction"
    if ( R0 ) {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R.at( prescribedEqs.at(i) ) = F.at( prescribedEqs.at(i) ) - R0->at( prescribedEqs.at(i) ) - R.at( prescribedEqs.at(i) );
        }
    } else {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R.at( prescribedEqs.at(i) ) = F.at( prescribedEqs.at(i) ) - R.at( prescribedEqs.at(i) );
        }
    }

    this->lastReactions.resize(numberOfPrescribedDofs);

#ifdef VERBOSE
    if ( numberOfPrescribedDofs ) {
        // print quasi reactions if direct displacement control used
        OOFEM_LOG_INFO("\n");
        OOFEM_LOG_INFO("NRSolver:     Quasi reaction table                                 \n");
        OOFEM_LOG_INFO("NRSolver:     Node            Dof             Displacement    Force\n");
        double reaction;
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            reaction = R.at( prescribedEqs.at(i) );
            if ( R0 ) {
                reaction += R0->at( prescribedEqs.at(i) );
            }
            lastReactions.at(i) = reaction;
            OOFEM_LOG_INFO("NRSolver:     %-15d %-15d %-+15.5e %-+15.5e\n", prescribedDofs.at(2 * i - 1), prescribedDofs.at(2 * i),
                           X.at( prescribedEqs.at(i) ), reaction);
        }
        OOFEM_LOG_INFO("\n");
    }
#endif

    return status;
}
#endif

bool
TrustRegionSolver :: checkConvergence(FloatArray &RT, FloatArray &F, FloatArray &rhs,  FloatArray &ddX, FloatArray &X,
                             double RRT, const FloatArray &internalForcesEBENorm,
                             int nite, bool &errorOutOfRange, bool printToScreen)
{
//	printf("Entering TrustRegionSolver :: checkConvergence.\n");

    double forceErr, dispErr;
    FloatArray dg_forceErr, dg_dispErr, dg_totalLoadLevel, dg_totalDisp;
    bool answer;
    EModelDefaultEquationNumbering dn;
    ParallelContext *parallel_context = engngModel->giveParallelContext( this->domain->giveNumber() );

    /*
     * The force errors are (if possible) evaluated as relative errors.
     * If the norm of applied load vector is zero (one may load by temperature, etc)
     * then the norm of reaction forces is used in relative norm evaluation.
     *
     * Note: This is done only when all dofs are included (nccdg = 0). Not implemented if
     * multiple convergence criteria are used.
     *
     */

    answer = true;
    errorOutOfRange = false;

    // Store the errors associated with the dof groups
    if ( this->constrainedNRFlag || true) {
        this->forceErrVecOld = this->forceErrVec; // copy the old values
        this->forceErrVec.resize( internalForcesEBENorm.giveSize() );
        forceErrVec.zero();
    }

    if ( internalForcesEBENorm.giveSize() > 1 ) { // Special treatment when just one norm is given; No grouping
        int nccdg = this->domain->giveMaxDofID();
        // Keeps tracks of which dof IDs are actually in use;
        IntArray idsInUse(nccdg);
        idsInUse.zero();
        // zero error norms per group
        dg_forceErr.resize(nccdg);
        dg_forceErr.zero();
        dg_dispErr.resize(nccdg);
        dg_dispErr.zero();
        dg_totalLoadLevel.resize(nccdg);
        dg_totalLoadLevel.zero();
        dg_totalDisp.resize(nccdg);
        dg_totalDisp.zero();
        // loop over dof managers
        for ( auto &dofman : domain->giveDofManagers() ) {
            if ( !parallel_context->isLocal(dofman.get()) ) {
                continue;
            }

            // loop over individual dofs
            for ( Dof *dof: *dofman ) {
                if ( !dof->isPrimaryDof() ) {
                    continue;
                }
                int eq = dof->giveEquationNumber(dn);
                int dofid = dof->giveDofID();
                if ( !eq ) {
                    continue;
                }

                dg_forceErr.at(dofid) += rhs.at(eq) * rhs.at(eq);
                dg_dispErr.at(dofid) += ddX.at(eq) * ddX.at(eq);
                dg_totalLoadLevel.at(dofid) += RT.at(eq) * RT.at(eq);
                dg_totalDisp.at(dofid) += X.at(eq) * X.at(eq);
                idsInUse.at(dofid)++;
            } // end loop over DOFs
        } // end loop over dof managers

        // loop over elements and their DOFs
        for ( auto &elem : domain->giveElements() ) {
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

            // loop over element internal Dofs
            for ( int idofman = 1; idofman <= elem->giveNumberOfInternalDofManagers(); idofman++ ) {
                DofManager *dofman = elem->giveInternalDofManager(idofman);
                // loop over individual dofs
                for ( Dof *dof: *dofman ) {
                    if ( !dof->isPrimaryDof() ) {
                        continue;
                    }
                    int eq = dof->giveEquationNumber(dn);
                    int dofid = dof->giveDofID();

                    if ( !eq ) {
                        continue;
                    }

                    dg_forceErr.at(dofid) += rhs.at(eq) * rhs.at(eq);
                    dg_dispErr.at(dofid) += ddX.at(eq) * ddX.at(eq);
                    dg_totalLoadLevel.at(dofid) += RT.at(eq) * RT.at(eq);
                    dg_totalDisp.at(dofid) += X.at(eq) * X.at(eq);
                    idsInUse.at(dofid)++;
                } // end loop over DOFs
            } // end loop over element internal dofmans
        } // end loop over elements

        // loop over boundary conditions and their internal DOFs
        for ( auto &bc : domain->giveBcs() ) {
            // loop over element internal Dofs
            for ( int idofman = 1; idofman <= bc->giveNumberOfInternalDofManagers(); idofman++ ) {
                DofManager *dofman = bc->giveInternalDofManager(idofman);
                // loop over individual dofs
                for ( Dof *dof: *dofman ) {
                    if ( !dof->isPrimaryDof() ) {
                        continue;
                    }
                    int eq = dof->giveEquationNumber(dn);
                    int dofid = dof->giveDofID();

                    if ( !eq ) {
                        continue;
                    }

                    dg_forceErr.at(dofid) += rhs.at(eq) * rhs.at(eq);
                    dg_dispErr.at(dofid) += ddX.at(eq) * ddX.at(eq);
                    dg_totalLoadLevel.at(dofid) += RT.at(eq) * RT.at(eq);
                    dg_totalDisp.at(dofid) += X.at(eq) * X.at(eq);
                    idsInUse.at(dofid)++;
                } // end loop over DOFs
            } // end loop over element internal dofmans
        } // end loop over elements

        // exchange individual partition contributions (simultaneously for all groups)
        FloatArray collectiveErr(nccdg);
        parallel_context->accumulate(dg_forceErr,       collectiveErr);
        dg_forceErr       = collectiveErr;
        parallel_context->accumulate(dg_dispErr,        collectiveErr);
        dg_dispErr        = collectiveErr;
        parallel_context->accumulate(dg_totalLoadLevel, collectiveErr);
        dg_totalLoadLevel = collectiveErr;
        parallel_context->accumulate(dg_totalDisp,      collectiveErr);
        dg_totalDisp      = collectiveErr;

        if ( engngModel->giveProblemScale() == macroScale && printToScreen) {
            OOFEM_LOG_INFO("NRSolver: %-5d", nite);
        }

        int maxNumPrintouts = 6;
        int numPrintouts = 0;

        //bool zeroNorm = false;
        // loop over dof groups and check convergence individually
        for ( int dg = 1; dg <= nccdg; dg++ ) {

            bool zeroFNorm = false, zeroDNorm = false;
            // Skips the ones which aren't used in this problem (the residual will be zero for these anyway, but it is annoying to print them all)
            if ( !idsInUse.at(dg) ) {
                continue;
            }

        	numPrintouts++;

            if ( (engngModel->giveProblemScale() == macroScale && numPrintouts <= maxNumPrintouts) && printToScreen) {
                OOFEM_LOG_INFO( "  %s:", __DofIDItemToString( ( DofIDItem ) dg ).c_str() );
            }

            if ( rtolf.at(1) > 0.0 ) {
                //  compute a relative error norm
                if ( dg_forceScale.find(dg) != dg_forceScale.end() ) {
                    forceErr = sqrt( dg_forceErr.at(dg) / ( dg_totalLoadLevel.at(dg) + internalForcesEBENorm.at(dg) +
                        idsInUse.at(dg)*dg_forceScale[dg]*dg_forceScale[dg] ) );
                } else if ( ( dg_totalLoadLevel.at(dg) + internalForcesEBENorm.at(dg) ) >= nrsolver_ERROR_NORM_SMALL_NUM ) {
                    forceErr = sqrt( dg_forceErr.at(dg) / ( dg_totalLoadLevel.at(dg) + internalForcesEBENorm.at(dg) ) );
                } else {
                    // If both external forces and internal ebe norms are zero, then the residual must be zero.
                    //zeroNorm = true; // Warning about this afterwards.
                    zeroFNorm = true;
                    forceErr = sqrt( dg_forceErr.at(dg) );
                }

                if ( forceErr > rtolf.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) {
                    errorOutOfRange = true;
                }
                if ( forceErr > rtolf.at(1) ) {
                    answer = false;
                }

                if ( (engngModel->giveProblemScale() == macroScale  && numPrintouts <= maxNumPrintouts) && printToScreen ) {
                    OOFEM_LOG_INFO(zeroFNorm ? " *%.3e" : "  %.3e", forceErr);
                }

                // Store the errors from the current iteration
                if ( this->constrainedNRFlag || true) {
                    forceErrVec.at(dg) = forceErr;
                }
            }

            if ( rtold.at(1) > 0.0 ) {
                // compute displacement error
                if ( dg_totalDisp.at(dg) >  nrsolver_ERROR_NORM_SMALL_NUM ) {
                    dispErr = sqrt( dg_dispErr.at(dg) / dg_totalDisp.at(dg) );
                } else {
                    ///@todo This is almost always the case for displacement error. nrsolveR_ERROR_NORM_SMALL_NUM is no good.
                    //zeroNorm = true; // Warning about this afterwards.
                    //zeroDNorm = true;
                    dispErr = sqrt( dg_dispErr.at(dg) );
                }
                if ( dispErr  > rtold.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) {
                    errorOutOfRange = true;
                }
                if ( dispErr > rtold.at(1) ) {
                    answer = false;
                }

                if ( (engngModel->giveProblemScale() == macroScale  && numPrintouts <= maxNumPrintouts) && printToScreen ) {
                    OOFEM_LOG_INFO(zeroDNorm ? " *%.3e" : "  %.3e", dispErr);
                }
            }
        }


        if ( engngModel->giveProblemScale() == macroScale && printToScreen) {
            OOFEM_LOG_INFO("\n");
        }

        //if ( zeroNorm ) OOFEM_WARNING("Had to resort to absolute error measure (marked by *)");
    } else { // No dof grouping
        double dXX, dXdX;

        if ( engngModel->giveProblemScale() == macroScale && printToScreen ) {
            OOFEM_LOG_INFO("NRSolver:     %-15d", nite);
        } else {
//            OOFEM_LOG_INFO("  NRSolver:     %-15d", nite);
        }


        forceErr = parallel_context->localNorm(rhs);
        forceErr *= forceErr;
        dXX = parallel_context->localNorm(X);
        dXX *= dXX;                                       // Note: Solutions are always total global values (natural distribution makes little sense for the solution)
        dXdX = parallel_context->localNorm(ddX);
        dXdX *= dXdX;

        if ( rtolf.at(1) > 0.0 ) {
            // we compute a relative error norm
            if ( ( RRT + internalForcesEBENorm.at(1) ) > nrsolver_ERROR_NORM_SMALL_NUM ) {
                forceErr = sqrt( forceErr / ( RRT + internalForcesEBENorm.at(1) ) );
            } else {
                forceErr = sqrt(forceErr);   // absolute norm as last resort
            }
            if ( fabs(forceErr) > rtolf.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) {
                errorOutOfRange = true;
            }
            if ( fabs(forceErr) > rtolf.at(1) ) {
                answer = false;
            }

            if ( engngModel->giveProblemScale() == macroScale && printToScreen ) {
                OOFEM_LOG_INFO(" %-15e", forceErr);
            }

            if ( this->constrainedNRFlag ) {
                // store the errors from the current iteration for use in the next
                forceErrVec.at(1) = forceErr;
            }
        }

        if ( rtold.at(1) > 0.0 ) {
            // compute displacement error
            // err is relative displacement change
            if ( dXX > nrsolver_ERROR_NORM_SMALL_NUM ) {
                dispErr = sqrt(dXdX / dXX);
            } else {
                dispErr = sqrt(dXdX);
            }
            if ( fabs(dispErr)  > rtold.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) {
                errorOutOfRange = true;
            }
            if ( fabs(dispErr)  > rtold.at(1) ) {
                answer = false;
            }

            if ( engngModel->giveProblemScale() == macroScale && printToScreen ) {
                OOFEM_LOG_INFO(" %-15e", dispErr);
            }
        }

        if ( engngModel->giveProblemScale() == macroScale && printToScreen ) {
            OOFEM_LOG_INFO("\n");
        }
    } // end default case (all dofs contributing)

    return answer;
}


} /* namespace oofem */
