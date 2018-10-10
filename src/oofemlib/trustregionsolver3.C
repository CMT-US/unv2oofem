/*
 * trustregionsolver3.C
 *
 *  Created on: Mar 20, 2017
 *      Author: svennine
 */
#if 1
#include "trustregionsolver3.h"

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

#include "prescribedgradientbcweak.h"

#include "petscsolver.h"
#include "petscsparsemtrx.h"
#include <petscvec.h>

namespace oofem {
#define nrsolver_ERROR_NORM_SMALL_NUM 1.e-6
#define NRSOLVER_MAX_REL_ERROR_BOUND 1.e20
#define NRSOLVER_MAX_RESTARTS 4
#define NRSOLVER_RESET_STEP_REDUCE 0.25
#define NRSOLVER_DEFAULT_NRM_TICKS 10


REGISTER_SparseNonLinearSystemNM(TrustRegionSolver3)


TrustRegionSolver3::TrustRegionSolver3(Domain * d, EngngModel * m) :
NRSolver(d,m),
mEta1(0.01),
mEta2(10.0),
mGamma1(0.5),
mGamma2(0.5),
mTrustRegionSize(1.0e-3),
mBeta(0.5),
mEigVecRecalc(100),
epsInit(false)
{

}

TrustRegionSolver3::~TrustRegionSolver3() {

	if ( epsInit ) {
        EPSDestroy(& eps);
    }
}


IRResultType
TrustRegionSolver3 :: initializeFrom(InputRecord *ir) {

    IRResultType result;                // Required by IR_GIVE_FIELD macro


    IR_GIVE_OPTIONAL_FIELD(ir, mTrustRegionSize, _IFT_TrustRegionSolver3_InitialSize);

	if ( engngModel->giveProblemScale() == macroScale ) {
		printf("mTrustRegionSize: %e\n", mTrustRegionSize);
	}


    IR_GIVE_OPTIONAL_FIELD(ir, mBeta, _IFT_TrustRegionSolver3_Beta);

	if ( engngModel->giveProblemScale() == macroScale ) {
		printf("mBeta: %e\n", mBeta);
	}


    IR_GIVE_OPTIONAL_FIELD(ir, mEigVecRecalc, _IFT_TrustRegionSolver3_EigVecRecompute);

	if ( engngModel->giveProblemScale() == macroScale ) {
		printf("mEigVecRecalc: %d\n", mEigVecRecalc);
	}

    return NRSolver :: initializeFrom(ir);

}



NM_Status
TrustRegionSolver3 :: solve(SparseMtrx &k, FloatArray &R, FloatArray *R0,
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


    double old_res = 0.0;
    double trial_res = 0.0;

    bool first_perturbation = true;
    FloatArray eig_vec, pert_eig_vec;
    double pert_tol = 0.0e1;


    nite = 0;
    for ( nite = 0; ; ++nite ) {

        // Compute the residual
        engngModel->updateComponent(tStep, InternalRhs, domain);
        rhs.beDifferenceOf(RT, F);

        old_res = rhs.computeNorm();

        // convergence check
        converged = this->checkConvergence(RT, F, rhs, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag);

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

        engngModel->updateComponent(tStep, NonLinearLhs, domain);



    	////////////////////////////////////////////////////////////////////////////
        // Step calculation: Solve trust-region subproblem
        PetscSparseMtrx &A = dynamic_cast< PetscSparseMtrx& >(k);
        IntArray loc_u;

        // Check if k is positive definite
        double smallest_eig_val = 0.0;

        // Dirty hack for weakly periodic boundary conditions
        PrescribedGradientBCWeak *bc = dynamic_cast<PrescribedGradientBCWeak*>(domain->giveBc(1));

        if( bc ) {

	    	if ( engngModel->giveProblemScale() == macroScale ) {
	    		printf("Found PrescribedGradientBCWeak.\n");
	    	}

	    	auto Kuu = std::dynamic_pointer_cast<PetscSparseMtrx>( bc->giveKuu(loc_u, tStep) );

            calcSmallestEigVal(smallest_eig_val, eig_vec, *Kuu);

            if ( engngModel->giveProblemScale() == macroScale ) {
            	printf("smallest_eig_val: %e\n", smallest_eig_val);
            }
        }
        else {
            calcSmallestEigVal(smallest_eig_val, eig_vec, A);
        }

        double lambda = 0.0;
        if(smallest_eig_val < 0.0) {
        	lambda = -smallest_eig_val;
        	A.addDiagonal(lambda, loc_u);
        }

        linSolver->solve(k, rhs, ddX);

        // Remove lambda from the diagonal again
        if(smallest_eig_val < 0.0) {
        	A.addDiagonal(-lambda, loc_u);
        }


        // Constrain the increment to stay within the trust-region
        double increment_ratio = 1.0;
        double maxInc = 0.0;

        for ( double inc : ddX ) {
            if(fabs(inc) > maxInc) {
                maxInc = fabs(inc);
            }
        }

        if(maxInc > mTrustRegionSize) {
            if ( engngModel->giveProblemScale() == macroScale ) {
            	printf("Restricting increment. maxInc: %e\n", maxInc);
            }
        	ddX.times(mTrustRegionSize/maxInc);
        	increment_ratio = mTrustRegionSize/maxInc;
        }



        if( smallest_eig_val < pert_tol ) {

            if ( engngModel->giveProblemScale() == macroScale ) {
				printf("Negative eigenvalue detected.\n");
				printf("Perturbing in lowest eigenvector direction.\n");
            }

        	if(first_perturbation || (nite%mEigVecRecalc==0) ) {
				pert_eig_vec.resize( ddX.giveSize() );

				for(int i = 0; i < loc_u.giveSize(); i++) {
					pert_eig_vec( loc_u(i)-1 ) = eig_vec(i);
				}


				// Rescale eigenvector such that the L_inf norm is 1.
				double max_eig_vec = 0.0;
				for ( double inc : pert_eig_vec ) {
					if(fabs(inc) > max_eig_vec) {
						max_eig_vec = fabs(inc);
					}
				}

				pert_eig_vec.times(1./max_eig_vec);

				first_perturbation = false;
        	}



        	double c = maxInc;
        	if(c > mTrustRegionSize) {
        		c = mTrustRegionSize;
        	}

        	if( ddX.dotProduct(pert_eig_vec) < 0.0 ) {
        		c *= -1.0;
        	}

        	ddX.add( c*mBeta, pert_eig_vec );

        }



    	if ( engngModel->giveProblemScale() == macroScale ) {
			printf("smallest_eig_val: %e increment_ratio: %e\n", smallest_eig_val, increment_ratio );
    	}


        X.add(ddX);
        dX.add(ddX);


    	////////////////////////////////////////////////////////////////////////////
        // Acceptance of trial point
        engngModel->updateComponent(tStep, InternalRhs, domain);
        rhs.beDifferenceOf(RT, F);

        trial_res = rhs.computeNorm();

        double rho_k = 1.0;

        if(old_res > 1.0e-12) {
        	rho_k = ( old_res - trial_res )/( 0.99*increment_ratio*old_res );
        }

    	////////////////////////////////////////////////////////////////////////////
        // Trust-region radius update
		if( rho_k >= mEta2  ) {

//				printf("Very successful update.\n");

			// Parameter on p.782 in Conn et al.
			double alpha1 = 2.5;

			if ( alpha1*maxInc > mTrustRegionSize ) {
				mTrustRegionSize = alpha1*mTrustRegionSize;
			}
		}
		else {

			if( !(rho_k >= mEta1 && rho_k < mEta2) ) {

				if(nite > 1000) { // Only contract trust-region size in case of emergency

					// Parameter on p.782 in Conn et al.
					double alpha2 = 0.5;

					mTrustRegionSize = alpha2*mTrustRegionSize;
				}

			}

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
}





bool
TrustRegionSolver3 :: checkConvergence(FloatArray &RT, FloatArray &F, FloatArray &rhs,  FloatArray &ddX, FloatArray &X,
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

void TrustRegionSolver3::checkPetscError(PetscErrorCode iErrorCode) const {

	if( iErrorCode != 0 ) {
		printf("In TrustRegionSolver3::checkPetscError: iErrorCode %d\n", int(iErrorCode) );
	}

}

void TrustRegionSolver3::calcSmallestEigVal(double &oEigVal, FloatArray &oEigVec, PetscSparseMtrx &K) {
    PetscErrorCode ierr;
    ST st;

    double eig_rtol = 1.0e-3;
    int max_iter = 10000;
    int nroot = 1;

    if ( !epsInit ) {
        /*
         * Create eigensolver context
         */
#ifdef __PARALLEL_MODE
        MPI_Comm comm = engngModel->giveParallelComm();
#else
        MPI_Comm comm = PETSC_COMM_SELF;
#endif
        ierr = EPSCreate(comm, & eps);
        checkPetscError(ierr);
        epsInit = true;
    }

    ierr = EPSSetOperators( eps, * K.giveMtrx(), NULL );
    checkPetscError(ierr);

    ierr = EPSSetProblemType(eps, EPS_NHEP);
    checkPetscError(ierr);

    ierr = EPSGetST(eps, & st);
    checkPetscError(ierr);

//        ierr = STSetType(st, STCAYLEY);
//        ierr = STSetType(st, STSINVERT);
        ierr = STSetType(st, STSHIFT);
        checkPetscError(ierr);

    ierr = STSetMatStructure(st, SAME_NONZERO_PATTERN);
    checkPetscError(ierr);

    ierr = EPSSetTolerances(eps, ( PetscReal ) eig_rtol, max_iter);
    checkPetscError(ierr);

    ierr = EPSSetDimensions(eps, ( PetscInt ) nroot, PETSC_DECIDE, PETSC_DECIDE);
    checkPetscError(ierr);

    ierr = EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
    checkPetscError(ierr);


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     *                   Solve the eigensystem
     *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    EPSConvergedReason eig_reason;
    int eig_nconv, eig_nite;

    ierr = EPSSolve(eps);
    checkPetscError(ierr);

    ierr = EPSGetConvergedReason(eps, & eig_reason);
    checkPetscError(ierr);

    ierr = EPSGetIterationNumber(eps, & eig_nite);
    checkPetscError(ierr);
//    printf("SLEPcSolver::solve EPSConvergedReason: %d, number of iterations: %d\n", eig_reason, eig_nite);

    ierr = EPSGetConverged(eps, & eig_nconv);
    checkPetscError(ierr);

    double smallest_eig_val = 1.0e20;

    if ( eig_nconv > 0 ) {
//        printf("SLEPcSolver :: solveYourselfAt: Convergence reached for RTOL=%20.15f\n", eig_rtol);

    	FloatArray eig_vals(nroot);
        PetscScalar kr;
        Vec Vr;

        K.createVecGlobal(& Vr);


            FloatArray Vr_loc;

        for ( int i = 0; i < eig_nconv && i < nroot; i++ ) {
        	// PetscErrorCode EPSGetEigenpair(EPS eps,PetscInt i,PetscScalar *eigr,PetscScalar *eigi,Vec Vr,Vec Vi)
            ierr = EPSGetEigenpair(eps, i, & kr, PETSC_NULL, Vr, PETSC_NULL);
            checkPetscError(ierr);

            //Store the eigenvalue
            eig_vals(i) = kr;

            if(kr < smallest_eig_val) {
            	smallest_eig_val = kr;

            	K.scatterG2L(Vr, Vr_loc);
            	oEigVec = Vr_loc;
            }

        }

        ierr = VecDestroy(& Vr);
        checkPetscError(ierr);

    } else {
//        OOFEM_ERROR("No converged eigenpairs.\n");
    	printf("Warning: No converged eigenpairs.\n");
    	smallest_eig_val = 1.0;
    }

    oEigVal = smallest_eig_val;
}


SparseLinearSystemNM *
TrustRegionSolver3 :: giveLinearSolver()
{
	solverType = ST_Petsc;

    if ( linSolver ) {
        if ( linSolver->giveLinSystSolverType() == solverType ) {
            return linSolver.get();
        }
    }

    linSolver = classFactory.createSparseLinSolver(solverType, domain, engngModel);
    if ( !linSolver ) {
        OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
    }

    return linSolver.get();
}

} /* namespace oofem */
#endif
