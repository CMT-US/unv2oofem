/*
 * trustregionsolver4.C
 *
 *  Created on: Jun 20, 2017
 *      Author: svennine
 */

#include "trustregionsolver4.h"

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

#include "petscsparsemtrx.h"
#include <petscvec.h>

namespace oofem {

#define nrsolver_ERROR_NORM_SMALL_NUM 1.e-6
#define NRSOLVER_MAX_REL_ERROR_BOUND 1.e20
#define NRSOLVER_MAX_RESTARTS 4
#define NRSOLVER_RESET_STEP_REDUCE 0.25
#define NRSOLVER_DEFAULT_NRM_TICKS 10


REGISTER_SparseNonLinearSystemNM(TrustRegionSolver4)

TrustRegionSolver4::TrustRegionSolver4(Domain * d, EngngModel * m) :
NRSolver(d,m),
mTrustRegionSize(1.0e-3),
epsInit(false)
{



}

TrustRegionSolver4::~TrustRegionSolver4()
{
	if ( epsInit ) {
        EPSDestroy(& eps);
    }
}

IRResultType
TrustRegionSolver4 :: initializeFrom(InputRecord *ir) {

    IRResultType result;                // Required by IR_GIVE_FIELD macro


    IR_GIVE_OPTIONAL_FIELD(ir, mTrustRegionSize, _IFT_TrustRegionSolver4_InitialSize);

	if ( engngModel->giveProblemScale() == macroScale ) {
		printf("mTrustRegionSize: %e\n", mTrustRegionSize);
	}

    return NRSolver :: initializeFrom(ir);
}


NM_Status
TrustRegionSolver4 :: solve(SparseMtrx &k, FloatArray &R, FloatArray *R0,
                  FloatArray &X, FloatArray &dX, FloatArray &F,
                  const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                  int &nite, TimeStep *tStep)
{

    // residual, iteration increment of solution, total external force
    FloatArray rhs, ddX, RT;
    FloatArray rhs_old;
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


//    double initial_res = 0.0;
    double old_res = 0.0;
//    double sd_trial_res = 0.0;
    double newton_trial_res = 0.0;
    double eig_trial_res = 0.0;

//    bool first_perturbation = true;
    FloatArray eig_vec, pert_eig_vec;
//    double pert_tol = 1.0e2; // Worked well for mesh 1
    double pert_tol = 2.0e2;
//    double pert_tol = 5.0e2;

//    bool recompute_eig_vec = false;

    nite = 0;
    for ( nite = 0; ; ++nite ) {


        // Compute the residual
        engngModel->updateComponent(tStep, InternalRhs, domain);
        rhs.beDifferenceOf(RT, F);

        old_res = rhs.computeNorm();
        rhs_old = rhs;

        // convergence check
        converged = this->checkConvergence(RT, F, rhs, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag);


        if(nite == 0) {
//        	initial_res = old_res;
//	    	if ( engngModel->giveProblemScale() == macroScale ) {
//	    		printf("initial_res: %e\n", initial_res);
//	    	}
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

        engngModel->updateComponent(tStep, NonLinearLhs, domain);



    	////////////////////////////////////////////////////////////////////////////
        // Step calculation: Solve trust-region subproblem
        A = dynamic_cast< PetscSparseMtrx * >(&k);



        // Check if k is positive definite



        double smallest_eig_val = 0.0;

        // Dirty hack for weakly periodic boundary conditions
        PrescribedGradientBCWeak *bc = dynamic_cast<PrescribedGradientBCWeak*>(domain->giveBc(1));
        IntArray loc_u;
        if( bc ) {

	    	if ( engngModel->giveProblemScale() == macroScale ) {
	    		printf("Found PrescribedGradientBCWeak.\n");
	    	}

	    	auto B = bc->giveKuu(loc_u, tStep);
	        PetscSparseMtrx *Kuu;
	    	Kuu = dynamic_cast< PetscSparseMtrx * >( B.get() );
            calcSmallestEigVal(smallest_eig_val, eig_vec, *Kuu);

            if ( engngModel->giveProblemScale() == macroScale ) {
            	printf("smallest_eig_val: %e\n", smallest_eig_val);
            }

        }
        else {
            calcSmallestEigVal(smallest_eig_val, eig_vec, *A);
        }

        double lambda = 0.0;
        if(smallest_eig_val < pert_tol) {
        	lambda = 3.0*( pert_tol-1.0*smallest_eig_val );
        	addOnDiagonal(lambda, *A);
        }


        linSolver->solve(k, rhs, ddX);







    	if ( engngModel->giveProblemScale() == macroScale ) {
			printf("smallest_eig_val: %e\n", smallest_eig_val );
    	}


    	////////////////////////////////////////////////////////////////////////////
        // Check residual of Newton update
//        calcTrialRes(newton_trial_res, X, dX, ddX_newton, tStep, nite, iR, R, F, RT);


//        double maxIncPrel = giveMaxAbs(ddX);



    	////////////////////////////////////////////////////////////////////////////
        // Check residual of eigenpoint
    	FloatArray ddX_eig;
    	eig_trial_res = 1.0e20;
        if( smallest_eig_val < pert_tol ) {

        	printf("Negative eigenvalue detected.\n");

        	if(bc) {

        		pert_eig_vec.resize( ddX.giveSize() );

				for(int i = 0; i < loc_u.giveSize(); i++) {
					pert_eig_vec( loc_u(i)-1 ) = eig_vec(i);
				}
        	}
        	else {
        		pert_eig_vec = eig_vec;
        	}

        	// Rescale eigenvector such that the L_inf norm is 1.
        	double max_eig_vec = giveMaxAbs(pert_eig_vec);
        	printf("max_eig_vec: %e\n", max_eig_vec );
        	pert_eig_vec.times(1./max_eig_vec);


        	double c = mTrustRegionSize;

        	if( rhs_old.dotProduct(pert_eig_vec) > 0.0 ) {
        		c *= -1.0;
        	}

//        	if(maxIncPrel >= mTrustRegionSize) {
//
//        		ddX.zero();
//        		ddX.add( c, pert_eig_vec );
//        		ddX_eig = ddX;
////        		calcTrialRes(eig_trial_res, X, dX, ddX_eig, tStep, nite, iR, R, F, RT);
//
//
////        		X.add(ddX_eig);
////        		dX.add(ddX_eig);
//
//        	}
//        	else {

	//            ddX.zero();
				ddX.add( c*1.0e-2, pert_eig_vec );
	//        	ddX_eig = ddX;
	//        	calcTrialRes(eig_trial_res, X, dX, ddX_eig, tStep, nite, iR, R, F, RT);
	//
	//
	//            X.add(ddX_eig);
	//            dX.add(ddX_eig);
//        	}
        }


            double maxInc = giveMaxAbs(ddX);
            double increment_ratio_newton = 1.0;
            clipToLimit(ddX, mTrustRegionSize, increment_ratio_newton);
        	FloatArray ddX_newton = ddX;

            if ( engngModel->giveProblemScale() == macroScale ) {
            	if( increment_ratio_newton < 0.999 ) {
            		printf("Restricting increment. maxInc: %e increment ratio: %e\n", maxInc, increment_ratio_newton);
            	}
            }


            X.add(ddX_newton);
            dX.add(ddX_newton);

//            calcTrialRes(newton_trial_res, X, dX, ddX_newton, tStep, nite, iR, R, F, RT);

//        printf("old_res: %e newton_trial_res: %e eig_trial_res: %e sd_trial_res: %e\n", old_res, newton_trial_res, eig_trial_res, sd_trial_res );



        updateTrustRegionSize(old_res, newton_trial_res, eig_trial_res);

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

    return status;
}

void TrustRegionSolver4 :: updateTrustRegionSize(const double &iOldRes, const double &iNewtonTrialRes, const double &iEigTrialRes)
{

	// Keep the trus-region size fixed
	return;

#if 0
	////////////////////////////////////////////////////////////////////////////
    // Trust-region radius update
    if(smallest_eig_val > 0.0) {

    	double mEta1 = 0.0;
		double mEta2 = 1.0;

		if( rho_k >= mEta2  ) {

			if ( engngModel->giveProblemScale() == macroScale ) {
//				printf("rho_k >= mEta2.\n");
				printf("Very successful update.\n");
			}

			// Parameter on p.782 in Conn et al.
			double alpha1 = 1.5;

//				if(nite > mFixIterSize) {
				if ( alpha1*maxInc > mTrustRegionSize ) {
					mTrustRegionSize = alpha1*mTrustRegionSize;
					if ( engngModel->giveProblemScale() == macroScale ) {
						printf("Expanding trust-region size.\n");
					}
				}
//				}

			if ( engngModel->giveProblemScale() == macroScale ) {
				printf("mTrustRegionSize: %e\n", mTrustRegionSize );
			}

		}
		else {
			if( rho_k >= mEta1 && rho_k < mEta2  ) {

				if ( engngModel->giveProblemScale() == macroScale ) {
//					printf("rho_k >= mEta1 && rho_k < mEta2.\n");

					printf("Successful update.\n");
					printf("Keeping trust-region size.\n");
//
					printf("mTrustRegionSize: %e\n", mTrustRegionSize );
				}
			}
			else {

//				if(lambda < 1.0e-3 || true) {
				int mFixIterSize = 0;
				if(nite > mFixIterSize) {

					if ( engngModel->giveProblemScale() == macroScale ) {

						printf("Unsuccessful update.\n");
						printf("Contracting trust-region.\n");
					}

					// Parameter on p.782 in Conn et al.
					double alpha2 = 0.5;

					mTrustRegionSize = alpha2*mTrustRegionSize;
				}

				if ( engngModel->giveProblemScale() == macroScale ) {
					printf("mTrustRegionSize: %e\n", mTrustRegionSize );
				}
			}

		}
    }
#endif
}

double TrustRegionSolver4 :: giveMaxAbs(const FloatArray &iVec) const
{
	double max_eig_vec = 0.0;
	for ( double inc : iVec ) {
		if(fabs(inc) > max_eig_vec) {
			max_eig_vec = fabs(inc);
		}
	}

	return max_eig_vec;
}

void TrustRegionSolver4 :: clipToLimit(FloatArray &ioVec, const double &iLimit, double &oIncrementRatio)
{

    double maxInc = giveMaxAbs(ioVec);

    if(maxInc > iLimit) {
        ioVec.times(mTrustRegionSize/maxInc);
        oIncrementRatio = mTrustRegionSize/maxInc;

    }

}

void TrustRegionSolver4 :: calcTrialRes(double &oTrialRes, FloatArray &iX, FloatArray &idX, FloatArray &iddX_trial, TimeStep *tStep, int nite, FloatArray *iR, FloatArray &R, FloatArray &F, FloatArray &RT)
{

	FloatArray rhs;

	iX.add(1.0, iddX_trial);
    idX.add(1.0, iddX_trial);


	////////////////////////////////////////////////////////////////////////////
    // Check residual of Newton update
    engngModel->updateComponent(tStep, InternalRhs, domain);
    if (nite || iR == NULL) {
        rhs.beDifferenceOf(RT, F);
    } else {
        rhs = R;
        if (iR) {
            rhs.add(*iR); // add initial guess
        }
    }

    oTrialRes = rhs.computeNorm();


    iX.add(-1.0, iddX_trial);
    idX.add(-1.0, iddX_trial);

}


bool
TrustRegionSolver4 :: checkConvergence(FloatArray &RT, FloatArray &F, FloatArray &rhs,  FloatArray &ddX, FloatArray &X,
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

void TrustRegionSolver4::checkPetscError(PetscErrorCode iErrorCode) const {

	if( iErrorCode != 0 ) {
		printf("In TrustRegionSolver4::checkPetscError: iErrorCode %d\n", int(iErrorCode) );
	}

}

void TrustRegionSolver4::calcSmallestEigVal(double &oEigVal, FloatArray &oEigVec, PetscSparseMtrx &K) {
    PetscErrorCode ierr;
    ST st;

    double eig_rtol = 1.0e-4;
    int max_iter = 1000;
    int nroot = 10;
//    int size = K.giveNumberOfRows();

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
//        CHKERRQ(ierr);
        epsInit = true;
    }

    ierr = EPSSetOperators( eps, * K.giveMtrx(), NULL );
    checkPetscError(ierr);
//    CHKERRQ(ierr);

    ierr = EPSSetProblemType(eps, EPS_NHEP);
    checkPetscError(ierr);
//    CHKERRQ(ierr);

    ierr = EPSGetST(eps, & st);
    checkPetscError(ierr);
//    CHKERRQ(ierr);

//    ierr = STSetType(st, STSINVERT);
        ierr = STSetType(st, STSHIFT);
        checkPetscError(ierr);
//        ierr = STSetType(st, STCAYLEY);
//    CHKERRQ(ierr);
    ierr = STSetMatStructure(st, SAME_NONZERO_PATTERN);
    checkPetscError(ierr);
//    CHKERRQ(ierr);
    ierr = EPSSetTolerances(eps, ( PetscReal ) eig_rtol, max_iter);
    checkPetscError(ierr);
//    CHKERRQ(ierr);
    ierr = EPSSetDimensions(eps, ( PetscInt ) nroot, PETSC_DECIDE, PETSC_DECIDE);
    checkPetscError(ierr);
//    CHKERRQ(ierr);
    ierr = EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
    checkPetscError(ierr);
//    CHKERRQ(ierr);


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     *                   Solve the eigensystem
     *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    EPSConvergedReason eig_reason;
    int eig_nconv, eig_nite;

    ierr = EPSSolve(eps);
    checkPetscError(ierr);
//    CHKERRQ(ierr);

    ierr = EPSGetConvergedReason(eps, & eig_reason);
    checkPetscError(ierr);
//    CHKERRQ(ierr);
    ierr = EPSGetIterationNumber(eps, & eig_nite);
    checkPetscError(ierr);
//    CHKERRQ(ierr);
//    printf("SLEPcSolver::solve EPSConvergedReason: %d, number of iterations: %d\n", eig_reason, eig_nite);

    ierr = EPSGetConverged(eps, & eig_nconv);
    checkPetscError(ierr);
//    CHKERRQ(ierr);

    double smallest_eig_val = 1.0e20;

    if ( eig_nconv > 0 ) {
//        printf("SLEPcSolver :: solveYourselfAt: Convergence reached for RTOL=%20.15f\n", eig_rtol);

#if 1
        FloatArray eig_vals(nroot);
        PetscScalar kr;
        Vec Vr;

        K.createVecGlobal(& Vr);
//        ierr = MatGetVecs(* K.giveMtrx(), PETSC_NULL, & Vr);
//        checkPetscError(ierr);
//            CHKERRQ(ierr);

            FloatArray Vr_loc;

//        printf("\n\n");
        for ( int i = 0; i < eig_nconv && i < nroot; i++ ) {
        	// PetscErrorCode EPSGetEigenpair(EPS eps,PetscInt i,PetscScalar *eigr,PetscScalar *eigi,Vec Vr,Vec Vi)
//            ierr = EPSGetEigenpair(eps, eig_nconv - i - 1, & kr, PETSC_NULL, Vr, PETSC_NULL);
            ierr = EPSGetEigenpair(eps, i, & kr, PETSC_NULL, Vr, PETSC_NULL);
            checkPetscError(ierr);
//            CHKERRQ(ierr);

            //Store the eigenvalue
            eig_vals(i) = kr;

            if(kr < smallest_eig_val) {
            	smallest_eig_val = kr;

            	K.scatterG2L(Vr, Vr_loc);
            	oEigVec = Vr_loc;
            }

//            printf("i: %d ev: %e\n", i, kr);
        }

        //Store the eigenvector
//        for ( int j = 0; j < size; j++ ) {
//            _r.at(j + 1, i + 1) = Vr_loc.at(j + 1);
//        }

        ierr = VecDestroy(& Vr);
        checkPetscError(ierr);
//        printf("Vr_loc: "); Vr_loc.printYourself();

//        printf("\n\n");

#endif
    } else {
//        OOFEM_ERROR("No converged eigenpairs.\n");
    	printf("Warning: No converged eigenpairs.\n");
    }

    oEigVal = smallest_eig_val;

}


void TrustRegionSolver4::addOnDiagonal(const double &iVal, PetscSparseMtrx &K) {

#if 1

    	int N = K.giveNumberOfRows();

    	Vec petsc_mat_diag;
        VecCreate(PETSC_COMM_SELF, & petsc_mat_diag);
        VecSetType(petsc_mat_diag, VECSEQ);
        VecSetSizes(petsc_mat_diag, PETSC_DECIDE, N);

    	MatGetDiagonal(K.mtrx, petsc_mat_diag);


//    	VecScale(petsc_mat_diag, coeff);
//            	VecScale(petsc_mat_diag, 1.0 + coeff);
//    	double coeff = 0.0e3;
    	for(int i = 0; i < N; i++) {


    		// More dirty hacking...
    		// If the diagonal value is very close to zero, it is
    		// probably a Lagrange multiplier row. Don't add anything
    		// to such rows.

    		double a = 0.0;

    		VecGetValues(petsc_mat_diag, 1, &i, &a);

//    		if( fabs(a) > 1.0e-3 ) {
    			VecSetValue(petsc_mat_diag, i, iVal, ADD_VALUES);
//    		}
    	}
//    	printf("coeff: %e\n", coeff );

    	MatDiagonalSet(K.mtrx, petsc_mat_diag, INSERT_VALUES);

        VecDestroy(& petsc_mat_diag);

#endif


}


} /* namespace oofem */
