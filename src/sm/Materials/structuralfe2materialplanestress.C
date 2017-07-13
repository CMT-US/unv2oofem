/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "structuralfe2materialplanestress.h"
#include "gausspoint.h"
#include "engngm.h"
#include "oofemtxtdatareader.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "util.h"
#include "contextioerr.h"
#include "generalboundarycondition.h"
#include "prescribedgradienthomogenization.h"
#include "exportmodulemanager.h"
#include "vtkxmlexportmodule.h"
#include "nummet.h"
#include "EngineeringModels/xfemsolverinterface.h"
#include "EngineeringModels/staticstructural.h"
#include "unknownnumberingscheme.h"
#include "xfem/xfemstructuremanager.h"
#include "mathfem.h"

#include "dynamicdatareader.h"

#include <sstream>

namespace oofem {
REGISTER_Material(StructuralFE2MaterialPlaneStress);

int StructuralFE2MaterialPlaneStress :: n = 1;

StructuralFE2MaterialPlaneStress :: StructuralFE2MaterialPlaneStress(int n, Domain *d) : StructuralMaterial(n, d),
useNumTangent(true)
{}

StructuralFE2MaterialPlaneStress :: ~StructuralFE2MaterialPlaneStress()
{}


IRResultType
StructuralFE2MaterialPlaneStress :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // Required by IR_GIVE_FIELD macro
    IR_GIVE_FIELD(ir, this->inputfile, _IFT_StructuralFE2MaterialPlaneStress_fileName);

    useNumTangent = ir->hasField(_IFT_StructuralFE2MaterialPlaneStress_useNumericalTangent);

    mRegCoeff = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, mRegCoeff, _IFT_StructuralFE2MaterialPlaneStress_RegCoeff);
    printf("mRegCoeff: %e\n", mRegCoeff );

    useExtStiff = ir->hasField(_IFT_StructuralFE2MaterialPlaneStress_useExternalStiffness);
    IR_GIVE_OPTIONAL_FIELD(ir, givenTangent, _IFT_StructuralFE2MaterialPlaneStress_useExternalStiffness);

    allGPRes = ir->hasField(_IFT_StructuralFE2MaterialPlaneStress_allGPResults);

    return StructuralMaterial :: initializeFrom(ir);
}


void
StructuralFE2MaterialPlaneStress :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);
    input.setField(this->inputfile, _IFT_StructuralFE2MaterialPlaneStress_fileName);

    if ( useNumTangent ) {
        input.setField(_IFT_StructuralFE2MaterialPlaneStress_useNumericalTangent);
    }

    if ( useExtStiff ) {
        input.setField(_IFT_StructuralFE2MaterialPlaneStress_useExternalStiffness);
        input.setField(givenTangent, _IFT_StructuralFE2MaterialPlaneStress_useExternalStiffness);
    }

    input.setField(mRegCoeff, _IFT_StructuralFE2MaterialPlaneStress_RegCoeff);
    input.setField(allGPRes, _IFT_StructuralFE2MaterialPlaneStress_allGPResults);
}


MaterialStatus *
StructuralFE2MaterialPlaneStress :: CreateStatus(GaussPoint *gp) const
{
    if ( allGPRes ) {
        int nel = gp->giveElement()->giveGlobalNumber();
        int gpn = gp->giveNumber();
        return new StructuralFE2MaterialPlaneStressStatus(gpn, nel, this->giveDomain(), gp, this->inputfile);
    } else {
        return new StructuralFE2MaterialPlaneStressStatus(1, 1, this->giveDomain(), gp, this->inputfile);
    }
}


void
StructuralFE2MaterialPlaneStress :: giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                 const FloatArray &totalStrain, TimeStep *tStep)
{
    FloatArray stress;
    StructuralFE2MaterialPlaneStressStatus *ms = static_cast< StructuralFE2MaterialPlaneStressStatus * >( this->giveStatus(gp) );

#if 0
	XfemStructureManager *xMan = dynamic_cast<XfemStructureManager*>( ms->giveRVE()->giveDomain(1)->giveXfemManager() );
	if(xMan) {
		printf("Total crack length in RVE: %e\n", xMan->computeTotalCrackLength() );
	}
#endif

    ms->setTimeStep(tStep);
    // Set input
    ms->giveBC()->setPrescribedGradientVoigt(totalStrain);
    // Solve subscale problem
    ms->giveRVE()->solveYourselfAt(tStep);
    // Post-process the stress
    ms->giveBC()->computeField(stress, tStep);

    if (stress.giveSize() == 4 ) {
        answer = {stress[0], stress[1], 0.5*(stress[2]+stress[3])};
    } else {
        StructuralMaterial::giveFullSymVectorForm(answer, stress, gp->giveMaterialMode() );
    }


    ////////////////////////////////////////////////////////
    // Viscous regularization

	const FloatArray &oldStrain = ms->giveStrainVector();

	FloatArray oldStrainFull = oldStrain;
	if(oldStrainFull.giveSize() < 6) {
		StructuralMaterial::giveFullSymVectorForm(oldStrainFull, oldStrain, gp->giveMaterialMode());
	}

	const FloatArray &newStrain = ms->giveTempStrainVector();
	FloatArray newStrainFull = newStrain;
	if(newStrainFull.giveSize() < 6) {
		StructuralMaterial::giveFullSymVectorForm(newStrainFull, newStrain, gp->giveMaterialMode());
	}

	const double dt = tStep->giveTimeIncrement();
//	printf("dt: %e\n", dt);

	FloatArray strainRate;
	strainRate.beDifferenceOf(newStrainFull, oldStrainFull);
	strainRate.times(1./dt);


	FloatArray strainRateFull = strainRate;
	if(strainRate.giveSize() < answer.giveSize()) {
		StructuralMaterial::giveFullSymVectorForm(strainRateFull, strainRate, gp->giveMaterialMode());
	}

	for( int i = 0;  i < answer.giveSize(); i++ ) {
		answer(i) += (mRegCoeff)*strainRateFull(i);
	}


    // Update the material status variables
    ms->letTempStressVectorBe(answer);
    ms->letTempStrainVectorBe(totalStrain);
    ms->markOldTangent(); // Mark this so that tangent is reevaluated if they are needed.
}

void StructuralFE2MaterialPlaneStress :: giveRealStressVector_StressControl(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, const IntArray &strainControl, TimeStep *tStep)
{
    this->giveRealStressVector_3d(answer, gp, reducedStrain, tStep);

    return;
}

void
StructuralFE2MaterialPlaneStress :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    if ( useNumTangent ) {
        // Numerical tangent
        StructuralFE2MaterialPlaneStressStatus *status = static_cast<StructuralFE2MaterialPlaneStressStatus*>( this->giveStatus( gp ) );
        double h = 1.0e-9;

        const FloatArray &epsRed = status->giveTempStrainVector();
        FloatArray eps;
        //StructuralMaterial::giveFullSymVectorForm(eps, epsRed, gp->giveMaterialMode() );


//         int dim = eps.giveSize();
        int dim = epsRed.giveSize();
        eps.resize(dim);
        answer.resize(dim, dim);
        answer.zero();

        FloatArray sig, sigPert, epsPert;

        for(int i = 1; i <= dim; i++) {
            // Add a small perturbation to the strain
            epsPert = eps;
            epsPert.at(i) += h;

            giveRealStressVector_3d(sigPert, gp, epsPert, tStep);
            answer.setColumn(sigPert, i);
        }

        giveRealStressVector_3d(sig, gp, eps, tStep);

        for(int i = 1; i <= dim; i++) {
            for(int j = 1; j <= dim; j++) {
                answer.at(j,i) -= sig.at(j);
                answer.at(j,i) /= h;
            }
        }

        status->setTangent(answer);

    } else if ( useExtStiff ) {

        StructuralFE2MaterialPlaneStressStatus *ms = static_cast< StructuralFE2MaterialPlaneStressStatus * >( this->giveStatus(gp) );
        answer = this->givenTangent;

        ms->setTangent(answer);

    } else {

        StructuralFE2MaterialPlaneStressStatus *ms = static_cast< StructuralFE2MaterialPlaneStressStatus * >( this->giveStatus(gp) );
        ms->computeTangent(tStep);
        const FloatMatrix &ans9 = ms->giveTangent();

        StructuralMaterial::giveReducedSymMatrixForm(answer, ans9, _3dMat);

//        const FloatMatrix &ans9 = ms->giveTangent();
//        printf("ans9: "); ans9.printYourself();
//
//        // Compute the (minor) symmetrized tangent:
//        answer.resize(6, 6);
//        for ( int i = 0; i < 6; ++i ) {
//            for ( int j = 0; j < 6; ++j ) {
//                answer(i, j) = ans9(i, j);
//            }
//        }
//        for ( int i = 0; i < 6; ++i ) {
//            for ( int j = 6; j < 9; ++j ) {
//                answer(i, j-3) += ans9(i, j);
//                answer(j-3, i) += ans9(j, i);
//            }
//        }
//        for ( int i = 6; i < 9; ++i ) {
//            for ( int j = 6; j < 9; ++j ) {
//                answer(j-3, i-3) += ans9(j, i);
//            }
//        }
//        for ( int i = 0; i < 6; ++i ) {
//            for ( int j = 3; j < 6; ++j ) {
//                answer(j, i) *= 0.5;
//                answer(i, j) *= 0.5;
//            }
//        }

    	const double dt = tStep->giveTimeIncrement();

    	for( int i = 0;  i < answer.giveNumberOfColumns(); i++ ) {
    		answer(i,i) += (mRegCoeff/dt);
    	}


#if 0
        // Numerical ATS for debugging
        FloatMatrix numericalATS(6, 6);
        FloatArray dsig;
        // Note! We need a copy of the temp strain, since the pertubations might change it.
        FloatArray tempStrain = ms->giveTempStrainVector();

        FloatArray sig, strain, sigPert;
        giveRealStressVector_3d(sig, gp, tempStrain, tStep);
        double hh = 1e-6;
        for ( int k = 1; k <= 6; ++k ) {
            strain = tempStrain;
            strain.at(k) += hh;
            giveRealStressVector_3d(sigPert, gp, strain, tStep);
            dsig.beDifferenceOf(sigPert, sig);
            numericalATS.setColumn(dsig, k);
        }
        numericalATS.times(1. / hh);
        giveRealStressVector_3d(sig, gp, tempStrain, tStep); // Reset

        //answer.printYourself("Analytical deviatoric tangent");
        //numericalATS.printYourself("Numerical deviatoric tangent");

        numericalATS.subtract(answer);
        double norm = numericalATS.computeFrobeniusNorm();
        if ( norm > answer.computeFrobeniusNorm() * 1e-3 && norm > 0.0 ) {
            OOFEM_ERROR("Error in deviatoric tangent");
        }
#endif
    }
}

void StructuralFE2MaterialPlaneStress :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    this->give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
}
//=============================================================================

StructuralFE2MaterialPlaneStressStatus :: StructuralFE2MaterialPlaneStressStatus(int n, int j, Domain * d, GaussPoint * g,  const std :: string & inputfile) :
StructuralMaterialStatus(n, d, g),
mNewlyInitialized(true)
{
	mInputFile = inputfile;

    this->oldTangent = true;

    if ( !this->createRVE(n, j, gp, inputfile) ) {
        OOFEM_ERROR("Couldn't create RVE");
    }

}

PrescribedGradientHomogenization* StructuralFE2MaterialPlaneStressStatus::giveBC()
{
	this->bc = dynamic_cast< PrescribedGradientHomogenization * >( this->rve->giveDomain(1)->giveBc(1) );
	return this->bc;
}


bool
StructuralFE2MaterialPlaneStressStatus :: createRVE(int n, int j, GaussPoint *gp, const std :: string &inputfile)
{
    OOFEMTXTDataReader dr( inputfile.c_str() );
    EngngModel *em = InstanciateProblem(& dr, _processor, 0); // Everything but nrsolver is updated.
    dr.finish();
    em->setProblemScale(microScale);
    em->checkProblemConsistency();
    em->initMetaStepAttributes( em->giveMetaStep(1) );
    em->giveNextStep(); // Makes sure there is a timestep (which we will modify before solving a step)
    em->init();

    this->rve.reset( em );

    std :: ostringstream name;
    name << this->rve->giveOutputBaseFileName() << "-gp" << n << "-el" << j;
    if ( this->domain->giveEngngModel()->isParallel() && this->domain->giveEngngModel()->giveNumberOfProcesses() > 1 ) {
        name << "." << this->domain->giveEngngModel()->giveRank();
    }

    this->rve->letOutputBaseFileNameBe( name.str() );

    this->bc = dynamic_cast< PrescribedGradientHomogenization * >( this->rve->giveDomain(1)->giveBc(1) );
    if ( !this->bc ) {
        OOFEM_ERROR("RVE doesn't have necessary boundary condition; should have a type of PrescribedGradientHomogenization as first b.c.");
    }

    return true;
}

void
StructuralFE2MaterialPlaneStressStatus :: setTimeStep(TimeStep *tStep)
{
    TimeStep *rveTStep = this->rve->giveCurrentStep(); // Should i create a new one if it is empty?
    rveTStep->setNumber( tStep->giveNumber() );
    rveTStep->setTime( tStep->giveTargetTime() );
    rveTStep->setTimeIncrement( tStep->giveTimeIncrement() );
}

void
StructuralFE2MaterialPlaneStressStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
}

void
StructuralFE2MaterialPlaneStressStatus :: markOldTangent() { this->oldTangent = true; }

void
StructuralFE2MaterialPlaneStressStatus :: computeTangent(TimeStep *tStep)
{
    if ( !tStep->isTheCurrentTimeStep() ) {
        OOFEM_ERROR("Only current timestep supported.");
    }

    if ( this->oldTangent ) {
        bc->computeTangent(this->giveTangent(), tStep);
    }

    this->oldTangent = false;
}

void 
StructuralFE2MaterialPlaneStressStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    this->rve->updateYourself(tStep);
    this->rve->terminate(tStep);

    mNewlyInitialized = false;
}


contextIOResultType
StructuralFE2MaterialPlaneStressStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return this->rve->saveContext(&stream, mode, obj);
}


contextIOResultType
StructuralFE2MaterialPlaneStressStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return this->rve->restoreContext(&stream, mode, obj);
}

double StructuralFE2MaterialPlaneStressStatus :: giveRveLength()
{
	double rveLength = sqrt( bc->domainSize() );
	return rveLength;
}

void StructuralFE2MaterialPlaneStressStatus :: copyStateVariables(const MaterialStatus &iStatus)
{
//	static int num = 0;
//	printf("Entering StructuralFE2MaterialPlaneStressStatus :: copyStateVariables.\n");

    this->oldTangent = true;

//    if ( !this->createRVE(this->giveNumber(), gp, mInputFile) ) {
//        OOFEM_ERROR("Couldn't create RVE");
//    }


	StructuralMaterialStatus::copyStateVariables(iStatus);


	//////////////////////////////
	MaterialStatus &tmpStat = const_cast< MaterialStatus & >(iStatus);
	StructuralFE2MaterialPlaneStressStatus *fe2ms = dynamic_cast<StructuralFE2MaterialPlaneStressStatus*>(&tmpStat);

	if(!fe2ms) {
		OOFEM_ERROR("Failed to cast StructuralFE2MaterialPlaneStressStatus.")
	}


	this->mNewlyInitialized = fe2ms->mNewlyInitialized;

	// The proper way to do this would be to clone the RVE from iStatus.
	// However, this is a mess due to all pointers that need to be tracked.
	// Therefore, we consider a simplified version: copy only the enrichment items.

	Domain *ext_domain = fe2ms->giveRVE()->giveDomain(1);
	Domain *rve_domain = rve->giveDomain(1);

	if( ext_domain->hasXfemManager() ) {


		XfemManager *ext_xMan = ext_domain->giveXfemManager();
		XfemManager *this_xMan = rve->giveDomain(1)->giveXfemManager();
		DynamicDataReader dataReader;
		if ( ext_xMan != NULL ) {

		    IRResultType result; // Required by IR_GIVE_FIELD macro
			std::vector<std::unique_ptr<EnrichmentItem>> eiList;

//			DynamicInputRecord *xmanRec = new DynamicInputRecord();
//			ext_xMan->giveInputRecord(* xmanRec);
//			dataReader.insertInputRecord(DataReader :: IR_xfemManRec, xmanRec);

			// Enrichment items
			int nEI = ext_xMan->giveNumberOfEnrichmentItems();
			for ( int i = 1; i <= nEI; i++ ) {
				EnrichmentItem *ext_ei = ext_xMan->giveEnrichmentItem(i);
				ext_ei->appendInputRecords(dataReader);


		        InputRecord *mir = dataReader.giveInputRecord(DataReader :: IR_enrichItemRec, i);
		        std :: string name;
		        result = mir->giveRecordKeywordField(name);

		        if ( result != IRRT_OK ) {
		            mir->report_error(this->giveClassName(), __func__, "", result, __FILE__, __LINE__);
		        }

		        std :: unique_ptr< EnrichmentItem >ei( classFactory.createEnrichmentItem( name.c_str(), i, this_xMan, rve_domain ) );
		        if ( ei.get() == NULL ) {
		            OOFEM_ERROR( "unknown enrichment item (%s)", name.c_str() );
		        }

		        ei->initializeFrom(mir);
		        ei->instanciateYourself(&dataReader);
		        eiList.push_back( std :: move(ei) );

			}

			this_xMan->clearEnrichmentItems();
			this_xMan->appendEnrichmentItems(eiList);

			rve_domain->postInitialize();
			rve->forceEquationNumbering();
		}


	}

	TimeStep *tStep = rve->giveCurrentStep();

	// Update state in bulk GPs
	for( auto &el : rve_domain->giveElements() ) {
//		printf("Mapping bulk state variables in fe2ms.\n");
		el->mapStateVariables(*ext_domain, *tStep);
	}


	// Update state in CZ elements

//	printf("done.\n");

#if 0
	Domain *newDomain = fe2ms->giveRVE()->giveDomain(1)->Clone();
	newDomain->SetEngngModel(rve.get());
	bool deallocateOld = true;
	rve->setDomain(1, newDomain, deallocateOld);

//	rve->giveDomain(1)->postInitialize();
	rve->giveNumericalMethod(NULL)->setDomain(newDomain);

	rve->postInitialize();
//	rve->forceEquationNumbering();

	rve->initMetaStepAttributes( rve->giveMetaStep(1) );
    rve->giveNextStep(); // Makes sure there is a timestep (which we will modify before solving a step)
    rve->init();


//    std :: ostringstream name;
//    name << this->rve->giveOutputBaseFileName() << "-gp" << n;
//    this->rve->letOutputBaseFileNameBe( name.str() );
//    n++;

    double crackLength = 0.0;
    XfemStructureManager *xMan = dynamic_cast<XfemStructureManager*>( rve->giveDomain(1)->giveXfemManager() );
    if(xMan) {
    	crackLength = xMan->computeTotalCrackLength();
    }

    std :: ostringstream name;
    name << this->rve->giveOutputBaseFileName() << "-gp" << num << "crackLength" << crackLength;
    if ( this->domain->giveEngngModel()->isParallel() && this->domain->giveEngngModel()->giveNumberOfProcesses() > 1 ) {
        name << "." << this->domain->giveEngngModel()->giveRank();
    }

    num++;

    this->rve->letOutputBaseFileNameBe( name.str() );


	// Update BC
	this->bc = dynamic_cast< PrescribedGradientHomogenization * >( this->rve->giveDomain(1)->giveBc(1) );

#if 1

    XfemSolverInterface *xfemSolInt = dynamic_cast<XfemSolverInterface*>(rve.get());
    StaticStructural *statStruct = dynamic_cast<StaticStructural*>(rve.get());
    if(xfemSolInt && statStruct) {
//    	printf("Successfully casted to XfemSolverInterface.\n");

    	TimeStep *tStep = rve->giveCurrentStep();

		EModelDefaultEquationNumbering num;
		int numDofsNew = rve->giveNumberOfDomainEquations( 1, num );
		FloatArray u;
		u.resize(numDofsNew);
		u.zero();

		xfemSolInt->xfemUpdatePrimaryField(*statStruct, tStep, u);

	    // Set domain pointer to various components ...
		rve->giveNumericalMethod(NULL)->setDomain(newDomain);
	//        ioEngngModel.nMethod->setDomain(domain);

    }

//    TimeStep *tStep = rve->giveNextStep();
//    setTimeStep(tStep);
//    rve->solveYourselfAt(tStep);



    int numExpModules = rve->giveExportModuleManager()->giveNumberOfModules();
    for ( int i = 1; i <= numExpModules; i++ ) {
        //  ... by diving deep into the hierarchies ... :-/
        VTKXMLExportModule *vtkxmlMod = dynamic_cast< VTKXMLExportModule * >( rve->giveExportModuleManager()->giveModule(i) );
        if ( vtkxmlMod != NULL ) {
            vtkxmlMod->giveSmoother()->setDomain(newDomain);
            vtkxmlMod->givePrimVarSmoother()->setDomain(newDomain);
        }
    }
#endif
#endif
}


} // end namespace oofem
