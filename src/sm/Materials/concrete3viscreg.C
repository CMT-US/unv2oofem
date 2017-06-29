/*
 * concrete3viscreg.C
 *
 *  Created on: Feb 13, 2017
 *      Author: svennine
 */

#include "concrete3viscreg.h"

#include "ConcreteMaterials/concrete3.h"
#include "classfactory.h"

#include "gausspoint.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(Concrete3ViscReg);

Concrete3ViscReg::Concrete3ViscReg(int n, Domain *d) : Concrete3(n,d),
mRegCoeff(0.0)
{


}

Concrete3ViscReg::~Concrete3ViscReg() {

}

IRResultType
Concrete3ViscReg :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, mRegCoeff, _IFT_Concrete3ViscReg_RegCoeff);
//    printf("mRegCoeff: %e\n", mRegCoeff );

    return Concrete3 :: initializeFrom(ir);
}

void
Concrete3ViscReg :: giveInputRecord(DynamicInputRecord &input)
{
	Concrete3 :: giveInputRecord(input);

    input.setField(mRegCoeff, _IFT_Concrete3ViscReg_RegCoeff);
}

void
Concrete3ViscReg :: giveMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep)
{
//	printf("Entering Concrete3ViscReg :: giveMaterialStiffnessMatrix\n");

	bool useNumTangent = true;
    if ( useNumTangent ) {
//    	printf("Using numerical tangent.\n");
        // Numerical tangent
    	RCM2MaterialStatus *status = static_cast<RCM2MaterialStatus*>( this->giveStatus( gp ) );
        double h = 1.0e-6;

        const FloatArray &epsRed = status->giveTempStrainVector();
        FloatArray eps = epsRed;
//        StructuralMaterial::giveFullSymVectorForm(eps, epsRed, gp->giveMaterialMode() );


        int dim = eps.giveSize();
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

    } else {



		Concrete3::giveMaterialStiffnessMatrix(answer, mode, gp, tStep);

		const double dt = tStep->giveTimeIncrement();

		for( int i = 0;  i < answer.giveNumberOfColumns(); i++ ) {
			answer(i,i) += (mRegCoeff/dt);
		}

    }
}

void
Concrete3ViscReg :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                     const FloatArray &totalStrain,
                                     TimeStep *tStep)
{
//	printf("Entering Concrete3ViscReg :: giveRealStressVector\n");

	Concrete3::giveRealStressVector(answer, gp, totalStrain, tStep);


	RCM2MaterialStatus *ms = static_cast< RCM2MaterialStatus * >( this->giveStatus(gp) );

	const FloatArray &oldStrain = ms->giveStrainVector();
//	printf("oldStrain: "); oldStrain.printYourself();

	const FloatArray &newStrain = ms->giveTempStrainVector();
//	printf("newStrain: "); newStrain.printYourself();

	const double dt = tStep->giveTimeIncrement();
//	printf("dt: %e\n", dt);

	FloatArray strainRate;
	strainRate.beDifferenceOf(newStrain, oldStrain);
	strainRate.times(1./dt);

//	printf("strainRate: "); strainRate.printYourself();


	for( int i = 0;  i < answer.giveSize(); i++ ) {
		answer(i) += (mRegCoeff)*strainRate(i);
	}

}

} /* namespace oofem */
