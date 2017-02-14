/*
 * concrete3viscreg.C
 *
 *  Created on: Feb 13, 2017
 *      Author: svennine
 */

#include "concrete3viscreg.h"

#include "ConcreteMaterials/concrete3.h"
#include "classfactory.h"


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
Concrete3ViscReg :: giveMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep)
{
//	printf("Entering Concrete3ViscReg :: giveMaterialStiffnessMatrix\n");

	Concrete3::giveMaterialStiffnessMatrix(answer, mode, gp, tStep);

	const double dt = tStep->giveTimeIncrement();

	for( int i = 0;  i < answer.giveNumberOfColumns(); i++ ) {
		answer(i,i) += (mRegCoeff/dt);
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
