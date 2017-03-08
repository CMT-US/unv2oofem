/*
 * prescribedgradientbcweakshear.h
 *
 *  Created on: Mar 6, 2017
 *      Author: svennine
 */

#ifndef PRESCRIBEDGRADIENTBCWEAKSHEAR_H_
#define PRESCRIBEDGRADIENTBCWEAKSHEAR_H_

#include "prescribedgradientbcweakperiodic.h"

#define _IFT_PrescribedGradientBCWeakShear_Name   "prescribedgradientbcweakshear"

namespace oofem {

class PrescribedGradientBCWeakShear: public PrescribedGradientBCWeakPeriodic {
public:
	PrescribedGradientBCWeakShear(int n, Domain *d);
	virtual ~PrescribedGradientBCWeakShear();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void postInitialize();

    virtual const char *giveClassName() const { return "PrescribedGradientBCWeakShear"; }
    virtual const char *giveInputRecordName() const { return _IFT_PrescribedGradientBCWeakShear_Name; }

    virtual void assembleVector(FloatArray &answer, TimeStep *tStep,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorm = NULL);

    virtual void computeExtForceElContrib(FloatArray &oContrib, TracSegArray &iEl, int iDim, TimeStep *tStep);

protected:
    virtual void giveBoundaryCoordVector(FloatArray &oX, const FloatArray &iPos) const;

    virtual bool boundaryPointIsOnActiveBoundary(const FloatArray &iPos) const { return pointIsOnGammaPlus(iPos); }

    bool pointIsOnGammaPlus(const FloatArray &iPos) const;

    bool pointIsOnGammaR(const FloatArray &iPos) const;

    // Remove segments on inactive parts of the boundary
    virtual void removeSupressedSegments(TracSegArray &ioTSeg, const double &iAbsTol);

};

} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBCWEAKSHEAR_H_ */
