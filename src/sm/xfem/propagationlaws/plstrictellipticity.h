/*
 * plstrictellipticity.h
 *
 *  Created on: Apr 7, 2017
 *      Author: svennine
 */

#ifndef PLSTRICTELLIPTICITY_H_
#define PLSTRICTELLIPTICITY_H_

#include "xfem/propagationlaw.h"

#define _IFT_PLStrictEllipticity_Name "propagationlawstrictellipticity"
#define _IFT_PLStrictEllipticity_Radius "radius" ///< Radius away from tip used when picking sampling point
#define _IFT_PLStrictEllipticity_IncLength "incrementlength" ///< Increment length per time step

namespace oofem {
class Domain;
class EnrichmentDomain;
class DynamicInputRecord;

class PLStrictEllipticity : public PropagationLaw {
public:
	PLStrictEllipticity();
	virtual ~PLStrictEllipticity();

    virtual const char *giveClassName() const { return "PLStrictEllipticity"; }
    virtual const char *giveInputRecordName() const { return _IFT_PLStrictEllipticity_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual bool hasPropagation() const { return true; }
    virtual bool propagateInterface(Domain &iDomain, EnrichmentFront &iEnrFront, TipPropagation &oTipProp);

    void setRadius(double iRadius) {mRadius = iRadius;}
    void setIncrementLength(double iIncrementLength) {mIncrementLength = iIncrementLength;}

    static bool findLocalizationDirectionFromStrain(FloatArray &oN, const FloatArray &iEps);

protected:
    double mRadius, mIncrementLength;

};

} /* namespace oofem */

#endif /* PLSTRICTELLIPTICITY_H_ */
