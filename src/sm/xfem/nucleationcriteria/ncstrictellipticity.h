/*
 * ncstrictellipticity.h
 *
 *  Created on: Mar 1, 2017
 *      Author: svennine
 */

#ifndef NCSTRICTELLIPTICITY_H_
#define NCSTRICTELLIPTICITY_H_

#define _IFT_NCStrictEllipticity_Name "ncstrictellipticity"
#define _IFT_NCStrictEllipticity_IncrementLength "incrementlength"
#define _IFT_NCStrictEllipticity_PropStrainThreshold "propagationstrainthreshold"
#define _IFT_NCStrictEllipticity_InitialCrackLength "initialcracklength"
#define _IFT_NCStrictEllipticity_CutOneEl "cut_one_el"
#define _IFT_NCStrictEllipticity_AllowInsertionInEnrEl "allow_ins_in_enr_el"

#include "xfem/nucleationcriterion.h"
#include <memory>

namespace oofem {

class FloatArray;

class NCStrictEllipticity : public NucleationCriterion {
public:
	NCStrictEllipticity(Domain *ipDomain);
	virtual ~NCStrictEllipticity();

	virtual std::vector<std::unique_ptr<EnrichmentItem>> nucleateEnrichmentItems();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void appendInputRecords(DynamicDataReader &oDR);

    /// @return Class name of the receiver.
    virtual const char *giveClassName() const {return "NCStrictEllipticity";}
    /// @return Input record name of the receiver.
    virtual const char *giveInputRecordName() const {return _IFT_NCStrictEllipticity_Name;};

protected:
    double mInitialCrackLength;
    double mIncrementLength;
    double mPropStrainThreshold;

    /// If the initiated crack should cut exactly one element.
    bool mCutOneEl;

    /// If nucleation of cracks in already enriched elements is allowed
    bool mAllowInsInEnrEl;

public:
    bool findLocalizationDirectionFromStrain(FloatArray &oN, const FloatArray &iEps);

};

} /* namespace oofem */

#endif /* NCSTRICTELLIPTICITY_H_ */
