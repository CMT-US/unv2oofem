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

#ifndef PRESCRIBEDFIELDSGRADIENTSSUBBCS_H_
#define PRESCRIBEDFIELDSGRADIENTSSUBBCS_H_

#include "prescribedfieldsgradientshomogenization.h"
#include "prescribedgradienthomogenization.h"
#include "activebc.h"
#include "timestep.h"

#include <memory>

#define _IFT_PrescribedFieldsGradientsSubBCs_Name   "prescribedfieldsgradientssubbcs"
#define _IFT_PrescribedFieldsGradientsSubBCs_sub_bcs "subbcs"

namespace oofem {
class Node;
class Element;
/**
 * Imposes a combination of PrescribedFieldsGradients boundary conditions
 * @author Adam Sciegaj
 */
class OOFEM_EXPORT PrescribedFieldsGradientsSubBCs : public ActiveBoundaryCondition, public PrescribedFieldsGradientsHomogenization
{
public:
    PrescribedFieldsGradientsSubBCs(int n, Domain *d);
    virtual ~PrescribedFieldsGradientsSubBCs();

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    bcType giveType() const override { return UnknownBT; }

    const char *giveClassName() const override { return "PrescribedFieldsGradientsSubBCs"; }
    const char *giveInputRecordName() const override { return _IFT_PrescribedFieldsGradientsSubBCs_Name; }

    void computeTangent(FloatMatrix &tangent, TimeStep *tStep) override;

    void setPrescribedGradient(const FloatMatrix &t) override;
    void setPrescribedGradientVoigt(const FloatArray &e) override;
    void setPrescribedField(const FloatArray &s) override;
    void setPrescribedGradientVoigtAsymmetric(const FloatArray &g) override;

    void computeStress(FloatArray &sigma, TimeStep *tStep) override;
    void computeTransferStress(FloatArray &bStress, TimeStep *tStep) override;
    void computeReinfStress(FloatArray &rStress, TimeStep *tStep) override;

protected:
    // Array containing numbers of combined boundary conditions
    IntArray sub_bcs;
};
} /* namespace oofem */

#endif /* PRESCRIBEDFIELDSGRADIENTSSUBBCS_H_ */
