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

#include "prescribedfieldsgradientssubbcs.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "engngm.h"
#include "domain.h"
#include "node.h"

namespace oofem {
REGISTER_BoundaryCondition(PrescribedFieldsGradientsSubBCs);

PrescribedFieldsGradientsSubBCs :: PrescribedFieldsGradientsSubBCs(int n, Domain *d) :
    ActiveBoundaryCondition(n, d),
    PrescribedFieldsGradientsHomogenization()
{
}

PrescribedFieldsGradientsSubBCs :: ~PrescribedFieldsGradientsSubBCs()
{
}


IRResultType PrescribedFieldsGradientsSubBCs :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, sub_bcs, _IFT_PrescribedFieldsGradientsSubBCs_sub_bcs)
    ActiveBoundaryCondition :: initializeFrom(ir);
    return PrescribedFieldsGradientsHomogenization :: initializeFrom(ir);
}


void PrescribedFieldsGradientsSubBCs :: giveInputRecord(DynamicInputRecord &input)
{
    ActiveBoundaryCondition :: giveInputRecord(input);
    PrescribedFieldsGradientsHomogenization :: giveInputRecord(input);
    input.setField(this->sub_bcs, _IFT_PrescribedFieldsGradientsSubBCs_sub_bcs);
}


void PrescribedFieldsGradientsSubBCs :: setPrescribedGradient(const FloatMatrix &t)
{
    this->mGradient1 = t;
    for (int bc_num : this->sub_bcs) {
        auto bc = dynamic_cast< PrescribedFieldsGradientsHomogenization* >(this->giveDomain()->giveBc(bc_num));
        bc->setPrescribedGradient(t);
    }
}


void PrescribedFieldsGradientsSubBCs :: setPrescribedGradientVoigt(const FloatArray &e)
{
    PrescribedFieldsGradientsHomogenization :: setPrescribedGradientVoigt(e);
    for (int bc_num : this->sub_bcs) {
        auto bc = dynamic_cast<PrescribedFieldsGradientsHomogenization*>(this->giveDomain()->giveBc(bc_num));
        bc->setPrescribedGradientVoigt(e);
    }
}

void PrescribedFieldsGradientsSubBCs :: setPrescribedField(const FloatArray &s)
{
    this->mField2 = s;
    for (int bc_num : this->sub_bcs) {
        auto bc = dynamic_cast<PrescribedFieldsGradientsHomogenization*>(this->giveDomain()->giveBc(bc_num));
        bc->setPrescribedField(s);
    };
}

void PrescribedFieldsGradientsSubBCs :: setPrescribedGradientVoigtAsymmetric(const FloatArray &g)
{
    this->mGradient2 = g;
    for (int bc_num : this->sub_bcs) {
        auto bc = dynamic_cast<PrescribedFieldsGradientsHomogenization*>(this->giveDomain()->giveBc(bc_num));
        bc->setPrescribedGradientVoigtAsymmetric(g);
    }
}

void PrescribedFieldsGradientsSubBCs :: computeStress(FloatArray &sigma, TimeStep *tStep)
{
    FloatArray tempStress;
    sigma.clear();
    for (int bc_num : this->sub_bcs) {
        auto bc = dynamic_cast<PrescribedFieldsGradientsHomogenization*>(this->giveDomain()->giveBc(bc_num));
        bc->computeStress(tempStress, tStep);
        sigma.add(tempStress);
        tempStress.clear();
    }
}

void PrescribedFieldsGradientsSubBCs :: computeTransferStress(FloatArray &bStress, TimeStep *tStep)
{
    FloatArray tempbStress;
    bStress.clear();
    for (int bc_num : this->sub_bcs) {
        auto bc = dynamic_cast<PrescribedFieldsGradientsHomogenization*>(this->giveDomain()->giveBc(bc_num));
        bc->computeTransferStress(tempbStress, tStep);
        bStress.add(tempbStress);
        tempbStress.clear();
    }
}

void PrescribedFieldsGradientsSubBCs :: computeReinfStress(FloatArray &rStress, TimeStep *tStep)
{
    FloatArray temprStress;
    rStress.clear();
    for (int bc_num : this->sub_bcs) {
        auto bc = dynamic_cast<PrescribedFieldsGradientsHomogenization*>(this->giveDomain()->giveBc(bc_num));
        bc->computeReinfStress(temprStress, tStep);
        rStress.add(temprStress);
        temprStress.clear();
    }
}

void PrescribedFieldsGradientsSubBCs :: computeTangent(FloatMatrix &tangent, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented");
}

} /* namespace oofem */
