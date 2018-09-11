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

#include "prescribedgradientsubbcs.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "engngm.h"
#include "domain.h"
#include "node.h"

namespace oofem {
REGISTER_BoundaryCondition(PrescribedGradientSubBCs);

PrescribedGradientSubBCs :: PrescribedGradientSubBCs(int n, Domain *d) :
    ActiveBoundaryCondition(n, d),
    PrescribedGradientHomogenization()
{
}

PrescribedGradientSubBCs :: ~PrescribedGradientSubBCs()
{
}


IRResultType PrescribedGradientSubBCs :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, sub_bcs, _IFT_PrescribedGradientSubBCs_sub_bcs)
    ActiveBoundaryCondition :: initializeFrom(ir);
    return PrescribedGradientHomogenization :: initializeFrom(ir);
}


void PrescribedGradientSubBCs :: giveInputRecord(DynamicInputRecord &input)
{
    ActiveBoundaryCondition :: giveInputRecord(input);
    PrescribedGradientHomogenization :: giveInputRecord(input);
    input.setField(this->sub_bcs, _IFT_PrescribedGradientSubBCs_sub_bcs);
}


void PrescribedGradientSubBCs :: setPrescribedGradient(const FloatMatrix &t)
{
    this->mGradient = t;
    for (int bc_num : this->sub_bcs) {
        auto bc = dynamic_cast< PrescribedGradientHomogenization* >(this->giveDomain()->giveBc(bc_num));
        bc->setPrescribedGradient(t);
    }
}


void PrescribedGradientSubBCs :: setPrescribedGradientVoigt(const FloatArray &t)
{
    PrescribedGradientHomogenization :: setPrescribedGradientVoigt(t);
    for (int bc_num : this->sub_bcs) {
        auto bc = dynamic_cast<PrescribedGradientHomogenization*>(this->giveDomain()->giveBc(bc_num));
        bc->setPrescribedGradientVoigt(t);
    }
}


void PrescribedGradientSubBCs :: setCenterCoordinate(FloatArray &x)
{
    this->mCenterCoord = x;
    for (int bc_num : this->sub_bcs) {
        auto bc = dynamic_cast<PrescribedGradientHomogenization*>(this->giveDomain()->giveBc(bc_num));
        bc->setCenterCoordinate(x);
    }
}


void PrescribedGradientSubBCs :: scale(double s)
{
    this->mGradient.times(s);

    for (int bc_num : this->sub_bcs) {
        auto bc = dynamic_cast<ActiveBoundaryCondition*>(this->giveDomain()->giveBc(bc_num));
        bc->scale(s);
    }
}

void PrescribedGradientSubBCs :: computeField(FloatArray &sigma, TimeStep *tStep)
{
    FloatArray tmp;
    sigma.clear();
    for (int bc_num : this->sub_bcs) {
        auto bc = dynamic_cast<PrescribedGradientHomogenization*>(this->giveDomain()->giveBc(bc_num));
        bc->computeField(tmp, tStep);
//         double size = bc->domainSize();
        sigma.add(tmp);
    }
//     double rve_size = this->domainSize();
//     sigma.times(1.0/rve_size);
}


void PrescribedGradientSubBCs :: computeTangent(FloatMatrix &tangent, TimeStep *tStep)
{
    OOFEM_ERROR("TODO");
}

} /* namespace oofem */
