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

#include "prescribedfieldsgradientshomogenization.h"
#include "domain.h"
#include "dynamicinputrecord.h"
#include "set.h"
#include "feinterpol.h"
#include "element.h"
#include "mathfem.h"

namespace oofem {

IRResultType PrescribedFieldsGradientsHomogenization :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, mGradient1, _IFT_PrescribedFieldsGradientsHomogenization_gradient1);
    IR_GIVE_FIELD(ir, mGradient2, _IFT_PrescribedFieldsGradientsHomogenization_gradient2);
    IR_GIVE_FIELD(ir, mField1, _IFT_PrescribedFieldsGradientsHomogenization_field1);
    IR_GIVE_FIELD(ir, mField2, _IFT_PrescribedFieldsGradientsHomogenization_field2);

    mCenterCoord.resize( mGradient1.giveNumberOfColumns() );
    mCenterCoord.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, mCenterCoord, _IFT_PrescribedFieldsGradientsHomogenization_centercoords)

    return IRRT_OK;
}

void PrescribedFieldsGradientsHomogenization :: giveInputRecord(DynamicInputRecord &input)
{
    input.setField(mGradient1, _IFT_PrescribedFieldsGradientsHomogenization_gradient1);
    input.setField(mGradient2, _IFT_PrescribedFieldsGradientsHomogenization_gradient2);
    input.setField(mField1, _IFT_PrescribedFieldsGradientsHomogenization_field1);
    input.setField(mField2, _IFT_PrescribedFieldsGradientsHomogenization_field2);
    input.setField(mCenterCoord, _IFT_PrescribedFieldsGradientsHomogenization_centercoords);
}


void PrescribedFieldsGradientsHomogenization :: setPrescribedGradientVoigt(const FloatArray &t)
{
    int n = t.giveSize();
    if ( n == 3 ) { // Then 2D
        this->mGradient1.resize(2, 2);
        this->mGradient1.at(1, 1) = t.at(1);
        this->mGradient1.at(2, 2) = t.at(2);
        // In voigt form, assuming the use of gamma_12 instead of eps_12
        this->mGradient1.at(1, 2) = this->mGradient1.at(2, 1) = t.at(3) * 0.5;
    } else if ( n == 6 ) { // Then 3D
        this->mGradient1.resize(3, 3);
        this->mGradient1.at(1, 1) = t.at(1);
        this->mGradient1.at(2, 2) = t.at(2);
        this->mGradient1.at(3, 3) = t.at(3);
        // In voigt form, assuming the use of gamma_12 instead of eps_12
        this->mGradient1.at(1, 2) = this->mGradient1.at(2, 1) = t.at(6) * 0.5;
        this->mGradient1.at(1, 3) = this->mGradient1.at(3, 1) = t.at(5) * 0.5;
        this->mGradient1.at(2, 3) = this->mGradient1.at(3, 2) = t.at(4) * 0.5;
    } else if ( n == 1 ) {
        this->mGradient1.resize(1, 1);
        this->mGradient1.at(1, 1) = t.at(1);
    } else {
        OOFEM_ERROR("Tensor is in strange voigt format. Should be 3 or 6. Use setPrescribedGradient directly if needed.");
    }
}

void PrescribedFieldsGradientsHomogenization :: setPrescribedGradientVoigtAsymmetric(const FloatArray &t)
{
    int n = t.giveSize();
    if ( n == 4 ) { // Then 2D
        this->mGradient2.resize(2, 2);
        this->mGradient2.at(1, 1) = t.at(1);
        this->mGradient2.at(2, 2) = t.at(2);
        this->mGradient2.at(1, 2) = t.at(3);
        this->mGradient2.at(2, 1) = t.at(4);
    } else {
        OOFEM_ERROR("Tensor is in strange format. Should be 4.");
    }
}

void PrescribedFieldsGradientsHomogenization :: setPrescribedField(const FloatArray &s)
{
    int n = s.giveSize();
    if ( n == 2 ) { //then 2D)
        this->mField2.resize(2);
        this->mField2.at(1) = s.at(1);
        this->mField2.at(2) = s.at(2);
    } else {
        OOFEM_ERROR("Vector is in strange format. Should be 2.");
    }
}


void PrescribedFieldsGradientsHomogenization :: giveFirstGradient(FloatArray &oGradient) const
{
    int numRows = mGradient1.giveNumberOfRows();
    switch ( numRows ) {
    case 1:
        oGradient = FloatArray {
            mGradient1.at(1, 1)
        };
        break;
    case 2:
        // Do not assume symmetry
        oGradient = {
            mGradient1.at(1, 1), mGradient1.at(2, 2), mGradient1.at(1, 2), mGradient1.at(2, 1)
        };
        break;
    case 3:
        // TODO: Fix this properly.
        oGradient = {
            mGradient1.at(1, 1), mGradient1.at(2, 2), mGradient1.at(1, 2), mGradient1.at(2, 1)
        };
//        OOFEM_ERROR("PrescribedFieldsGradientsHomogenization :: giveGradientVoigt() not implemented for 3 rows.\n")
        break;
    }
}

void PrescribedFieldsGradientsHomogenization :: giveSecondGradient(FloatArray &oGradient) const
{
    int numRows = mGradient2.giveNumberOfRows();
    switch ( numRows ) {
    case 1:
        oGradient = FloatArray {
            mGradient2.at(1, 1)
        };
        break;
    case 2:
        // Do not assume symmetry
        oGradient = {
            mGradient2.at(1, 1), mGradient2.at(2, 2), mGradient2.at(1, 2), mGradient2.at(2, 1)
        };
        break;
    case 3:
        // TODO: Fix this properly.
        oGradient = {
            mGradient2.at(1, 1), mGradient2.at(2, 2), mGradient2.at(1, 2), mGradient2.at(2, 1)
        };
//        OOFEM_ERROR("PrescribedFieldsGradientsHomogenization :: giveGradientVoigt() not implemented for 3 rows.\n")
        break;
    }
}

void PrescribedFieldsGradientsHomogenization :: giveSecondField(FloatArray &sField) const
{
    sField = { mField2.at(1), mField2.at(2) };
}


double PrescribedFieldsGradientsHomogenization :: domainSize(Domain *d, int setNum)
{
    int nsd = d->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;
    // This requires the boundary to be consistent and ordered correctly.
    Set *set = d->giveSet(setNum);
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = d->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);
        FEInterpolation *fei = e->giveInterpolation();
        domain_size += fei->evalNXIntegral( boundary, FEIElementGeometryWrapper(e) );
    }
    return fabs(domain_size / nsd);
}
} /* namespace oofem */
