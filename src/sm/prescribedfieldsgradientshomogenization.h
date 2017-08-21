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

#ifndef prescribedfieldsgradientshomogenization_h
#define prescribedfieldsgradientshomogenization_h

#include "inputrecord.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "dofmanager.h"

#include "error.h"

///@name Input fields for PrescribedFieldGradientsHomogenization
//@{
#define _IFT_PrescribedFieldsGradientsHomogenization_Name "prescribedfieldsgradients"
#define _IFT_PrescribedFieldsGradientsHomogenization_centercoords "ccoord"
#define _IFT_PrescribedFieldsGradientsHomogenization_field1 "field1"
#define _IFT_PrescribedFieldsGradientsHomogenization_gradient1 "gradient1"
#define _IFT_PrescribedFieldsGradientsHomogenization_field2 "field2"
#define _IFT_PrescribedFieldsGradientsHomogenization_gradient2 "gradient2"
//@}

namespace oofem {
class TimeStep;
class DynamicInputRecord;
class Domain;

/**
 * Class for homogenization of applied fields and their gradients.
 *
 * Extension to handle two different gradients and a field variable. Useful in multiscale analysis with
 * two macroscopic fields
 *
 * @author Mikael Öhman
 * @author Adam Sciegaj
 */
class OOFEM_EXPORT PrescribedFieldsGradientsHomogenization
{
protected:
    FloatMatrix mGradient1;
    FloatMatrix mGradient2;
    FloatArray mField1;
    FloatArray mField2;

    /// Center coordinate @f$ \bar{x}_i @f$
    FloatArray mCenterCoord;

    virtual double domainSize(Domain *d, int set);

public:
    PrescribedFieldsGradientsHomogenization() { }
    virtual ~PrescribedFieldsGradientsHomogenization() { }

    virtual double domainSize() { OOFEM_ERROR("Not implemented."); return 0.0; }

    /**
     * Initializes receiver according to object description stored in input record.
     * The input record contains two fields;
     * - gradient \#rows \#columns { d_11 d_12 ... ; d_21 ... } (required)
     * - cCoords \#columns x_1 x_2 ... (optional, default 0)
     * The prescribed gradients columns must be equal to the size of the center coordinates.
     * The size of the center coordinates must be equal to the size of the coordinates in the applied nodes.
     */
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual void computeStress(FloatArray &sigma, TimeStep *tStep) = 0;

    /**
     * Computes the homogenized, macroscopic effective bond stress
     * @param bStress Output quantity (bond stress * circumference / RVE area)
     * @param tStep Active time step.
     */
    virtual void computeBondStress(FloatArray &bStress, TimeStep *tStep) = 0;

    /**
     * Computes the homogenized, macroscopic effective reinforcement stress
     * (membrane stress acting on the reinforcement bars)
     * @param rStress Output quantity (membrane stress acting on reinforcement grid)
     * @param tStep Active time step.
     */
    virtual void computeReinfStress(FloatArray &rStress, TimeStep *tStep) = 0;

    /**
     * Computes the macroscopic tangent for homogenization problems through sensitivity analysis.
     * @param tangent Output tangent.
     * @param tStep Active time step.
     */
    virtual void computeTangent(FloatMatrix &tangent, TimeStep *tStep) = 0;

    /**
     * Set prescribed gradient.
     * @param t New prescribed gradient.
     */
    virtual void setPrescribedGradient(const FloatMatrix &t) { mGradient1 = t; }

    /**
     * Sets the prescribed gradient from the matrix from given voigt notation.
     * Assumes use of double values for off-diagonal, usually the way for strain in Voigt form.
     * @param t Vector in voigt format.
     */
    virtual void setPrescribedGradientVoigt(const FloatArray &t);
    /**
     * Gives back the applied gradient in Voigt form.
     * @param oGradient The applied gradient, in Voigt form.
     */
    void giveFirstGradient(FloatArray &oGradient) const;

    /**
     * Sets the prescribed gradient from the matrix from given voigt notation.
     * Symmetry is not assumed, i.e. off-diagonal values are not doubled.
     * @param t Vector in voigt format.
     */
    virtual void setPrescribedGradientVoigtAsymmetric(const FloatArray &t);
    /**
     * Gives back the applied gradient in Voigt form.
     * @param oGradient The applied gradient, in Voigt form.
     */
    void giveSecondGradient(FloatArray &oGradient) const;

    /**
     * Sets the precribed field from the vector
     * @param s Vector describing prescribed field at a location
     */
    virtual void setPrescribedField(const FloatArray &s);
    /**
     * Gives back the prescribed field vector.
     * @param sField The prescribed field
     */
    void giveSecondField(FloatArray &sField) const;

    /**
     * Set the center coordinate for the prescribed values to be set for.
     * @param x Center coordinate.
     */
    virtual void setCenterCoordinate(FloatArray &x) { mCenterCoord = x; }
    /// Returns the center coordinate
    FloatArray &giveCenterCoordinate() { return mCenterCoord; }

    virtual const char *giveClassName() const = 0;
};
} // end namespace oofem

#endif // prescribedfieldsgradientshomogenization_h
