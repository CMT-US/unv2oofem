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

#ifndef PRESCRIBEDGRADIENTSUBBCS_H_
#define PRESCRIBEDGRADIENTSUBBCS_H_

#include "prescribedgradienthomogenization.h"
#include "activebc.h"
#include "timestep.h"

#include <memory>

#define _IFT_PrescribedGradientSubBCs_Name   "prescribedgradientsubbcs"
#define _IFT_PrescribedGradientSubBCs_sub_bcs "subbcs"

namespace oofem {
class Node;
class Element;
/**
 * Imposes a combination of PrescribedGradient boundary conditions
 *
 */
class OOFEM_EXPORT PrescribedGradientSubBCs : public ActiveBoundaryCondition, public PrescribedGradientHomogenization
{
public:
    PrescribedGradientSubBCs(int n, Domain *d);
    virtual ~PrescribedGradientSubBCs();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);


    virtual bcType giveType() const { return UnknownBT; }

    virtual void scale(double s);

    virtual const char *giveClassName() const { return "PrescribedGradientSubBCs"; }
    virtual const char *giveInputRecordName() const { return _IFT_PrescribedGradientSubBCs_Name; }

    void computeField(FloatArray &sigma, TimeStep *tStep);
    void computeTangent(FloatMatrix &tangent, TimeStep *tStep);

    void setPrescribedGradient(const FloatMatrix &t);

    void setPrescribedGradientVoigt(const FloatArray &t);

    void setCenterCoordinate(FloatArray &x);

protected:
    /// DOF-manager containing the unknown homogenized stress.
    std :: unique_ptr< Node > mpSigmaHom;
    IntArray mSigmaIds;

    /// Help function that integrates the tangent contribution from a single element boundary.
    void integrateTangent(FloatMatrix &oTangent, Element *e, int iBndIndex);

    // Array containing numbers of combined boundary conditions
    IntArray sub_bcs;
};
} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBCNEUMANN_H_ */
