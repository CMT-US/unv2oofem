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

/*
 * prescribedshearbcneumann.h
 *
 *  Created on: Mar 7, 2017
 *      Author: svennine
 */

#ifndef PRESCRIBEDSHEARBCNEUMANN_H_
#define PRESCRIBEDSHEARBCNEUMANN_H_

#include "activebc.h"
#include "floatarray.h"

#include <memory>

#define _IFT_PrescribedShearBCNeumann_Name   "prescribedshearbcneumann"
#define _IFT_PrescribedShearBCNeumann_Shear   "shearstrain"

namespace oofem {

class Node;
class Element;

/**
 * Imposes a prescribed transverse shear in y-direction weakly on
 * the boundary with a Neumann boundary condition.
 *
 * @author Erik Svenning
 */
class OOFEM_EXPORT PrescribedShearBCNeumann : public ActiveBoundaryCondition {
public:
	PrescribedShearBCNeumann(int n, Domain *d);
	virtual ~PrescribedShearBCNeumann();

    virtual int giveNumberOfInternalDofManagers() { return 1; }
    virtual DofManager *giveInternalDofManager(int i);

    virtual const char *giveClassName() const { return "PrescribedShearBCNeumann"; }
    virtual const char *giveInputRecordName() const { return _IFT_PrescribedShearBCNeumann_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual bcType giveType() const { printf("Entering bcType giveType()\n"); return UnknownBT; }


    virtual void assembleVector(FloatArray &answer, TimeStep *tStep,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorm = NULL);

    virtual void assemble(SparseMtrx &answer, TimeStep *tStep,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    double giveL() const {return mUC[0] - mLC[0];}
    double giveH() const {return mUC[1] - mLC[1];}

    bool pointIsOnRightEdge(const FloatArray &iX, double iTol = 1.0e-12) const {return iX[0] > mUC[0]-iTol;}
    bool pointIsOnLeftEdge(const FloatArray &iX, double iTol = 1.0e-12) const {return iX[0] < mLC[0]+iTol;}

protected:

    /// Help function that integrates the tangent contribution from a single element.
    void integrateTangent(FloatMatrix &oTangent, Element *e, int iBndIndex);

    /// Prescribed shear strain
    double mShearStrain;

    /// DOF-manager for the unknown shear stress
    std :: unique_ptr< Node > mpShearStress;
    IntArray mShearStressDofIds;

    /// Lower and upper corner of the domain
    FloatArray mLC, mUC;
};

} /* namespace oofem */

#endif /* PRESCRIBEDSHEARBCNEUMANN_H_ */
