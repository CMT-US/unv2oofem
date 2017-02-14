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

#ifndef SRC_SM_MATERIALS_CONCRETE3VISCREG_H_
#define SRC_SM_MATERIALS_CONCRETE3VISCREG_H_

#include "ConcreteMaterials/concrete3.h"

///@name Input fields for Concrete3ViscReg
//@{
#define _IFT_Concrete3ViscReg_Name "concrete3viscreg"
#define _IFT_Concrete3ViscReg_RegCoeff "regcoeff"
//@}

namespace oofem {
/**
 * This class implements a Concrete3 material in a finite element problem,
 * with the addition of viscous regularization to improve convergence
 * in the Newton iterations.
 *
 */

class Concrete3ViscReg : public Concrete3 {
public:
	Concrete3ViscReg(int n, Domain * d);
	virtual ~Concrete3ViscReg();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveClassName() const { return "Concrete3ViscReg"; }
    virtual const char *giveInputRecordName() const { return _IFT_Concrete3ViscReg_Name; }

    virtual void giveMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode,
                                             GaussPoint *gp,
                                             TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

protected:
    double mRegCoeff;
};

} /* namespace oofem */

#endif /* SRC_SM_MATERIALS_CONCRETE3VISCREG_H_ */
