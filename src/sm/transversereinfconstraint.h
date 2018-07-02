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
 *               Copyright (C) 1993 - 2018   Borek Patzak
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

#ifndef TRANSVERSEREINFCONSTRAINT_H_
#define TRANSVERSEREINFCONSTRAINT_H_

#include "activebc.h"

#include <memory>

#define _IFT_TransverseReinfConstraint_Name   "transversereinforcementconstraint"
#define _IFT_TransverseReinfConstraint_SteelElSet "steelelset"
#define _IFT_TransverseReinfConstraint_ConElBoundSet "conelboundset"

namespace oofem {
class Node;
class Element;
/**
 * Imposes transverse constraint on the reinforcement weakly
 *
 * @author Adam Sciegaj
 */
class OOFEM_EXPORT TransverseReinfConstraint : public ActiveBoundaryCondition
{
public:
    TransverseReinfConstraint(int n, Domain *d);
    virtual ~TransverseReinfConstraint();

    int giveNumberOfInternalDofManagers() override { return 1; }
    DofManager *giveInternalDofManager(int i) override;

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    bcType giveType() const override { return UnknownBT; }

    void assembleVector(FloatArray &answer, TimeStep *tStep,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorm = NULL) override;

    void assemble(SparseMtrx &answer, TimeStep *tStep,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale = 1.0) override;

    void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override;

    const char *giveClassName() const override { return "TransverseReinfConstraint"; }
    const char *giveInputRecordName() const override { return _IFT_TransverseReinfConstraint_Name; }

    void computeField(FloatArray &lambda, TimeStep *tStep);
    void computeTangent(FloatMatrix &tangent, TimeStep *tStep);

    void giveLambdaLocationArray(IntArray &oCols, const UnknownNumberingScheme &r_s);

protected:
    /// DOF-manager containing the unknown Lagrange multiplier
    std :: unique_ptr< Node > lmLambda;
    IntArray lmLambdaIds;

    /// Reinforcement bar element set
    int steelElSet;
    /// Element boundary set of the concrete element along the reinforcement
    int conElBoundSet;

    /// Help functions that integrate the tangent contribution from steel and concrete elements.
    void integrateTangent(FloatMatrix &oTangent, Element *es, Element *ec, int iBndIndex);
    void integrateTangentOnSteel(FloatMatrix &oTangent, Element *e);
    void integrateTangentOnConcrete(FloatMatrix &oTangent, Element *e, int iBndIndex);
};
} /* namespace oofem */

#endif /* TRANSVERSEREINFCONSTRAINT_H_ */
