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

#ifndef PRESCRIBEDGRADIENTNN_H_
#define PRESCRIBEDGRADIENTNN_H_

#include "prescribedgradienthomogenization.h"
#include "activebc.h"

#include <memory>

#define _IFT_PrescribedGradientNN_Name   "prescribedgradientnn"
#define _IFT_PrescribedGradientNN_Thickness "thick"

namespace oofem {
class Node;
class Element;
/**
 * Prescribes a gradient as a Neumann-Neumann boundary condition, cf. Sciegaj et al. "Two-scale finite element
 * modelling of reinforced concrete structures: effective response and sub-scale fracture development. Works with 2D RVEs 
 * comprising solid elements (concrete) and reinforcement (beam elements). Useful in multi-scale analyses of reinforced concrete 
 * structures.
 *
 * Based on PrescribedGradientBCNeumann.
 * 
 * @author Adam Sciegaj
 */
class OOFEM_EXPORT PrescribedGradientNN : public ActiveBoundaryCondition, public PrescribedGradientHomogenization
{
public:
    PrescribedGradientNN(int n, Domain *d);
    virtual ~PrescribedGradientNN();

    int giveNumberOfInternalDofManagers() override { return 1; }
    DofManager *giveInternalDofManager(int i) override;

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    bcType giveType() const override { return UnknownBT; }

    void scale(double s) override;

    void assembleVector(FloatArray &answer, TimeStep *tStep,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorm = NULL) override;

    void assemble(SparseMtrx &answer, TimeStep *tStep,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale = 1.0) override;

    void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override;

    const char *giveClassName() const override { return "PrescribedGradientNN"; }
    const char *giveInputRecordName() const override { return _IFT_PrescribedGradientNN_Name; }

    void computeField(FloatArray &sigma, TimeStep *tStep) override;
    void computeTangent(FloatMatrix &tangent, TimeStep *tStep) override;

    void giveStressLocationArray(IntArray &oCols, const UnknownNumberingScheme &r_s);

protected:
    /// DOF-manager containing the unknown homogenized stress.
    std :: unique_ptr< Node > mpSigmaHom;
    IntArray mSigmaIds;

    /// Help function that integrates the tangent contribution from a single element boundary.
    void integrateTangent(FloatMatrix &oTangent, Element *e, int iBndIndex);
    
    /// Thickness of the 2D RVE
    double thick;
};
} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBCNEUMANN_H_ */
