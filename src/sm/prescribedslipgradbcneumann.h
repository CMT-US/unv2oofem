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

#ifndef PRESCRIBEDSLIPGRADBCNEUMANN_H_
#define PRESCRIBEDSLIPGRADBCNEUMANN_H_

#include "activebc.h"
#include "prescribedfieldsgradientshomogenization.h"

#include <memory>

#define _IFT_PrescribedSlipGradientBCNeumann_Name   "prescribedslipgradientbcneumann"
#define _IFT_PrescribedSlipGradientBCNeumann_SteelElSet "steelelset"
#define _IFT_PrescribedSlipGradientBCNeumann_RebarSets "rebarsets"
#define _IFT_PrescribedSlipGradientBCNeumann_Version "version"

namespace oofem {
class Node;
class Element;
class CrossSection;

/**
 * Imposes effective slip gradient weakly on the RVE
 *
 * @author Adam Sciegaj
 */
class OOFEM_EXPORT PrescribedSlipGradientBCNeumann : public ActiveBoundaryCondition, public PrescribedFieldsGradientsHomogenization 
{
public:
    PrescribedSlipGradientBCNeumann(int n, Domain *d);
    virtual ~PrescribedSlipGradientBCNeumann() { };

    int giveNumberOfInternalDofManagers() override { return 1; }
    DofManager *giveInternalDofManager(int i) override;

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    bcType giveType() const override { return UnknownBT; }

    void scale(double s) override;

    void assembleVector(FloatArray &answer, TimeStep *tStep,
                        CharType type, ValueModeType mode,
                        const UnknownNumberingScheme &s, FloatArray *eNorm=nullptr) override;

    void assemble(SparseMtrx &answer, TimeStep *tStep,
                  CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale = 1.0) override;

    void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                            const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override;

    void giveReinfStressLocationArray(IntArray &oCols, const UnknownNumberingScheme &r_s);

    const char *giveClassName() const override { return "PrescribedSlipGradientBCNeumann"; }
    const char *giveInputRecordName() const override { return _IFT_PrescribedSlipGradientBCNeumann_Name; }

    /**
     * Computes the homogenized, macroscopic, field (stress).
     * @param sigma Output quantity (typically stress).
     * @param tStep Active time step.
     */
    void computeStress(FloatArray &stress, TimeStep *tStep) override { };

    /**
     * Computes the homogenized, macroscopic effective transfer stress
     * @param bStress Output quantity (bond stress * circumference / RVE area)
     * @param tStep Active time step.
     */
    void computeTransferStress(FloatArray &bStress, TimeStep *tStep) override { };

    /**
     * Computes the homogenized, macroscopic effective reinforcement stress
     * (membrane stress acting on the reinforcement bars)
     * @param rStress Output quantityt (membrane stress acting on reinforcement grid)
     * @param tStep Active time step.
     */
    void computeReinfStress(FloatArray &rStress, TimeStep *tStep) override;

    void computeTangent(FloatMatrix &tangent, TimeStep *tStep) override;

protected:
    /// DOF-manager containing the unknown Lagrange multiplier (homogenised reinforcement stress)
    std :: unique_ptr< Node > lmSigmaSHom;
    IntArray lmSigmaSIds;

    /// IntArray containing sets of individual reinforcement bars
    IntArray rebarSets;

    ///Variants of the BC : shear components present (default), shear components omitted
    enum slipGradVariant { sg_shear, sg_noShear } modelVersion;

    /// Help functions that integrate the tangent contribution
    void integrateTangent(FloatMatrix &oTangent, Element *es, Element *ec, int iBndIndex);
    void integrateTangentOnSteel(FloatMatrix &oTangent, Element *e, const int rebarSet);
    void integrateTangentOnConcrete(FloatMatrix& oTangent, Element *e, int iBndIndex);
    void computeWeightMatrix(FloatMatrix &C, const IntArray reinfSets);
    void computeRebarDyad(FloatMatrix &dyad, const int reinfSet);
    double computeInterfaceLength(const IntArray reinfSets);
};
} /* namespace oofem */

#endif /* PRESCRIBEDSLIPGRADBCNEUMANN_H_ */

