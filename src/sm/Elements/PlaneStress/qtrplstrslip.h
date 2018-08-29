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

#ifndef qtrplstrdlip_h
#define qtrplstrslip_h

#include "Elements/structural2delement.h"
#include "ErrorEstimators/directerrorindicatorrc.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

#define _IFT_QTrPlaneStress2dSlip_Name "qtrplstrslip"

namespace oofem {
class FEI2dTrQuad;

/**
 * This class implements a quadratic triangular 6-node plane-
 * stress elasticity finite element with independent slip field. Each node has 2 degrees of freedom.
 * Works only in FE2 setting, with structuralfe2materialplanestress.
 */
class QTrPlaneStress2dSlip : public PlaneStressElement, public SpatialLocalizerInterface,
public SPRNodalRecoveryModelInterface
{
protected:
    static FEI2dTrQuad interpolation;

public:
    QTrPlaneStress2dSlip(int n, Domain * d);
    virtual ~QTrPlaneStress2dSlip() { }

    Interface *giveInterface(InterfaceType it) override;
    FEInterpolation *giveInterpolation() const override;
    double giveParentElSize() const override { return 0.5; }

    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;

    void computeHomogenizedFields(FloatArray &Stress, FloatArray &bStress, FloatArray &rStress, const FloatArray &strain, const FloatArray &slip, const FloatArray &slipGradient, GaussPoint *gp, TimeStep *tStep);

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_QTrPlaneStress2dSlip_Name; }
    const char *giveClassName() const override { return "QTrPlaneStress2dSlip"; }
    IRResultType initializeFrom(InputRecord *ir) override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override;
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

protected:
    int giveNumberOfIPForMassMtrxIntegration() override { return 4; }

};
} // end namespace oofem
#endif // qtrplstrslip_h
