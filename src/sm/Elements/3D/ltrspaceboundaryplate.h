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

#ifndef ltrspaceboundaryplate_h
#define ltrspaceboundaryplate_h

#include "sm/Elements/3D/ltrspaceboundary.h"

#define _IFT_LTRSpaceBoundaryPlate_Name "ltrspaceboundaryplate"
#define _IFT_LTRSpaceBoundaryPlate_Location "location"

namespace oofem {
class FEI3dTetLin;

/**
 * This class implements a linear tetrahedral four-node finite element.
 * Each node has 3 degrees of freedom. This element is used for 3D RVE analyses with Periodic Boundary Conditions.
 * At least one node is located at the image boundary.
 * These nodes are replaced with a periodic mirror nodes and a control node is used to impose the macroscopic (average) strain.
 * MACROSCOPIC INPUT: DEFORMATIONS AND CURVATURES (PLATE, 10 COMPONENTS: Hxx Hxy Hyx Hyy Hzx Hzy Gxx Gyy Gxy Gyx)
 * Hxx = du/dx, Hxy = du/dy, Hyx = dv/dx, Hyy = dv/dy (in-plane strains)
 * Hzx = dw/dx, Hzy = dw/dy (slopes)
 * Gxx = d^2(w)/dx^2 (Kirchhoff) or d(phi_x)/dx (Mindlin)  curvature
 * Gyy = d^2(w)/dy^2 (Kirchhoff) or d(phi_y)/dy (Mindlin)  curvature
 * Gxy = d^2(w)/dxdy (Kirchhoff) or d(phi_x)/dy (Mindlin)  curvature
 * Gyx = d^2(w)/dxdy (Kirchhoff) or d(phi_y)/dx (Mindlin)  curvature
 *
 * @author: Adam Sciegaj
 */
class LTRSpaceBoundaryPlate : public LTRSpaceBoundary
{
protected:
    void computeTransformationMatrix(FloatMatrix &answer, TimeStep *tStep) override;

public:
    LTRSpaceBoundaryPlate(int n, Domain * d);
    virtual ~LTRSpaceBoundaryPlate() { }

    int computeNumberOfDofs() override { return 22; };
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    // definition & identification
    IRResultType initializeFrom(InputRecord *ir) override;
    const char *giveInputRecordName() const override { return _IFT_LTRSpaceBoundaryPlate_Name; }
    const char *giveClassName() const override { return "LTRSpaceBoundaryPlate"; }
};
} // end namespace oofem
#endif // LTRSpaceBoundaryPlate_h
