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

#ifndef prescribedslipgradientsdd_h
#define prescribedslipgradientsdd_h

#include "prescribedfieldsgradientshomogenization.h"
#include "boundarycondition.h"
#include "dof.h"
#include "bctype.h"
#include "valuemodetype.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "dofmanager.h"

///@name Input fields for PrescribedGradientDD
//@{
#define _IFT_PrescribedSlipGradientsDD_Name "prescribedslipgradientsdd"
#define _IFT_PrescribedSlipGradientsDD_ConcreteBoundary "conboundset"
#define _IFT_PrescribedSlipGradientsDD_Thickness "thick"
#define _IFT_PrescribedSlipGradientsDD_ReinfXBound "reinfxbound"
#define _IFT_PrescribedSlipGradientsDD_ReinfYBound "reinfybound"

//@}

namespace oofem {
/**
 * Prescribes a gradient as a Dirichlet-Dirichlet boundary condition, cf. Sciegaj et al. "Two-scale finite element
 * modelling of reinforced concrete structures: effective response and sub-scale fracture development. Works with 2D RVEs
 * comprising solid elements (concrete) and reinforcement (beam elements). Useful in multi-scale analyses of reinforced concrete
 * structures. Support only for orthogonal reinforcement in X and Y direction.
 *
 * @author Adam Sciegaj
 */
class OOFEM_EXPORT PrescribedSlipGradientsDD : public BoundaryCondition, public PrescribedFieldsGradientsHomogenization
{
public:
    /**
     * Creates boundary condition with given number, belonging to given domain.
     * @param n Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    PrescribedSlipGradientsDD(int n, Domain *d) : BoundaryCondition(n, d) { }

    /// Destructor
    virtual ~PrescribedSlipGradientsDD() { }

    virtual double give(Dof *dof, ValueModeType mode, double time);

    virtual bcType giveType() const { return DirichletBT; }

    /**
     * Initializes receiver according to object description stored in input record.
     * The input record contains two fields;
     * - gradient \#rows \#columns { d_11 d_12 ... ; d_21 ... } (required)
     * - cCoords \#columns x_1 x_2 ... (optional, default 0)
     * The prescribed tensor's columns must be equal to the size of the center coordinates.
     * The size of the center coordinates must be equal to the size of the coordinates in the applied nodes.
     */
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    /**
     * Constructs a coefficient matrix for all prescribed unknowns.
     * Helper routine for computational homogenization.
     * @todo Perhaps this routine should only give results for the dof it prescribes?
     * @param C Coefficient matrix to fill.
     */
    void updateCoefficientMatrix(FloatMatrix &C);

    /**
     * Computes the homogenized, macroscopic, field (stress).
     * @param sigma Output quantity (typically stress).
     * @param tStep Active time step.
     */
    virtual void computeStress(FloatArray &stress, TimeStep *tStep);

    /**
     * Computes the homogenized, macroscopic effective bond stress
     * @param bStress Output quantity (bond stress * circumference / RVE area)
     * @param tStep Active time step.
     */
    virtual void computeTransferStress(FloatArray &bStress, TimeStep *tStep);

    /**
     * Computes the homogenized, macroscopic effective reinforcement stress
     * (membrane stress acting on the reinforcement bars)
     * @param rStress Output quantityt (membrane stress acting on reinforcement grid)
     * @param tStep Active time step.
     */
    virtual void computeReinfStress(FloatArray &rStress, TimeStep *tStep);

    /**
     * Computes the macroscopic tangent for homogenization problems through sensitivity analysis.
     * @param tangent Output tangent.
     * @param tStep Active time step.
     */
    virtual void computeTangent(FloatMatrix &tangent, TimeStep *tStep);

    virtual void scale(double s) { mGradient1.times(s); }

    virtual const char *giveClassName() const { return "PrescribedSlipGradientsDD"; }
    virtual const char *giveInputRecordName() const { return _IFT_PrescribedSlipGradientsDD_Name; }

protected:
    /**
     * Element boundaries set of the concrete solid
     */
    int conboundset;
    /**
     * Thickness of the RVE
     */
    double thick;
    /**
     * Node set containing end nodes of rebars in x-dir (horizontal).
     */
    int reinfxbound;
    /**
     * Node set contaning end nodes of rebars in y-dir (vertical).
     */
    int reinfybound;
};
} // end namespace oofem

#endif // prescribedgradient_h
