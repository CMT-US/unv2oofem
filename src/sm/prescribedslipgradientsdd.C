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

#include "prescribedslipgradientsdd.h"
#include "dofiditem.h"
#include "dofmanager.h"
#include "dof.h"
#include "valuemodetype.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "function.h"
#include "engngm.h"
#include "set.h"
#include "node.h"
#include "element.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "feinterpol.h"
#include "unknownnumberingscheme.h"
#include "sparsemtrx.h"
#include "sparselinsystemnm.h"
#include "assemblercallback.h"
#include "mathfem.h"

namespace oofem {
REGISTER_BoundaryCondition(PrescribedSlipGradientsDD);

double PrescribedSlipGradientsDD :: give(Dof *dof, ValueModeType mode, double time)
{
    DofIDItem id = dof->giveDofID();
    FloatArray *coords = dof->giveDofManager()->giveCoordinates();

    if ( coords->giveSize() != this->mCenterCoord.giveSize() ) {
        mCenterCoord.resizeWithValues(coords->giveSize());
//        OOFEM_ERROR("Size of coordinate system different from center coordinate in b.c.");
//        printf("Warning: Size of coordinate system different from center coordinate in b.c.\n");
    }

    double factor = 0;
    if ( mode == VM_Total ) {
        factor = this->giveTimeFunction()->evaluateAtTime(time);
    } else if ( mode == VM_Velocity ) {
        factor = this->giveTimeFunction()->evaluateVelocityAtTime(time);
    } else if ( mode == VM_Acceleration ) {
        factor = this->giveTimeFunction()->evaluateAccelerationAtTime(time);
    } else {
        OOFEM_ERROR("Should not be called for value mode type then total, velocity, or acceleration.");
    }
    // Reminder: u_i = d_ij . (x_j - xb_j) = d_ij . dx_j
    FloatArray dx;
    dx.beDifferenceOf(* coords, this->mCenterCoord);

    mGradient1.resizeWithData(coords->giveSize(), coords->giveSize());
    mField1.resizeWithValues(coords->giveSize());

    FloatArray u;
    u.beProductOf(mGradient1, dx);
    u.at(1)+=mField1.at(1);
    u.at(2)+=mField1.at(2);
    u.times( factor );

    ///@todo Use the user-specified dofs here instead:
    int pos = this->dofs.findFirstIndexOf(id);

    //Prescribing rotations at beams' ends. Defining sets and checking whether the node belongs to one of them
    //Works only in 2d, i.e. for a reinforcement with 3 degrees of freedom
    Set* setx = domain->giveSet(reinfxbound);
    Set* sety= domain->giveSet(reinfybound);
    int isXReinf=0;
    int isYReinf=0;
    FloatArray eL, ePerp, us;
    double sField=0.0;
    eL.resize(coords->giveSize());
    ePerp.resize(coords->giveSize());
    eL.zero();
    ePerp.zero();
    mGradient2.resizeWithData(coords->giveSize(), coords->giveSize());
    mField2.resizeWithValues(coords->giveSize());

    if (setx->giveNodeList().contains(dof->giveDofManager()->giveGlobalNumber()) ) {
        isXReinf=1;
        eL.at(1)=1;
        eL.at(2)=0;
        ePerp.at(1)=0;
        ePerp.at(2)=1;
        us.beProductOf(mGradient2, dx);
        sField = eL.dotProduct(mField2);
    } else if (sety->giveNodeList().contains(dof->giveDofManager()->giveGlobalNumber()) ) {
        isYReinf=1;
        eL.at(1)=0;
        eL.at(2)=1;
        ePerp.at(1)=-1;
        ePerp.at(2)=0;
        us.beProductOf(mGradient2, dx);
        sField= eL.dotProduct(mField2);
    } else {
        isXReinf=0;
        isYReinf=0;
    }

    if(pos > 0 && pos <= u.giveSize()) {
        if (isXReinf) {
            if (pos == 1) {return u.at(pos) += sField + us.at(pos);}
            else if (pos == 2 ) {return u.at(pos) += us.at(pos);}
            else if ( pos == 3) {return u.at(pos) = -mGradient1.at(2,1);} //@todo: contribution from slip gradient?
            else { return 0.0;}
        } else if (isYReinf) {
            if (pos == 1) {return u.at(pos) += us.at(pos);}
            else if (pos == 2) {return u.at(pos) += sField + us.at(pos);} //@todo: contribution from slip gradient?
            else if (pos == 3) {return u.at(pos) = mGradient1.at(2,1);}
            else { return 0.0;}
        } else {
        return u.at(pos);
        }
    } else {
    	// XFEM dofs
    	return 0.0;
    }
}

void PrescribedSlipGradientsDD :: updateCoefficientMatrix(FloatMatrix &C)
// This is written in a very general way, supporting both fm and sm problems.
// v_prescribed = C.d = (x-xbar).d;
// Modified by ES.
// C = [x 0 0 y]
//     [0 y x 0]
//     [ ... ] in 2D, voigt form [d_11, d_22, d_12 d_21]
// C = [x 0 0 0 z y 0 0 0]
//     [0 y 0 z 0 0 0 0 x]
//     [0 0 z 0 0 0 y x 0]
//     [ ............... ] in 3D, voigt form [d_11, d_22, d_33, d_23, d_13, d_12, d_32, d_31, d_21]
//
//Modified by AS:
//Include end moments from the reinforcement. Formula according to Sciegaj et al. Two-scale finite element modelling
//of reinforced concrete structures: effective response and sub-scale fracture.
//\sum (R_L e_l + R_perp e_{\perp} ) \outerp (x-\bar{x}) already included in C^T.R_c (in computeField)
//Added term \sum R_M e_{\perp} \outerp e_l
// C = [x                0              0               y  ]
//     [0                y              x               0  ]
//     [ePerp1*eL1   ePerp2*eL2    ePerp1*eL2    ePerp2*eL1]
//  for DofManagers with rotational degrees of freedom. For DofManagers with only translational dofs,
//  C matrix is as before (2 first rows)
{
    Domain *domain = this->giveDomain();

    int nsd = domain->giveNumberOfSpatialDimensions();
    int npeq = domain->giveEngngModel()->giveNumberOfDomainEquations( domain->giveNumber(), EModelDefaultPrescribedEquationNumbering() );
    C.resize(npeq, nsd * nsd);
    C.zero();

    FloatArray &cCoords = this->giveCenterCoordinate();
    double xbar = cCoords.at(1), ybar = cCoords.at(2), zbar = 0.0;
    if ( nsd == 3 ) {
        zbar = cCoords.at(3);
    }

    for ( auto &n : domain->giveDofManagers() ) {
        FloatArray *coords = n->giveCoordinates();
        Dof *d1 = n->giveDofWithID( this->dofs(0) );
        Dof *d2 = n->giveDofWithID( this->dofs(1) );
        //check if the DofManager has a rotational degree of freedom
        int krot=0;
        FloatArray eL, ePerp;
        eL.resize(nsd);
        ePerp.resize(nsd);
        eL.zero();
        ePerp.zero();
        if (n->giveNumberOfDofs() == 3) {
            Dof *drot = n->giveDofWithID( this->dofs(2) );
            krot = drot->__givePrescribedEquationNumber();

            //Compute tangential and normal basis vectors - supported only for orthogonal reinforcement in x and y direction
            Set* setx = domain->giveSet(reinfxbound);
            Set* sety = domain->giveSet(reinfybound);

            if ( setx->giveNodeList().contains(n->giveGlobalNumber()) ) {
                eL.at(1)=1;
                eL.at(2)=0;
                ePerp.at(1)=0;
                ePerp.at(2)=1;
            }

            if ( sety->giveNodeList().contains(n->giveGlobalNumber()) ) {
                eL.at(1)=0;
                eL.at(2)=1;
                ePerp.at(1)=-1;
                ePerp.at(2)=0;
            }

        }

        int k1 = d1->__givePrescribedEquationNumber();
        int k2 = d2->__givePrescribedEquationNumber();
        if ( nsd == 2 ) {
            if ( k1 ) {
                C.at(k1, 1) = coords->at(1) - xbar;
                C.at(k1, 4) = coords->at(2) - ybar;
            }

            if ( k2 ) {
                C.at(k2, 2) = coords->at(2) - ybar;
                C.at(k2, 3) = coords->at(1) - xbar;
            }

            if ( krot ) {
                C.at(krot,1) = ePerp.at(1)*eL.at(1);
                C.at(krot,2) = ePerp.at(2)*eL.at(2);
                C.at(krot,3) = ePerp.at(1)*eL.at(2);
                C.at(krot,4) = ePerp.at(2)*eL.at(1);
            }
        } else { // nsd == 3
            Dof *d3 = n->giveDofWithID( this->dofs(2) );
            int k3 = d3->__givePrescribedEquationNumber();

            if ( k1 ) {
                C.at(k1, 1) = coords->at(1) - xbar;
                C.at(k1, 6) = coords->at(2) - ybar;
                C.at(k1, 5) = coords->at(3) - zbar;
            }
            if ( k2 ) {
                C.at(k2, 2) = coords->at(2) - ybar;
                C.at(k2, 9) = coords->at(1) - xbar;
                C.at(k2, 4) = coords->at(3) - zbar;
            }
            if ( k3 ) {
                C.at(k3, 3) = coords->at(3) - zbar;
                C.at(k3, 8) = coords->at(1) - xbar;
                C.at(k3, 7) = coords->at(2) - ybar;
            }
        }
    }
}


void PrescribedSlipGradientsDD :: computeStress(FloatArray &sigma, TimeStep *tStep)
{
    EngngModel *emodel = this->domain->giveEngngModel();
    int npeq = emodel->giveNumberOfDomainEquations( this->giveDomain()->giveNumber(), EModelDefaultPrescribedEquationNumbering() );
    FloatArray R_c(npeq), R_ext(npeq);

    R_c.zero();
    R_ext.zero();
    emodel->assembleVector( R_c, tStep, InternalForceAssembler(), VM_Total,
                            EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );
    //R_c contains reactions at the boundary - coming from concrete tractions, beam end forces (normal and shear) and end moments

    emodel->assembleVector( R_ext, tStep, ExternalForceAssembler(), VM_Total,
                            EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );

    R_c.subtract(R_ext);

    // Condense it;
    FloatMatrix C;
    this->updateCoefficientMatrix(C);
    sigma.beTProductOf(C, R_c);
    sigma.times( 1. / (this->domainSize(this->giveDomain(), conboundset) * thick ) );
}

void PrescribedSlipGradientsDD :: computeBondStress(FloatArray &bStress, TimeStep *tStep)
{
     //According to 1/(\Omega_\Box) * \sum R_L \be{l}

    EngngModel *emodel = this->domain->giveEngngModel();
    Domain *domain = this->giveDomain();

    int nsd = domain->giveNumberOfSpatialDimensions();
    int npeq = domain->giveEngngModel()->giveNumberOfDomainEquations( domain->giveNumber(), EModelDefaultPrescribedEquationNumbering() );

    FloatArray R_c(npeq), R_ext(npeq);

    R_c.zero();
    R_ext.zero();
    emodel->assembleVector( R_c, tStep, InternalForceAssembler(), VM_Total,
                            EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );
    //R_c contains reactions at the boundary - coming from concrete tractions, beam end forces (normal and shear) and end moments

    emodel->assembleVector( R_ext, tStep, ExternalForceAssembler(), VM_Total,
                            EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );

    R_c.subtract(R_ext);

    //taking contribution only from end nodes of the reinforcement
    IntArray nodesX = domain->giveSet(reinfxbound)->giveNodeList();
    IntArray nodesY = domain->giveSet(reinfybound)->giveNodeList();
    FloatArray eL;
    double R_L;
    eL.resize(nsd); eL.zero();
    bStress.resize(nsd); bStress.zero();

    for (int i=1; i<=nodesX.giveSize(); i++) {
        DofManager *d1 = domain->giveDofManager( nodesX.at(i) );
        eL.at(1) = 1;
        eL.at(2) = 0;
        R_L = R_c.at( d1->giveDofWithID(1)->__givePrescribedEquationNumber() );
        bStress.at(1) += R_L*eL.at(1);
        bStress.at(2) += R_L*eL.at(2);
    }

    for (int i=1; i<=nodesY.giveSize(); i++) {
        DofManager *d2 = domain->giveDofManager( nodesY.at(i) );
        eL.at(1) = 0;
        eL.at(2) = 1;
        R_L = R_c.at( d2->giveDofWithID(2)->__givePrescribedEquationNumber() );
        bStress.at(1) += R_L*eL.at(1);
        bStress.at(2) += R_L*eL.at(2);
    }

    bStress.times( 1. / (this->domainSize(this->giveDomain(), conboundset) * thick ) );
}

void PrescribedSlipGradientsDD :: computeReinfStress(FloatArray &rStress, TimeStep *tStep)
{
    //According to 1/(\Omega_\Box) * \sum R_L \be{l} \outerp [x - \bar{x} ]

    EngngModel *emodel = this->domain->giveEngngModel();
    Domain *domain = this->giveDomain();

    int nsd = domain->giveNumberOfSpatialDimensions();
    int npeq = domain->giveEngngModel()->giveNumberOfDomainEquations( domain->giveNumber(), EModelDefaultPrescribedEquationNumbering() );

    FloatArray R_c(npeq), R_ext(npeq);

    R_c.zero();
    R_ext.zero();
    emodel->assembleVector( R_c, tStep, InternalForceAssembler(), VM_Total,
                            EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );
    //R_c contains reactions at the boundary - coming from concrete tractions, beam end forces (normal and shear) and end moments

    emodel->assembleVector( R_ext, tStep, ExternalForceAssembler(), VM_Total,
                            EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );

    R_c.subtract(R_ext);

    //taking contribution only from end nodes of the reinforcement
    IntArray nodesX = domain->giveSet(reinfxbound)->giveNodeList();
    IntArray nodesY = domain->giveSet(reinfybound)->giveNodeList();
    FloatArray eL, cCoords = this->giveCenterCoordinate();
    double xbar=cCoords.at(1), ybar=cCoords.at(2), R_L;
    eL.resize(nsd); eL.zero();
    rStress.resize(nsd*2); rStress.zero();

    for (int i=1; i<=nodesX.giveSize(); i++) {
        DofManager *d1 = domain->giveDofManager( nodesX.at(i) );
        FloatArray *coords = d1->giveCoordinates();
        eL.at(1) = 1;
        eL.at(2) = 0;
        R_L = R_c.at( d1->giveDofWithID(1)->__givePrescribedEquationNumber() );
        rStress.at(1) += R_L*eL.at(1) * (coords->at(1) - xbar);
        rStress.at(2) += R_L*eL.at(2) * (coords->at(2) - ybar);
        rStress.at(3) += R_L*eL.at(1) * (coords->at(2) - ybar);
        rStress.at(4) += R_L*eL.at(2) * (coords->at(1) - xbar);
    }

    for (int i=1; i<=nodesY.giveSize(); i++) {
        DofManager *d2 = domain->giveDofManager( nodesY.at(i) );
        FloatArray *coords = d2->giveCoordinates();
        eL.at(1) = 0;
        eL.at(2) = 1;
        R_L = R_c.at( d2->giveDofWithID(2)->__givePrescribedEquationNumber() );
        rStress.at(1) += R_L*eL.at(1) * (coords->at(1) - xbar);
        rStress.at(2) += R_L*eL.at(2) * (coords->at(2) - ybar);
        rStress.at(3) += R_L*eL.at(1) * (coords->at(2) - ybar);
        rStress.at(4) += R_L*eL.at(2) * (coords->at(1) - xbar);
    }

    rStress.times( 1. / (this->domainSize(this->giveDomain(), conboundset) * thick ) );
}

void PrescribedSlipGradientsDD :: computeTangent(FloatMatrix &tangent, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented.");
}


IRResultType PrescribedSlipGradientsDD :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    GeneralBoundaryCondition :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, conboundset , _IFT_PrescribedSlipGradientsDD_ConcreteBoundary);

    IR_GIVE_FIELD(ir, reinfxbound , _IFT_PrescribedSlipGradientsDD_ReinfXBound);

    IR_GIVE_FIELD(ir, reinfybound , _IFT_PrescribedSlipGradientsDD_ReinfYBound);

    IR_GIVE_FIELD(ir, thick , _IFT_PrescribedSlipGradientsDD_Thickness);

    return PrescribedFieldsGradientsHomogenization :: initializeFrom(ir);
}


void PrescribedSlipGradientsDD :: giveInputRecord(DynamicInputRecord &input)
{
    GeneralBoundaryCondition :: giveInputRecord(input);

    input.setField(conboundset, _IFT_PrescribedSlipGradientsDD_ConcreteBoundary);
    input.setField(thick, _IFT_PrescribedSlipGradientsDD_Thickness);
    input.setField(reinfxbound, _IFT_PrescribedSlipGradientsDD_ReinfXBound);
    input.setField(reinfybound, _IFT_PrescribedSlipGradientsDD_ReinfYBound);
    return PrescribedFieldsGradientsHomogenization :: giveInputRecord(input);
}
} // end namespace oofem
