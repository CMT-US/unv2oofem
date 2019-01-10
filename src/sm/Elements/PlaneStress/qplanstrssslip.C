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

#include "sm/Elements/PlaneStress/qplanstrssslip.h"
#include "fei2dquadquad.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"
#include "Materials/structuralmaterial.h"
#include "Materials/structuralms.h"
#include "Materials/structuralfe2materialplanestress.h"
#include "CrossSections/structuralcrosssection.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(QPlaneStress2dSlip);

FEI2dQuadQuad QPlaneStress2dSlip :: interpolation(1, 2);

QPlaneStress2dSlip :: QPlaneStress2dSlip(int n, Domain *aDomain) :
    PlaneStressElement(n, aDomain), ZZNodalRecoveryModelInterface(this), NodalAveragingRecoveryModelInterface()
    // Constructor.
{
    numberOfDofMans  = 8;
    numberOfGaussPoints = 4;
}

void QPlaneStress2dSlip :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, S_u, S_v
    };
}

Interface *
QPlaneStress2dSlip :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    }

    return NULL;
}

FEInterpolation *QPlaneStress2dSlip :: giveInterpolation() const { return & interpolation; }


void QPlaneStress2dSlip :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    int ndof=this->computeNumberOfDofs();

    FloatMatrix Kaa, Kab, Kba, Kbb, daa, dab1, dab2, dba1, dba2, dbb1, dbb2, dbb3, dbb4;

    for ( auto &gp : *this->giveDefaultIntegrationRulePtr() ) {
        FloatMatrix Ba, Bb, Nb;
        double dV;
        this->computeBmatrixAt(gp, Ba);
        this->computeBHmatrixAt(gp, Bb);
        this->computeNmatrixAt(gp->giveNaturalCoordinates(), Nb);
        dV = this->computeVolumeAround(gp);

        //Compute the sensitivities
        FloatMatrix dStressdEps, dStressdS, dStressdG, dBStressdEps, dBStressdS, dBStressdG, dRStressdEps, dRStressdS, dRStressdG;
        StructuralFE2MaterialPlaneStress *mat = dynamic_cast<StructuralFE2MaterialPlaneStress *>( this->giveStructuralCrossSection()->giveMaterial(gp) );
        if ( mat ) {
            mat->computeSensitivities(dStressdEps, dStressdS, dStressdG, dBStressdEps, dBStressdS, dBStressdG, dRStressdEps, dRStressdS, dRStressdG, rMode, gp, tStep);
//             printf("dSdE \n");
//             dStressdEps.printYourself();
//             printf("dBSdE \n");
//             dBStressdEps.printYourself();
//             printf("dRSdE \n");
//             dRStressdEps.printYourself();
//             printf("dSdS \n");
//             dStressdS.printYourself();
//             printf("dBSdS \n");
//             dBStressdS.printYourself();
//             printf("dRSdS \n");
//             dRStressdS.printYourself();
//             printf("dSdG \n");
//             dStressdG.printYourself();
//             printf("dBSdG \n");
//             dBStressdG.printYourself();
//             printf("dRSdG \n");
//             dRStressdG.printYourself();
        } else {
            OOFEM_ERROR("Wrong material, this element works only with StructuralFE2MaterialPlaneStress");
        }

        daa.beProductOf(dStressdEps, Ba);
        dab1.beProductOf(dStressdS,Nb);
        dab2.beProductOf(dStressdG,Bb);
        dba1.beProductOf(dBStressdEps,Ba);
        dba2.beProductOf(dRStressdEps,Ba);
        dbb1.beProductOf(dBStressdS,Nb);
        dbb2.beProductOf(dBStressdG,Bb);
        dbb3.beProductOf(dRStressdS,Nb);
        dbb4.beProductOf(dRStressdG,Bb);

        Kaa.plusProductUnsym(Ba, daa, dV);
        Kab.plusProductUnsym(Ba, dab1, dV);
        Kab.plusProductUnsym(Ba, dab2, dV);
        Kba.plusProductUnsym(Nb,dba1,dV);
        Kba.plusProductUnsym(Bb,dba2,dV);
        Kbb.plusProductUnsym(Nb,dbb1,dV);
        Kbb.plusProductUnsym(Nb,dbb2,dV);
        Kbb.plusProductUnsym(Bb,dbb3,dV);
        Kbb.plusProductUnsym(Bb,dbb4,dV);
    }

    //there's probably a cleaner way to split these fields
    IntArray aMask = {1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30};
    IntArray bMask = {3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32};

    answer.resize(ndof,ndof);
    answer.assemble(Kaa,aMask,aMask);
    answer.assemble(Kbb,bMask,bMask);
    answer.assemble(Kab,aMask,bMask);
    answer.assemble(Kba,bMask,aMask);
}

void QPlaneStress2dSlip :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    int ndof=this->computeNumberOfDofs();

    //Compute local fields at nodes
    FloatArray u, a, b;
    this->computeVectorOf(VM_Total, tStep, u);
    //Split into displacement and slip fields;
    IntArray aMask = {1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30};
    IntArray bMask = {3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32};
    a.beSubArrayOf(u,aMask);
    b.beSubArrayOf(u,bMask);

    FloatArray vStrain, vSlip, vSlipGradient;
    FloatArray Stress, bStress, rStress;
    FloatArray finta, fintb;

    for ( auto &gp : *this->giveDefaultIntegrationRulePtr() ) {
        FloatMatrix Ba, Bb, Nb;
        FloatArray f;
        double dV;
        this->computeBmatrixAt(gp, Ba);
        this->computeBHmatrixAt(gp, Bb);
        this->computeNmatrixAt(gp->giveNaturalCoordinates(), Nb);
        dV = this->computeVolumeAround(gp);

        vStrain.beProductOf(Ba, a);
        vSlip.beProductOf(Nb, b);
        vSlipGradient.beProductOf(Bb, b);

        this->computeHomogenizedFields(Stress, bStress, rStress, vStrain, vSlip, vSlipGradient, gp , tStep);

        finta.plusProduct(Ba,Stress,dV);
        fintb.plusProduct(Nb,bStress,dV);
        fintb.plusProduct(Bb,rStress,dV);
    }

    answer.resize(ndof);
    for (int i=1; i <= ndof/2 ; i++) {
        if (i % 2 == 1)  {
            answer.at(2*i-1) = finta.at(i);
            answer.at(2*i+1) = fintb.at(i);
        } else if (i % 2 == 0) {
            answer.at(2*i-2) = finta.at(i);
            answer.at(2*i) = fintb.at(i);
        }
    }
}

void QPlaneStress2dSlip :: computeHomogenizedFields(FloatArray &Stress, FloatArray &bStress, FloatArray &rStress, const FloatArray &strain, const FloatArray &slip, const FloatArray &slipGradient, GaussPoint *gp, TimeStep *tStep)
{
    //Homogenize stress, transfer stress and reinforcement stress from the RVE
    StructuralFE2MaterialPlaneStress *mat = dynamic_cast<StructuralFE2MaterialPlaneStress *>( this->giveStructuralCrossSection()->giveMaterial(gp) );
    if ( mat ) {
        mat->giveHomogenizedFields(Stress, bStress, rStress, strain, slip, slipGradient, gp, tStep);
    } else {
        OOFEM_ERROR("Check material, this element works only with StructuralFE2MaterialPlaneStress");
    }
}

void
QPlaneStress2dSlip :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                             InternalStateType type, TimeStep *tStep)
{
    if ( numberOfGaussPoints == 4 ) {

        GaussPoint *gp;

        if ( node < 5 ) {
            int i = 0;
            switch ( node ) {
            case 1: i = 4;
                break;
            case 2: i = 2;
                break;
            case 3: i = 1;
                break;
            case 4: i = 3;
                break;
            }

            gp = integrationRulesArray [ 0 ]->getIntegrationPoint(i - 1);
            this->giveIPValue(answer, gp, type, tStep);
        } else {
            int i1 = 0, i2 = 0;
            switch ( node ) {
            case 5: i1 = 4;
                i2 = 2;
                break;
            case 6: i1 = 2;
                i2 = 1;
                break;
            case 7: i1 = 1;
                i2 = 3;
                break;
            case 8: i1 = 3;
                i2 = 4;
                break;
            }

            FloatArray contrib;
            gp = integrationRulesArray [ 0 ]->getIntegrationPoint(i1 - 1);
            this->giveIPValue(contrib, gp, type, tStep);
            gp = integrationRulesArray [ 0 ]->getIntegrationPoint(i2 - 1);
            this->giveIPValue(answer, gp, type, tStep);
            answer.add(contrib);
            answer.times(0.5);
        }
    } else if ( numberOfGaussPoints == 9 ) {

        GaussPoint *gp;

        int i=0;
        switch ( node ) {
        case 1: i = 9;
            break;
        case 2: i = 3;
            break;
        case 3: i = 1;
            break;
        case 4: i = 7;
            break;
        case 5: i = 6;
            break;
        case 6: i = 2;
            break;
        case 7: i = 4;
            break;
        case 8: i = 8;
            break;
        }

        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(i - 1);
        this->giveIPValue(answer, gp, type, tStep);
    } else {
        return;
    }
}

int QPlaneStress2dSlip :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_DisplacementVector ) {
        FloatArray u, a;
        FloatMatrix N;
        IntArray aMask = {1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30};
        this->computeVectorOf(VM_Total, tStep, u);
        a.beSubArrayOf(u, aMask);
        this->computeNmatrixAt(gp->giveSubPatchCoordinates(), N);
        answer.beProductOf(N, a);
        return 1;
    } else if ( type == IST_ShearSlip ) {
        FloatArray u, b;
        FloatMatrix N;
        IntArray bMask = {3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32};
        this->computeVectorOf(VM_Total, tStep, u);
        b.beSubArrayOf(u,bMask);
        this->computeNmatrixAt(gp->giveSubPatchCoordinates(), N);
        answer.beProductOf(N, b);
        return 1;
    } else if ( type == IST_TransferStress ) {
        StructuralFE2MaterialPlaneStressStatus* status=static_cast< StructuralFE2MaterialPlaneStressStatus * >( gp->giveMaterialStatus() );
        answer = status->giveTransferStressVector();
        return 1;
    } else if ( type == IST_ShearSlipGradient ) {
        FloatArray u, b;
        FloatMatrix B;
        IntArray bMask = {3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32};
        this->computeVectorOf(VM_Total, tStep, u);
        b.beSubArrayOf(u,bMask);
        this->computeBHmatrixAt(gp, B);
        answer.beProductOf(B, b);
        return 1;
    } else if (type == IST_ReinforcementMembraneStress ) {
        StructuralFE2MaterialPlaneStressStatus* status=static_cast< StructuralFE2MaterialPlaneStressStatus *>( gp->giveMaterialStatus() );
        answer = status->giveReinfStressVector();
        return 1;
    } else {
        return Element :: giveIPValue(answer, gp, type, tStep);
    }
}


} // end namespace oofem
