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

#include "prescribedslipgradbcneumann.h"
#include "classfactory.h"
#include "node.h"
#include "masterdof.h"
#include "element.h"
#include "Elements/structuralelement.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "gausspoint.h"
#include "sparsemtrx.h"
#include "xfem/xfemelementinterface.h"
#include "xfem/integrationrules/discsegintegrationrule.h"
#include "timestep.h"
#include "function.h"
#include "sparselinsystemnm.h"
#include "unknownnumberingscheme.h"
#include "engngm.h"
#include "mathfem.h"
#include "crosssection.h"

namespace oofem {
REGISTER_BoundaryCondition(PrescribedSlipGradientBCNeumann);

PrescribedSlipGradientBCNeumann :: PrescribedSlipGradientBCNeumann(int n, Domain *d) :
    ActiveBoundaryCondition(n, d),
    PrescribedFieldsGradientsHomogenization(),
    lmSigmaSHom( new Node(0, d) )
{
}

DofManager* PrescribedSlipGradientBCNeumann :: giveInternalDofManager(int i)
{
    return lmSigmaSHom.get();
}

IRResultType PrescribedSlipGradientBCNeumann :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    int ver;

    IR_GIVE_FIELD(ir, rebarSets, _IFT_PrescribedSlipGradientBCNeumann_RebarSets);
    Domain *d = this->giveDomain();
    int nsd = d->giveNumberOfSpatialDimensions();

    ver = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, ver, _IFT_PrescribedSlipGradientBCNeumann_Version);
    if ( ver == 1 ) {
        this->modelVersion = sg_noShear;
        for ( int i = 0; i < nsd; i++) {
            int dofId = d->giveNextFreeDofID();
            lmSigmaSIds.followedBy(dofId);
            lmSigmaSHom->appendDof( new MasterDof ( lmSigmaSHom.get(), ( DofIDItem ) ( dofId ) ) );
        }
    } else if ( ver == 0 ) {
        this->modelVersion = sg_shear;
        for ( int i = 0; i < nsd*nsd; i++) {
            int dofId = d->giveNextFreeDofID();
            lmSigmaSIds.followedBy(dofId);
            lmSigmaSHom->appendDof( new MasterDof ( lmSigmaSHom.get(), ( DofIDItem ) ( dofId ) ) );
        }
    } else {
        OOFEM_WARNING("unknown version");
        return IRRT_BAD_FORMAT;
    }

    ActiveBoundaryCondition :: initializeFrom(ir);
    return PrescribedFieldsGradientsHomogenization :: initializeFrom(ir);
}

void PrescribedSlipGradientBCNeumann :: scale(double s)
{
    this->mGradient2.times(s);
}

void PrescribedSlipGradientBCNeumann :: giveInputRecord(DynamicInputRecord &input)
{
    ActiveBoundaryCondition :: giveInputRecord(input);
    PrescribedFieldsGradientsHomogenization :: giveInputRecord(input);
}

void PrescribedSlipGradientBCNeumann :: assembleVector(FloatArray &answer, TimeStep *tStep,
                        CharType type, ValueModeType mode,
                        const UnknownNumberingScheme &s, FloatArray *eNorm)
{
    IntArray sigmas_loc;
    lmSigmaSHom->giveLocationArray(lmSigmaSIds, sigmas_loc, s);

    if ( type == ExternalForcesVector ) {
        // The external forces have two contributions. On the additional equations for sigmaS, the load is equal to the prescribed slip gradient value.
        FloatArray stressLoad;
        FloatArray slipGrad;
        this->giveSecondGradient(slipGrad);

        if ( this->modelVersion == sg_noShear ) {
            slipGrad.resizeWithValues(2);
        }

        double loadLevel = this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
        stressLoad.beScaled(loadLevel, slipGrad);

        answer.assemble(stressLoad, sigmas_loc);
    } else if ( type == InternalForcesVector ) {
        FloatMatrix Ksig, Ksig_el;
        FloatArray fe_sig, fe_u;
        FloatArray sigmaS, u_s, u_c;
        IntArray loc, loc_s, loc_c, masterDofIDs;

        // Fetch the current values of the Lagrange multiplier;
        lmSigmaSHom->giveUnknownVector(sigmaS, lmSigmaSIds, mode, tStep);

        // CONTRIBUTION FROM REINFORCEMENT
        FloatMatrix C;
        double gammaBoxInt = computeInterfaceLength(rebarSets);
        this->computeWeightMatrix(C, rebarSets);
        Ksig.clear();
        Ksig_el.clear();
        fe_sig.clear();
        fe_u.clear();

        for ( int i = 0; i < rebarSets.giveSize(); ++i ) {
            Set *steelSet = this->giveDomain()->giveSet(rebarSets.at(i+1));

            for (int pos = 1; pos <= steelSet->giveElementList().giveSize(); ++pos) {
                Element *es = this->giveDomain()->giveElement( steelSet->giveElementList().at(pos) );

                // Fetch the element information;
                es->giveLocationArray(loc_s, s, &masterDofIDs);

                // Fetch the nodal displacements
                // since computeVectorOf gives the unknown vector in local coordinate system, it needs to be rotated to global coordinates
                FloatMatrix G2L;
                es->computeGtoLRotationMatrix(G2L);
                es->computeVectorOf(mode, tStep, u_s);
                u_s.rotatedWith(G2L, 't');

                //Compute the stiffness matrix expansion
                this->integrateTangentOnSteel(Ksig_el, es, rebarSets.at(i+1));
                Ksig.beTProductOf(C, Ksig_el);
                Ksig.times(1/gammaBoxInt);

                //Compute the contribution to internal force vector
                fe_u.beTProductOf(Ksig, sigmaS);
                fe_sig.beProductOf(Ksig, u_s);

                //Assemble
                answer.assemble(fe_u, loc_s);
                answer.assemble(fe_sig, sigmas_loc);
            }
        }

        Ksig.clear();
        fe_sig.clear();
        fe_u.clear();

        // CONTRIBUTION FROM CONCRETE
        double omegaBox = this->domainSize(this->giveDomain(), this->giveSetNumber());
        Set *concreteBoundSet = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = concreteBoundSet->giveBoundaryList();

        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *ec = this->giveDomain()->giveElement( boundaries.at(pos*2-1) );
            int boundary = boundaries.at(pos*2);

            // Fetch the element information;
            ec->giveLocationArray(loc_c, s, &masterDofIDs);

            // Fetch the nodal displacements
            ec->computeVectorOf(mode, tStep, u_c);

            //Compute the stiffness matrix expansion
            this->integrateTangentOnConcrete(Ksig, ec, boundary);
            Ksig.times(1/omegaBox);

            //Compute the contribution to internal force vector
            fe_u.beTProductOf(Ksig, sigmaS);
            fe_sig.beProductOf(Ksig, u_c);

            // Note: The terms appear negative in the equations:
            fe_u.negated();
            fe_sig.negated();

            //Assemble
            answer.assemble(fe_u, loc_c);
            answer.assemble(fe_sig, sigmas_loc);
        }
    }
}

void PrescribedSlipGradientBCNeumann :: assemble(SparseMtrx &answer, TimeStep *tStep,
                  CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale)
{
    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == ElasticStiffnessMatrix ) {
        FloatMatrix Ke, KeT;
        FloatMatrix Ksig;
        IntArray loc_r, loc_c, sigmaS_loc_r, sigmaS_loc_c;

        // Fetch the columns/rows for the stress contributions;
        lmSigmaSHom->giveLocationArray(lmSigmaSIds, sigmaS_loc_r, r_s);
        lmSigmaSHom->giveLocationArray(lmSigmaSIds, sigmaS_loc_c, c_s);

        //Contribution from reinforcement
        FloatMatrix C;
        double gammaBoxInt = computeInterfaceLength(rebarSets);
        this->computeWeightMatrix(C, rebarSets);

        for ( int i = 0; i < rebarSets.giveSize(); ++i ) {
            Set *steelSet = this->giveDomain()->giveSet(rebarSets.at(i+1));

            for ( int pos = 1; pos <= steelSet->giveElementList().giveSize(); ++pos ) {
                Element *es = this->giveDomain()->giveElement( steelSet->giveElementList().at(pos) );

                es->giveLocationArray(loc_r, r_s);
                es->giveLocationArray(loc_c, c_s);

                this->integrateTangentOnSteel(Ksig, es, rebarSets.at(i+1));
                Ke.beTProductOf(C, Ksig);
                Ke.times(1/gammaBoxInt);
                KeT.beTranspositionOf(Ke);

                answer.assemble(sigmaS_loc_r, loc_c, Ke);
                answer.assemble(loc_r, sigmaS_loc_c, KeT);
            }
        }

        Ksig.clear();
        Ke.clear();
        KeT.clear();

        //Contribution from concrete
        Set *concreteBoundSet = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = concreteBoundSet->giveBoundaryList();
        double omegaBox = this->domainSize(this->giveDomain(), this->giveSetNumber());

        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *ec = this->giveDomain()->giveElement( boundaries.at(pos*2-1) );
            int boundary = boundaries.at(pos*2);

            ec->giveLocationArray(loc_r, r_s);
            ec->giveLocationArray(loc_c, c_s);

            this->integrateTangentOnConcrete(Ke, ec, boundary);
            Ke.times(1/omegaBox);
            Ke.negated();
            KeT.beTranspositionOf(Ke);

            answer.assemble(sigmaS_loc_r, loc_c, Ke);
            answer.assemble(loc_r, sigmaS_loc_c, KeT);
        }
    } else   {
        OOFEM_LOG_DEBUG("Skipping assembly in PrescribedSlipGradientBCNeumann::assemble().");
    }
}

void PrescribedSlipGradientBCNeumann :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    IntArray loc_r, loc_c, sigmaS_loc_r, sigmaS_loc_c;

    // Fetch the columns/rows for the reinforcement stress contributions;
    lmSigmaSHom->giveLocationArray(lmSigmaSIds, sigmaS_loc_r, r_s);
    lmSigmaSHom->giveLocationArray(lmSigmaSIds, sigmaS_loc_c, c_s);

    Set *concreteSet = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = concreteSet->giveBoundaryList();
    IntArray steelEl;

    for ( int i = 1; i <= this->rebarSets.giveSize(); ++i ) {
        steelEl.followedBy(this->giveDomain()->giveSet(rebarSets.at(i))->giveElementList());
    }

    rows.resize( 2*steelEl.giveSize() + boundaries.giveSize() );
    cols.resize( 2*steelEl.giveSize() + boundaries.giveSize() );
    int i = 0;
    //Contribution from steel
    for ( int pos = 1; pos <= steelEl.giveSize(); ++pos ) {
        Element *e = this->giveDomain()->giveElement(steelEl.at(pos));

        e->giveLocationArray(loc_r, r_s);
        e->giveLocationArray(loc_c, c_s);

        // For most uses, loc_r == loc_c, and sigmaS_loc_r == sigmaS_loc_c.
        rows [ i ] = loc_r;
        cols [ i ] = sigmaS_loc_c;
        i++;
        // and the symmetric part (usually the transpose of above)
        rows [ i ] = sigmaS_loc_r;
        cols [ i ] = loc_c;
        i++;
    }

    //Contribution from concrete
    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos*2-1) );

        e->giveLocationArray(loc_r, r_s);
        e->giveLocationArray(loc_c, c_s);

        // For most uses, loc_r == loc_c, and sigmaS_loc_r == sigmaS_loc_c.
        rows [ i ] = loc_r;
        cols [ i ] = sigmaS_loc_c;
        i++;
        // and the symmetric part (usually the transpose of above)
        rows [ i ] = sigmaS_loc_r;
        cols [ i ] = loc_c;
        i++;
    }
}

void PrescribedSlipGradientBCNeumann :: giveReinfStressLocationArray(IntArray &oCols, const UnknownNumberingScheme &r_s)
{
    lmSigmaSHom->giveLocationArray(lmSigmaSIds, oCols, r_s);
}

void PrescribedSlipGradientBCNeumann :: computeReinfStress(FloatArray &rStress, TimeStep *tStep)
{
    double volRVE = 1.;

    if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 2 ) {
        //assuming that the RVE thickness is constant (2D)
        Element *e = this->giveDomain()->giveElement(this->giveDomain()->giveSet(this->set)->giveBoundaryList().at(1));
        std :: unique_ptr< IntegrationRule > ir = e->giveInterpolation()->giveIntegrationRule(e->giveInterpolation()->giveInterpolationOrder());
        CrossSection *cs = e->giveCrossSection();
        GaussPoint *gp = ir->getIntegrationPoint(1);
        double thick = cs->give(CS_Thickness, gp);
        double omegaBox = this->domainSize(this->giveDomain(), this->giveSetNumber());
        volRVE = omegaBox * thick;
    } else if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 3 ) {
        volRVE = this->domainSize(this->giveDomain(), this->giveSetNumber());
    }

    lmSigmaSHom->giveUnknownVector(rStress, lmSigmaSIds, VM_Total, tStep);
    rStress.times(1 / volRVE );
    rStress.negated();

    if ( modelVersion == sg_noShear ) {
        rStress.resizeWithValues(4);
    }
}

void PrescribedSlipGradientBCNeumann :: computeTangent(FloatMatrix &tangent, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented");
}

void PrescribedSlipGradientBCNeumann :: computeRebarDyad(FloatMatrix &dyad, const int reinfSet)
{
    ///Assuming that the reinforcement bar remains straight throughout the RVE
    ///it is then enough to compute the tangent vector for the first element

    FloatArray e_l;
    Set *set = this->giveDomain()->giveSet(reinfSet);
    IntArray ellist = set->giveElementList();
    Element *element = this->giveDomain()->giveElement(ellist.at(1));
    FEIElementGeometryWrapper cellgeo(element);

    e_l.resize(2);
    e_l.at(1) = cellgeo.giveVertexCoordinates(2).at(1) - cellgeo.giveVertexCoordinates(1).at(1);
    e_l.at(2) = cellgeo.giveVertexCoordinates(2).at(2) - cellgeo.giveVertexCoordinates(1).at(2);
    e_l.normalize();

    dyad.beDyadicProductOf(e_l,e_l);
}

void PrescribedSlipGradientBCNeumann :: computeWeightMatrix(FloatMatrix& C, const IntArray reinfSets)
{
    ///Assumption: the rebars cross the entire RVE, i.e. they do not end within the RVE
    /// - the rebars have circular cross section
    double gammaBoxInt = computeInterfaceLength(reinfSets);
    FloatMatrix dyadMatrix;
    for ( int i = 0; i < reinfSets.giveSize(); ++i ) {
        FloatMatrix dya;
        double perimeterCoeff = 0;

        Set *set = this->giveDomain()->giveSet(reinfSets.at(i+1));
        for ( int j = 0; j < set->giveElementList().giveSize(); ++j ) {
            Element *e = this->giveDomain()->giveElement(set->giveElementList().at(j+1));
            FEInterpolation *interp = e->giveInterpolation();
            int order = interp->giveInterpolationOrder();
            std :: unique_ptr< IntegrationRule > ir;
            ir = interp->giveIntegrationRule(order);

            for ( auto &gp : *ir ) {
                const FloatArray &lcoords = gp->giveNaturalCoordinates();
                FEIElementGeometryWrapper cellgeo(e);
                double detJ = interp->giveTransformationJacobian(lcoords, cellgeo);
                CrossSection *cs = e->giveCrossSection();
                double perim = sqrt(4 * cs->give(CS_Area, gp) * M_PI);
                perimeterCoeff = perimeterCoeff + perim * detJ * gp->giveWeight();
            }
        }

        this->computeRebarDyad(dya, reinfSets.at(i+1));
        dyadMatrix.add(perimeterCoeff, dya);
    }

    C.beInverseOf(dyadMatrix);

    if ( this->modelVersion == sg_shear ) {
        C.resizeWithData(4,4);
    }

    C.times(gammaBoxInt);
}

double PrescribedSlipGradientBCNeumann :: computeInterfaceLength(const IntArray reinfSets)
{
    double gamma = 1.;
    for ( int i=0; i < reinfSets.giveSize(); ++i ) {
        Set *set = this->giveDomain()->giveSet(reinfSets.at(i+1));
        for ( int j=0; j < set->giveElementList().giveSize(); ++j ) {
            Element *e = this->giveDomain()->giveElement(set->giveElementList().at(j+1));
            FEInterpolation *interp = e->giveInterpolation();
            int order = interp->giveInterpolationOrder();
            std :: unique_ptr< IntegrationRule > ir;
            ir = interp->giveIntegrationRule(order);

            for ( auto &gp : *ir ) {
                const FloatArray &lcoords = gp->giveNaturalCoordinates();
                FEIElementGeometryWrapper cellgeo(e);
                double detJ = interp->giveTransformationJacobian(lcoords, cellgeo);
                gamma = gamma + detJ * gp->giveWeight();
            }
        }
    }

    return gamma;
}

void PrescribedSlipGradientBCNeumann :: integrateTangentOnSteel(FloatMatrix &oTangent, Element *e, const int rebarSet)
{
    //reminder: Ksig must be 4x6
    FloatMatrix Nssig, Bs, contrib;
    double Ss;
    IntArray DofMask;

    //Get number of DOFs per DofManager
    e->giveDofManDofIDMask(1, DofMask);

    FEInterpolation *interp = e->giveInterpolation();
    int order = interp->giveInterpolationOrder();
    std :: unique_ptr< IntegrationRule > ir;
    ir = interp->giveIntegrationRule(order);
    //Cast into StructuralElement
    StructuralElement *se = dynamic_cast<StructuralElement*>(e);

    oTangent.clear();

    for ( auto &gp : *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        //Compute shape function derivatives
        FloatMatrix b;
        se->computeBmatrixAt(gp, b);
        //TODO: make it element independent; now hardcoded for libeam2d elements
        Bs.resize(2,6);
        Bs.at(1,1) = b.at(1,1);
        Bs.at(2,2) = b.at(1,1);
        Bs.at(1,4) = b.at(1,4);
        Bs.at(2,5) = b.at(1,4);

        //Constant traction interpolation
        this->computeRebarDyad(Nssig, rebarSet);
        if ( this->modelVersion == sg_shear ) {
            Nssig.resizeWithData(2,4);
        }

        //Compute rebar perimeter
        CrossSection *cs = e->giveCrossSection();
        Ss = sqrt(4 * cs->give(CS_Area, gp) * M_PI);

        //Compute the integral
        contrib.beTProductOf(Nssig, Bs);
        contrib.times(Ss);
        double detJ = interp->giveTransformationJacobian(lcoords, cellgeo);
        oTangent.add(detJ * gp->giveWeight(), contrib);
    }
}


void PrescribedSlipGradientBCNeumann :: integrateTangentOnConcrete(FloatMatrix &oTangent, Element *e, int iBndIndex)
{
    //reminder: Ksig must be 4x8
    FloatArray normal, n;
    FloatMatrix nMatrix, E_n;
    FloatMatrix contrib;

    FEInterpolation *interp = e->giveInterpolation();

    int nsd = e->giveDomain()->giveNumberOfSpatialDimensions();

    //Interpolation order
    int order = interp->giveInterpolationOrder();
    std :: unique_ptr< IntegrationRule > ir = interp->giveBoundaryEdgeIntegrationRule(order, iBndIndex);

    oTangent.clear();

    for ( auto &gp : *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        //Evaluate the normal
        double detJ = interp->boundaryEvalNormal(normal, iBndIndex, lcoords, cellgeo);

        //Compute global coordinates of Gauss point
        FloatArray globalCoord;
        interp->boundaryLocal2Global(globalCoord, iBndIndex, lcoords, cellgeo);

        //Compute local coordinates on the element
        FloatArray bulkElLocCoords;
        e->computeLocalCoordinates(bulkElLocCoords, globalCoord);

        //Evaluate the shape functions
        interp->evalN(n, bulkElLocCoords, cellgeo);
        nMatrix.beNMatrixOf(n, nsd);

        if ( this->modelVersion == sg_shear ) {
            E_n.resize(nsd*nsd, nsd);
            E_n.zero();
            if ( nsd == 3 ) {
                E_n.at(1, 1) = normal.at(1);
                E_n.at(2, 2) = normal.at(2);
                E_n.at(3, 3) = normal.at(3);
                E_n.at(4, 1) = normal.at(2);
                E_n.at(5, 1) = normal.at(3);
                E_n.at(6, 2) = normal.at(3);
                E_n.at(7, 2) = normal.at(1);
                E_n.at(8, 3) = normal.at(1);
                E_n.at(9, 3) = normal.at(2);
            } else if ( nsd == 2 ) {
                E_n.at(1, 1) = normal.at(1);
                E_n.at(2, 2) = normal.at(2);
                E_n.at(3, 1) = normal.at(2);
                E_n.at(4, 2) = normal.at(1);
            } else {
                E_n.at(1, 1) = normal.at(1);
            }
        } else {  //sg_noShear
            E_n.resize(nsd, nsd);
            E_n.zero();
            if ( nsd == 3 ) {
                E_n.at(1, 1) = normal.at(1);
                E_n.at(2, 2) = normal.at(2);
                E_n.at(3, 3) = normal.at(3);
            } else if ( nsd == 2 ) {
                E_n.at(1, 1) = normal.at(1);
                E_n.at(2,2) = normal.at(2);
            } else {
                E_n.at(1, 1) = normal.at(1);
            }
        }
        //AS: According to OOFEM convention, the normal is computed as positive inwards, it needs to be negated here to make it consistent with the RVE boundary
        E_n.negated();
        contrib.beProductOf(E_n, nMatrix);

        oTangent.add(detJ * gp->giveWeight(), contrib);
    }
}

} /* namespace oofem */
