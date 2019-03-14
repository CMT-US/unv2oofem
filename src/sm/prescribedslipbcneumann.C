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

#include "prescribedslipbcneumann.h"
#include "classfactory.h"
#include "node.h"
#include "masterdof.h"
#include "element.h"
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
REGISTER_BoundaryCondition(PrescribedSlipBCNeumann);

PrescribedSlipBCNeumann :: PrescribedSlipBCNeumann(int n, Domain *d) :
    ActiveBoundaryCondition(n, d),
    PrescribedFieldsGradientsHomogenization(),
    lmTauHom( new Node(0, d) )
{
    int nsd = d->giveNumberOfSpatialDimensions();

    for ( int i = 0; i < nsd; i++) {
        int dofId = d->giveNextFreeDofID();
        lmTauIds.followedBy(dofId);
        lmTauHom->appendDof( new MasterDof ( lmTauHom.get(), ( DofIDItem ) ( dofId ) ) );
    }
}

DofManager* PrescribedSlipBCNeumann :: giveInternalDofManager(int i)
{
    return lmTauHom.get();
}

IRResultType PrescribedSlipBCNeumann :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, concreteVolElSet, _IFT_PrescribedSlipBCNeumann_ConcreteVolElSet);
    IR_GIVE_FIELD(ir, rebarSets, _IFT_PrescribedSlipBCNeumann_RebarSets);

    ActiveBoundaryCondition :: initializeFrom(ir);
    return PrescribedFieldsGradientsHomogenization :: initializeFrom(ir);
}

void PrescribedSlipBCNeumann :: scale(double s)
{
    this->mField2.times(s);
}

void PrescribedSlipBCNeumann :: giveInputRecord(DynamicInputRecord &input)
{
    ActiveBoundaryCondition :: giveInputRecord(input);
    PrescribedFieldsGradientsHomogenization :: giveInputRecord(input);
}

void PrescribedSlipBCNeumann :: assembleVector(FloatArray &answer, TimeStep *tStep,
                        CharType type, ValueModeType mode,
                        const UnknownNumberingScheme &s, FloatArray *eNorm)
{
    IntArray tau_loc;
    lmTauHom->giveLocationArray(lmTauIds, tau_loc, s);

    if ( type == ExternalForcesVector ) {
        // The external forces have two contributions. On the additional equations for tau, the load is equal to the prescribed slip value.
        FloatArray stressLoad;
        FloatArray slipVec;
        this->giveSecondField(slipVec);

        double loadLevel = this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
        stressLoad.beScaled(loadLevel, slipVec);

        answer.assemble(stressLoad, tau_loc);
    } else if ( type == InternalForcesVector ) {
        FloatMatrix Ktau, Ktau_el;
        FloatArray fe_tau, fe_u;
        FloatArray tau, u_s, u_c;
        IntArray loc, loc_s, loc_c, masterDofIDs;

        // Fetch the current values of the Lagrange multiplier;
        lmTauHom->giveUnknownVector(tau, lmTauIds, mode, tStep);

        // CONTRIBUTION FROM REINFORCEMENT
        FloatMatrix C;
        double gammaBoxInt = computeInterfaceLength(rebarSets);
        this->computeWeightMatrix(C, rebarSets);
        Ktau.clear();
        Ktau_el.clear();
        fe_tau.clear();
        fe_u.clear();

        for (int i = 0; i < rebarSets.giveSize(); ++i) {
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
                this->integrateTangentOnSteel(Ktau_el, es, rebarSets.at(i+1));
                Ktau.beTProductOf(C, Ktau_el);
                Ktau.times(1/gammaBoxInt);

                //Compute the contribution to internal force vector
                fe_u.beTProductOf(Ktau, tau);
                fe_tau.beProductOf(Ktau, u_s);

                //Assemble
                answer.assemble(fe_u, loc_s);
                answer.assemble(fe_tau, tau_loc);
            }
        }
        Ktau.clear();
        fe_tau.clear();
        fe_u.clear();

        // CONTRIBUTION FROM CONCRETE
        Set *concreteSet = this->giveDomain()->giveSet(concreteVolElSet);
        double omegaBox = this->domainSize(this->giveDomain(), this->giveSetNumber());

        for ( int pos = 1; pos <= concreteSet->giveElementList().giveSize(); ++pos) {
            Element* ec = this->giveDomain()->giveElement( concreteSet->giveElementList().at(pos) );

            // Fetch the element information;
            ec->giveLocationArray(loc_c, s, &masterDofIDs);

            // Fetch the nodal displacements
            ec->computeVectorOf(mode, tStep, u_c);

            //Compute the stiffness matrix expansion
            this->integrateTangentOnConcrete(Ktau, ec);
            Ktau.times(1/omegaBox);

            //Compute the contribution to internal force vector
            fe_u.beTProductOf(Ktau, tau);
            fe_tau.beProductOf(Ktau, u_c);

            // Note: The terms appear negative in the equations:
            fe_u.negated();
            fe_tau.negated();

            //Assemble
            answer.assemble(fe_u, loc_c);
            answer.assemble(fe_tau, tau_loc);
        }
    }
}

void PrescribedSlipBCNeumann :: assemble(SparseMtrx &answer, TimeStep *tStep,
                  CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale)
{
    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == ElasticStiffnessMatrix ) {
        FloatMatrix Ke, KeT;
        FloatMatrix Ktau;
        IntArray loc_r, loc_c, tau_loc_r, tau_loc_c;

        // Fetch the columns/rows for the stress contributions;
        lmTauHom->giveLocationArray(lmTauIds, tau_loc_r, r_s);
        lmTauHom->giveLocationArray(lmTauIds, tau_loc_c, c_s);

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

                this->integrateTangentOnSteel(Ktau, es, rebarSets.at(i+1));
                Ke.beTProductOf(C, Ktau);
                Ke.times(1/gammaBoxInt);
                KeT.beTranspositionOf(Ke);

                answer.assemble(tau_loc_r, loc_c, Ke);
                answer.assemble(loc_r, tau_loc_c, KeT);
            }
        }

        Ktau.clear();
        Ke.clear();
        KeT.clear();

        //Contribution from concrete
        Set *concreteSet = this->giveDomain()->giveSet(concreteVolElSet);
        double omegaBox = this->domainSize(this->giveDomain(), this->giveSetNumber());

        for ( int pos = 1; pos <= concreteSet->giveElementList().giveSize(); ++pos ) {
            Element *ec = this->giveDomain()->giveElement( concreteSet->giveElementList().at(pos) );

            ec->giveLocationArray(loc_r, r_s);
            ec->giveLocationArray(loc_c, c_s);

            this->integrateTangentOnConcrete(Ke, ec);
            Ke.times(1/omegaBox);
            Ke.negated();
            KeT.beTranspositionOf(Ke);

            answer.assemble(tau_loc_r, loc_c, Ke);
            answer.assemble(loc_r, tau_loc_c, KeT);
        }
    } else   {
        OOFEM_LOG_DEBUG("Skipping assembly in PrescribedSlipBCNeumann::assemble().");
    }
}

void PrescribedSlipBCNeumann :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    IntArray loc_r, loc_c, tau_loc_r, tau_loc_c;

    // Fetch the columns/rows for the transfer stress contributions;
    lmTauHom->giveLocationArray(lmTauIds, tau_loc_r, r_s);
    lmTauHom->giveLocationArray(lmTauIds, tau_loc_c, c_s);

    Set *concreteSet = this->giveDomain()->giveSet(this->concreteVolElSet);
    const IntArray &conEl = concreteSet->giveElementList();
    IntArray steelEl;

    for ( int i = 1; i <= this->rebarSets.giveSize(); ++i ) {
        steelEl.followedBy(this->giveDomain()->giveSet(rebarSets.at(i))->giveElementList());
    }

    rows.resize( 2*steelEl.giveSize() + 2*conEl.giveSize() );
    cols.resize( 2*steelEl.giveSize() + 2*conEl.giveSize() );
    int i = 0;
    //Contribution from steel
    for ( int pos = 1; pos <= steelEl.giveSize(); ++pos ) {
        Element *e = this->giveDomain()->giveElement(steelEl.at(pos));

        e->giveLocationArray(loc_r, r_s);
        e->giveLocationArray(loc_c, c_s);

        // For most uses, loc_r == loc_c, and tau_loc_r == tau_loc_c.
        rows [ i ] = loc_r;
        cols [ i ] = tau_loc_c;
        i++;
        // and the symmetric part (usually the transpose of above)
        rows [ i ] = tau_loc_r;
        cols [ i ] = loc_c;
        i++;
    }

    //Contribution from concrete
    for ( int pos = 1; pos <= conEl.giveSize(); ++pos ) {
        Element *e = this->giveDomain()->giveElement(conEl.at(pos));

        e->giveLocationArray(loc_r, r_s);
        e->giveLocationArray(loc_c, c_s);

        // For most uses, loc_r == loc_c, and tau_loc_r == tau_loc_c.
        rows [ i ] = loc_r;
        cols [ i ] = tau_loc_c;
        i++;
        // and the symmetric part (usually the transpose of above)
        rows [ i ] = tau_loc_r;
        cols [ i ] = loc_c;
        i++;
    }
}

void PrescribedSlipBCNeumann :: giveTransferStressLocationArray(IntArray &oCols, const UnknownNumberingScheme &r_s)
{
    lmTauHom->giveLocationArray(lmTauIds, oCols, r_s);
}

void PrescribedSlipBCNeumann :: computeTransferStress(FloatArray &bStress, TimeStep *tStep)
{
    //scaling needed - divide by the volume of the RVE
    double volRVE = 1.;

    if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 2 ) {
        //assuming that the RVE thickness is constant (2D)
        Element *e = this->giveDomain()->giveElement(this->giveDomain()->giveSet(concreteVolElSet)->giveElementList().at(1));
        std :: unique_ptr< IntegrationRule > ir = e->giveInterpolation()->giveIntegrationRule(e->giveInterpolation()->giveInterpolationOrder());
        CrossSection *cs = e->giveCrossSection();
        GaussPoint *gp = ir->getIntegrationPoint(1);
        double thick = cs->give(CS_Thickness, gp);
        double omegaBox = this->domainSize(this->giveDomain(), this->giveSetNumber());
        volRVE = omegaBox * thick;
    } else if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 3 ) {
        volRVE = this->domainSize(this->giveDomain(), this->giveSetNumber());
    }

    // Transfer stress is computed from the tractions IN THE INTERFACE
    // Lagrange multiplier here corresponds to distributed load ALONG THE REBARS / IN THE CONCRETE
    // therefore need to change sign
    lmTauHom->giveUnknownVector(bStress, lmTauIds, VM_Total, tStep);
    bStress.times( 1/ volRVE );
    bStress.negated();
}

void PrescribedSlipBCNeumann :: computeTangent(FloatMatrix &tangent, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented");
}

void PrescribedSlipBCNeumann :: computeRebarDyad(FloatMatrix &dyad, const int reinfSet)
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

void PrescribedSlipBCNeumann :: computeWeightMatrix(FloatMatrix& C, const IntArray reinfSets)
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
    C.times(gammaBoxInt);
}

double PrescribedSlipBCNeumann :: computeInterfaceLength(const IntArray reinfSets)
{
    double gamma = 0.;
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

void PrescribedSlipBCNeumann :: integrateTangentOnSteel(FloatMatrix &oTangent, Element *e, const int rebarSet)
{
    //reminder: Ktau must be 2x6
    FloatMatrix Nstau, Ns, contrib;
    double Ss;
    IntArray DofMask;

    //Get number of DOFs per DofManager
    e->giveDofManDofIDMask(1, DofMask);
    int ndof = DofMask.giveSize();

    FEInterpolation *interp = e->giveInterpolation();
    int order = interp->giveInterpolationOrder();
    std :: unique_ptr< IntegrationRule > ir;
    ir = interp->giveIntegrationRule(order);

    oTangent.clear();

    for ( auto &gp : *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        //Compute shape functions
        FloatArray nu;
        interp->evalN(nu, lcoords, cellgeo);
        Ns.beNMatrixOf(nu, ndof);

        //Constant traction interpolation
        this->computeRebarDyad(Nstau, rebarSet);
        Nstau.resizeWithData(3,2);

        //Compute rebar perimeter
        CrossSection *cs = e->giveCrossSection();
        Ss = sqrt(4 * cs->give(CS_Area, gp) * M_PI);

        //Compute the integral
        contrib.beTProductOf(Nstau, Ns);
        contrib.times(Ss);
        double detJ = interp->giveTransformationJacobian(lcoords, cellgeo);
        oTangent.add(detJ * gp->giveWeight(), contrib);
    }
}

void PrescribedSlipBCNeumann :: integrateTangentOnConcrete(FloatMatrix &oTangent, Element *e)
{
    //reminder: Ktau must be 2x8
    FloatMatrix Nctau, Nc, contrib;
    IntArray DofMask;

    //Get number of DOFs per DofManager
    e->giveDofManDofIDMask(1, DofMask);
    int ndof = DofMask.giveSize();

    FEInterpolation *interp = e->giveInterpolation();
    int order = interp->giveInterpolationOrder();
    std :: unique_ptr< IntegrationRule > ir;
    ir = interp->giveIntegrationRule(order);

    oTangent.clear();
    Nctau.resize(ndof,ndof);

    for ( auto &gp : *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        //Compute shape functions
        FloatArray nu;
        interp->evalN(nu, lcoords, cellgeo);
        Nc.beNMatrixOf(nu, ndof);

        //Constant traction interpolation
        Nctau.beUnitMatrix();

        //Compute the integral
        contrib.beTProductOf(Nctau, Nc);
        double detJ = interp->giveTransformationJacobian(lcoords, cellgeo);
        oTangent.add(detJ * gp->giveWeight(), contrib);
    }
}

} /* namespace oofem */
