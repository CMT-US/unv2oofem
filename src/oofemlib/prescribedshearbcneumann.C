/*
 * prescribedshearbcneumann.C
 *
 *  Created on: Mar 7, 2017
 *      Author: svennine
 */

#include "prescribedshearbcneumann.h"

#include "node.h"
#include "classfactory.h"
#include "domain.h"
#include "masterdof.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "function.h"
#include "timestep.h"
#include "element.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "xfem/xfemelementinterface.h"
#include "xfem/integrationrules/discsegintegrationrule.h"
#include "gausspoint.h"
#include "sparsemtrx.h"

namespace oofem {
REGISTER_BoundaryCondition(PrescribedShearBCNeumann);

PrescribedShearBCNeumann::PrescribedShearBCNeumann(int n, Domain *d) :
ActiveBoundaryCondition(n, d),
mShearStrain(0.),
mpShearStress( new Node(0, d) )
{

	printf("Entering PrescribedShearBCNeumann::PrescribedShearBCNeumann(int n, Domain *d)\n");

	// Add new dof for the Lagrange multiplier
    int dofId = d->giveNextFreeDofID();
    mShearStressDofIds.followedBy(dofId);
    mpShearStress->appendDof( new MasterDof( mpShearStress.get(), ( DofIDItem ) ( dofId ) ) );

    if(d) {
    	d->computeDomainBoundingBox(mLC, mUC);
    }
}

PrescribedShearBCNeumann::~PrescribedShearBCNeumann() {

}

IRResultType PrescribedShearBCNeumann::initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, mShearStrain, _IFT_PrescribedShearBCNeumann_Shear);
    printf("mShearStrain: %e\n", mShearStrain);

	return ActiveBoundaryCondition :: initializeFrom(ir);
}

void PrescribedShearBCNeumann::giveInputRecord(DynamicInputRecord &input)
{

    ActiveBoundaryCondition :: giveInputRecord(input);
}

DofManager *PrescribedShearBCNeumann :: giveInternalDofManager(int i)
{
    return mpShearStress.get();
}

void PrescribedShearBCNeumann :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                                   CharType type, ValueModeType mode,
                                                   const UnknownNumberingScheme &s, FloatArray *eNorm)
{
    IntArray stress_loc;  // For the displacements and stress respectively
    mpShearStress->giveLocationArray(mShearStressDofIds, stress_loc, s);

	double L = giveL();
	double H = giveH();

    if ( type == ExternalForcesVector ) {
        // External force contribution for the shear stress LMP.
        FloatArray stressLoad;

        double loadLevel = this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
        stressLoad = {-H*L*loadLevel*mShearStrain};

        answer.assemble(stressLoad, stress_loc);
    } else if ( type == InternalForcesVector ) {

    	FloatMatrix Ke_us;
        FloatArray f_int_u, f_int_s;
        FloatArray stress_y, e_u;
        IntArray loc, masterDofIDs, sigmaMasterDofIDs;

        // Fetch the current value of the shear stress;
        mpShearStress->giveUnknownVector(stress_y, mShearStressDofIds, mode, tStep);

        // Assemble by looping over all boundary segments
        Set *setPointer = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = setPointer->giveBoundaryList();
        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            // Fetch the element information;
            e->giveLocationArray(loc, s, & masterDofIDs);

            // Here, we could use only the nodes actually located on the boundary, but we don't.
            // Instead, we use all nodes belonging to the element, which is allowed because the
            // basis functions related to the interior nodes will be zero on the boundary.
            e->computeVectorOf(mode, tStep, e_u);

            this->integrateTangent(Ke_us, e, boundary);

            // We just use the tangent, less duplicated code (the LMP contribution is linear).
            f_int_s.beTProductOf(Ke_us, e_u);
            f_int_u.beProductOf(Ke_us, stress_y);

            // Note: The terms appear negative in the equations:
            f_int_u.negated();
            f_int_s.negated();

            answer.assemble(f_int_u, loc); // Contributions to delta_v equations
            answer.assemble(f_int_s, stress_loc); // Contribution to delta_s_i equations
        }
    }

}

void PrescribedShearBCNeumann :: assemble(SparseMtrx &answer, TimeStep *tStep,
                                             CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == ElasticStiffnessMatrix ) {
        FloatMatrix Ke, KeT;
        IntArray disp_loc_r, disp_loc_c, stress_loc_r, stress_loc_c;
        Set *set = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = set->giveBoundaryList();

        // Fetch the columns/rows for the stress contributions;
        mpShearStress->giveLocationArray(mShearStressDofIds, stress_loc_r, r_s);
        mpShearStress->giveLocationArray(mShearStressDofIds, stress_loc_c, c_s);

        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            e->giveLocationArray(disp_loc_r, r_s);
            e->giveLocationArray(disp_loc_c, c_s);

            this->integrateTangent(Ke, e, boundary);
            Ke.negated();
            KeT.beTranspositionOf(Ke);

            answer.assemble(disp_loc_r, stress_loc_c, Ke); // Contributions to displacement equations
            answer.assemble(stress_loc_r, disp_loc_c, KeT); // Contribution to stress equations
        }


        // Make sure that all diagonal entries are present to make the linear solver happy.
        FloatMatrix Ke_diag(1,1);
        Ke_diag(0,0) = 0.;
        answer.assemble(stress_loc_r, stress_loc_c, Ke_diag);

    } else   {
        OOFEM_LOG_DEBUG("Skipping assembly in PrescribedShearBCNeumann::assemble().");
    }
}

void PrescribedShearBCNeumann :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                                       const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    IntArray loc_r, loc_c, sigma_loc_r, sigma_loc_c;

    // Fetch the columns/rows for the stress contributions;
    mpShearStress->giveLocationArray(mShearStressDofIds, sigma_loc_r, r_s);
    mpShearStress->giveLocationArray(mShearStressDofIds, sigma_loc_c, c_s);

    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    rows.resize( boundaries.giveSize() );
    cols.resize( boundaries.giveSize() );
    int i = 0;
    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );

        e->giveLocationArray(loc_r, r_s);
        e->giveLocationArray(loc_c, c_s);

        // For most uses, loc_r == loc_c, and sigma_loc_r == sigma_loc_c.
        rows [ i ] = loc_r;
        cols [ i ] = sigma_loc_c;
        i++;
        // and the symmetric part (usually the transpose of above)
        rows [ i ] = sigma_loc_r;
        cols [ i ] = loc_c;
        i++;
    }
}



void PrescribedShearBCNeumann :: integrateTangent(FloatMatrix &oTangent, Element *e, int iBndIndex)
{
    FloatArray normal, nu;
    FloatMatrix Nu, E_n;
    FloatMatrix contrib;

    // Constant traction interpolation
    FloatMatrix Nt(2,1);
    Nt(0,0) = 0.;
    Nt(1,0) = 1.;

    Domain *domain = e->giveDomain();

    FEInterpolation *interp = e->giveInterpolation(); // Geometry interpolation

    int nsd = 2;

    // Interpolation order
    int order = interp->giveInterpolationOrder();
    std :: unique_ptr< IntegrationRule > ir;

    // Fetch integration rule
    XfemElementInterface *xfemElInt = dynamic_cast< XfemElementInterface * >( e );
    if ( xfemElInt != NULL && domain->hasXfemManager() ) {
        IntArray edgeNodes;
        FEInterpolation2d *interp2d = dynamic_cast< FEInterpolation2d * >( interp );
        if ( interp2d == NULL ) {
            OOFEM_ERROR("failed to cast to FEInterpolation2d.")
        }
        interp2d->computeLocalEdgeMapping(edgeNodes, iBndIndex);

        std :: vector< Line >segments;
        std :: vector< FloatArray >intersecPoints;
        xfemElInt->partitionEdgeSegment(iBndIndex, segments, intersecPoints);
        MaterialMode matMode = e->giveMaterialMode();
        ir.reset( new DiscontinuousSegmentIntegrationRule(1, e, segments) );
        int numPointsPerSeg = 1;
        ir->SetUpPointsOnLine(numPointsPerSeg, matMode);
    } else {
        ir = interp->giveBoundaryIntegrationRule(order, iBndIndex);
    }

    oTangent.clear();

    // Loop over Gauss points
    for ( auto &gp: *ir ) {

    	// Check if the current GP is on Gamma^+, Gamma^-, top, or bottom.
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        FloatArray gcoords;
        interp->boundaryLocal2Global(gcoords, iBndIndex, lcoords, cellgeo);

        FloatArray el_lcoords;
        interp->global2local(el_lcoords, gcoords, cellgeo);


    	double c = 0.0;
    	if( pointIsOnRightEdge(gcoords) ) {
    		c = 1.0;
    	}

    	if( pointIsOnLeftEdge(gcoords) ) {
    		c = -1.0;
    	}



        // Evaluate the normal;
        double detJ = interp->boundaryEvalNormal(normal, iBndIndex, lcoords, cellgeo);

        interp->evalN(nu, el_lcoords, cellgeo);
        // If cracks cross the edge, special treatment is necessary.
        // Exploit the XfemElementInterface to minimize duplication of code.
        if ( xfemElInt != NULL && domain->hasXfemManager() ) {
            // Compute global coordinates of Gauss point
            FloatArray globalCoord;

            interp->boundaryLocal2Global(globalCoord, iBndIndex, lcoords, cellgeo);

            // Compute local coordinates on the element
            FloatArray locCoord;
            e->computeLocalCoordinates(locCoord, globalCoord);

            xfemElInt->XfemElementInterface_createEnrNmatrixAt(Nu, locCoord, * e, false);
        } else {
            // Evaluate the velocity/displacement coefficients
            Nu.beNMatrixOf(nu, nsd);
        }

        // Clear first row in Nu, because only the y-component should be included.
        for( int i = 0; i < Nu.giveNumberOfColumns(); i++ ) {
        	Nu(0,i) = 0.;
        }

        contrib.beTProductOf(Nu, Nt);
        oTangent.add(c*detJ * gp->giveWeight(), contrib);
    }

}

} /* namespace oofem */
