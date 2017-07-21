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

#ifndef structuralfe2materialplanestress_h
#define structuralfe2materialplanestress_h

#include "Materials/structuralmaterial.h"
#include "Materials/structuralms.h"

#include <memory>

///@name Input fields for StructuralFE2MaterialPlaneStress
//@{
#define _IFT_StructuralFE2MaterialPlaneStress_Name "structfe2materialplanestress"
#define _IFT_StructuralFE2MaterialPlaneStress_fileName "filename"
#define _IFT_StructuralFE2MaterialPlaneStress_useNumericalTangent "use_num_tangent"
#define _IFT_StructuralFE2MaterialPlaneStress_RegCoeff "regcoeff"
#define _IFT_StructuralFE2MaterialPlaneStress_useExternalStiffness "use_ext_stiffness"
#define _IFT_StructuralFE2MaterialPlaneStress_allGPResults "allgpres"
//@}

namespace oofem {
class EngngModel;
class PrescribedGradientHomogenization;

class StructuralFE2MaterialPlaneStressStatus : public StructuralMaterialStatus
{
protected:
    /// The RVE
    std :: unique_ptr< EngngModel > rve;
    /// Boundary condition in RVE that performs the computational homogenization.
    PrescribedGradientHomogenization *bc;

    FloatMatrix tangent;
    bool oldTangent;


    /// Interface normal direction
    FloatArray mNormalDir;

    std :: string mInputFile;

public:
    StructuralFE2MaterialPlaneStressStatus(int n, int j, Domain * d, GaussPoint * g,  const std :: string & inputfile);
    virtual ~StructuralFE2MaterialPlaneStressStatus() {}

    EngngModel *giveRVE() { return this->rve.get(); }
    PrescribedGradientHomogenization *giveBC();// { return this->bc; }

    void markOldTangent();
    void computeTangent(TimeStep *tStep);

    /// Creates/Initiates the RVE problem.
    bool createRVE(int n, int j, GaussPoint *gp, const std :: string &inputfile);

    /// Copies time step data to RVE.
    void setTimeStep(TimeStep *tStep);

    FloatMatrix &giveTangent() { return tangent; }
    void setTangent(const FloatMatrix &iTangent) {tangent = iTangent; oldTangent = false;}

    virtual const char *giveClassName() const { return "StructuralFE2MaterialPlaneStressStatus"; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    const FloatArray &giveNormal() const { return mNormalDir; }
    void letNormalBe(FloatArray iN) { mNormalDir = std :: move(iN); }

    double giveRveLength();

    /// Functions for MaterialStatusMapperInterface
    virtual void copyStateVariables(const MaterialStatus &iStatus);
    virtual void addStateVariables(const MaterialStatus &iStatus) {OOFEM_ERROR("Not implemented.")};

    // For debugging only
    bool mNewlyInitialized;

};


/**
 * Multiscale constitutive model for subscale structural problems.
 *
 * The material uses the PrescribedGradient boundary condition to perform computational homogenization.
 * The requirement for the supplied subscale problem is:
 * - It must have a PrescribedGradient boundary condition.
 * - It must be the first boundary condition
 *
 * @author Mikael Öhman
 */
class StructuralFE2MaterialPlaneStress : public StructuralMaterial
{
protected:
    std :: string inputfile;
    static int n;
    bool useNumTangent;
    bool useExtStiff;
    bool allGPRes;

    double mRegCoeff;
    FloatMatrix givenTangent; //if the tangent is specified by user

public:
    StructuralFE2MaterialPlaneStress(int n, Domain * d);
    virtual ~StructuralFE2MaterialPlaneStress();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual const char *giveInputRecordName() const { return _IFT_StructuralFE2MaterialPlaneStress_Name; }
    virtual const char *giveClassName() const { return "StructuralFE2MaterialPlaneStress"; }
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return true; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
    // stress computation methods
    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);
    virtual void giveRealStressVector_StressControl(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, const IntArray &strainControl, TimeStep *tStep);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep);
};

} // end namespace oofem
#endif // structuralfe2material_h
