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
#define _IFT_StructuralFE2MaterialPlaneStress_dStressdEps "dsde"
#define _IFT_StructuralFE2MaterialPlaneStress_dBStressdEps "dbsde"
#define _IFT_StructuralFE2MaterialPlaneStress_dRStressdEps "drsde"
#define _IFT_StructuralFE2MaterialPlaneStress_dStressdS "dsds"
#define _IFT_StructuralFE2MaterialPlaneStress_dBStressdS "dbsds"
#define _IFT_StructuralFE2MaterialPlaneStress_dRStressdS "drsds"
#define _IFT_StructuralFE2MaterialPlaneStress_dStressdG "dsdg"
#define _IFT_StructuralFE2MaterialPlaneStress_dBStressdG "dbsdg"
#define _IFT_StructuralFE2MaterialPlaneStress_dRStressdG "drsdg"
//@}

namespace oofem {
class EngngModel;
class PrescribedGradientHomogenization;
class PrescribedFieldsGradientsHomogenization;

class StructuralFE2MaterialPlaneStressStatus : public StructuralMaterialStatus
{
protected:
    /// The RVE
    std :: shared_ptr< EngngModel > rve;
    /// Boundary condition in RVE that performs the computational homogenization.
    PrescribedGradientHomogenization *bc;
    /// Boundary condition with two independent fields
    PrescribedFieldsGradientsHomogenization *bc2f;

    FloatMatrix tangent;
    bool oldTangent;

    //Tangents wrt to strain, slip and slip gradient
    FloatMatrix dStressdEpsTangent, dBStressdEpsTangent, dRStressdEpsTangent;
    bool olddSdETangent, olddBSdETangent, olddRSdETangent;

    FloatMatrix dStressdSTangent, dBStressdSTangent, dRStressdSTangent;
    bool olddSdSTangent, olddBSdSTangent, olddRSdSTangent;

    FloatMatrix dStressdGTangent, dBStressdGTangent, dRStressdGTangent;
    bool olddSdGTangent, olddBSdGTangent, olddRSdGTangent;

    /// Interface normal direction
    FloatArray mNormalDir;

    std :: string mInputFile;

    /// Equilibrated slip vector in reduced form
    FloatArray slipVector;
    /// Equilibrated transfer stress vector in reduced form
    FloatArray bStressVector;
    /// Temporary bond stress vector in reduced form (increments are used mainly in nonlinear analysis)
    FloatArray tempbStressVector;
    /// Temporary slip vector in reduced form (to find balanced state)
    FloatArray tempSlipVector;

    /// Equilibrated slip gradient vector in reduced form
    FloatArray slipGradVector;
    /// Equilibrated reinforcement stress vector in reduced form
    FloatArray rStressVector;
    /// Temporary reinforcement stress vector in reduced form (increments are used mainly in nonlinear analysis)
    FloatArray temprStressVector;
    /// Temporary slip gradient vector in reduced form (to find balanced state)
    FloatArray tempSlipGradVector;

public:
    StructuralFE2MaterialPlaneStressStatus(int n, int j, Domain * d, GaussPoint * g,  const std :: string & inputfile);
    virtual ~StructuralFE2MaterialPlaneStressStatus() {};

    EngngModel *giveRVE() { return this->rve.get(); }
    PrescribedGradientHomogenization *giveBC();// { return this->bc; }
    PrescribedFieldsGradientsHomogenization *giveBC2F(); // { return this->bc2f; }

    void markOldTangent();
    void computeTangent(TimeStep *tStep);

    /// Creates/Initiates the RVE problem.
    bool createRVE(int n, int j, GaussPoint *gp, const std :: string &inputfile);

    /// Copies time step data to RVE.
    void setTimeStep(TimeStep *tStep);

    FloatMatrix &giveTangent() { return tangent; }
    void setTangent(const FloatMatrix &iTangent) {tangent = iTangent; oldTangent = false;}

    FloatMatrix &givedStressdEpsTangent() { return dStressdEpsTangent; }
    void setdStressdEpsTangent(const FloatMatrix &iTangent) {dStressdEpsTangent = iTangent; olddSdETangent = false;}

    FloatMatrix &givedBStressdEpsTangent() { return dBStressdEpsTangent; }
    void setdBStressdEpsTangent(const FloatMatrix &iTangent) {dBStressdEpsTangent = iTangent; olddBSdETangent = false;}

    FloatMatrix &givedRStressdEpsTangent() { return dRStressdEpsTangent; }
    void setdRStressdEpsTangent(const FloatMatrix &iTangent) {dRStressdEpsTangent = iTangent; olddRSdETangent = false;}

    FloatMatrix &givedStressdSTangent() { return dStressdSTangent; }
    void setdStressdSTangent(const FloatMatrix &iTangent) {dStressdSTangent = iTangent; olddSdSTangent = false;}

    FloatMatrix &givedBStressdSTangent() { return dBStressdSTangent; }
    void setdBStressdSTangent(const FloatMatrix &iTangent) {dBStressdSTangent = iTangent; olddBSdSTangent = false;}

    FloatMatrix &givedRStressdSTangent() { return dRStressdSTangent; }
    void setdRStressdSTangent(const FloatMatrix &iTangent) {dRStressdSTangent = iTangent; olddRSdSTangent = false;}

    FloatMatrix &givedStressdGTangent() { return dStressdGTangent; }
    void setdStressdGTangent(const FloatMatrix &iTangent) {dStressdGTangent = iTangent; olddSdGTangent = false;}

    FloatMatrix &givedBStressdGTangent() { return dBStressdGTangent; }
    void setdBStressdGTangent(const FloatMatrix &iTangent) {dBStressdGTangent = iTangent; olddBSdGTangent = false;}

    FloatMatrix &givedRStressdGTangent() { return dRStressdGTangent; }
    void setdRStressdGTangent(const FloatMatrix &iTangent) {dRStressdGTangent = iTangent; olddRSdGTangent = false;}

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

    /// Returns the const pointer to receiver's slip vector.
    const FloatArray &giveSlipVector() const { return slipVector; }
    /// Returns the const pointer to receiver's bond stress vector.
    const FloatArray &giveTransferStressVector() const { return bStressVector; }
    /// Returns the const pointer to receiver's temporary slip vector.
    const FloatArray &giveTempSlipVector() const { return tempSlipVector; }
    /// Returns the const pointer to receiver's temporary bond stress vector.
    const FloatArray &giveTempTransferStressVector() const { return tempbStressVector; }
    /// Assigns slip vector to given vector v.
    void letSlipVectorBe(const FloatArray &v) { slipVector = v; }
    /// Assigns bond stress vector to given vector v.
    void letTransferStressVectorBe(const FloatArray &v) { bStressVector = v; }
    /// Assigns tempbStressVector to given vector v.
    void letTempTransferStressVectorBe(const FloatArray &v) { tempbStressVector = v; }
    /// Assigns tempSlipVector to given vector v
    void letTempSlipVectorBe(const FloatArray &v) { tempSlipVector = v; }

    /// Returns the const pointer to receiver's slip vector.
    const FloatArray &giveSlipGradVector() const { return slipGradVector; }
    /// Returns the const pointer to receiver's bond stress vector.
    const FloatArray &giveReinfStressVector() const { return rStressVector; }
    /// Returns the const pointer to receiver's temporary slip vector.
    const FloatArray &giveTempSlipGradVector() const { return tempSlipGradVector; }
    /// Returns the const pointer to receiver's temporary bond stress vector.
    const FloatArray &giveTempReinfStressVector() const { return temprStressVector; }
    /// Assigns slip vector to given vector v.
    void letSlipGradVectorBe(const FloatArray &v) { slipGradVector = v; }
    /// Assigns bond stress vector to given vector v.
    void letReinfStressVectorBe(const FloatArray &v) { rStressVector = v; }
    /// Assigns tempbStressVector to given vector v.
    void letTempReinfStressVectorBe(const FloatArray &v) { temprStressVector = v; }
    /// Assigns tempSlipVector to given vector v
    void letTempSlipGradVectorBe(const FloatArray &v) { tempSlipGradVector = v; }

};


/**
 * Multiscale constitutive model for subscale structural problems.
 *
 * The material uses the PrescribedGradient (or PrescribedSlipGradients) condition to perform computational homogenization.
 * The requirement for the supplied subscale problem is:
 * - It must have a PrescribedGradient boundary condition.
 * - It must be the first boundary condition
 *
 * @author Mikael Öhman
 * @author Adam Sciegaj
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
    FloatMatrix givendStressdEpsTangent, givendBStressdEpsTangent, givendRStressdEpsTangent;
    FloatMatrix givendStressdSTangent, givendBStressdSTangent, givendRStressdSTangent;
    FloatMatrix givendStressdGTangent, givendBStressdGTangent, givendRStressdGTangent;

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

    //for fe2 analyses with independent slip field
    /**
     * Computes the homogenized stress, homogenized bond stress and reinforcement stress for the createRVE
     */
    virtual void giveHomogenizedFields(FloatArray &Stress, FloatArray &bStress, FloatArray &rStress, const FloatArray &strain, const FloatArray &slip,
                                       const FloatArray &slipGradient, GaussPoint *gp, TimeStep *tStep);
    /**
     * Computes the sensitivity matrices for the RVE. Necessary for building the stiffness matrix.
     */
    virtual void computeSensitivities(FloatMatrix &dStressdEps, FloatMatrix &dStressdS, FloatMatrix &dStressdG, FloatMatrix &dBStressdEps,
                                      FloatMatrix &dBStressdS, FloatMatrix &dBStressdG, FloatMatrix &dRStressdEps, FloatMatrix &dRStressdS,
                                      FloatMatrix &dRStressdG, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
};

} // end namespace oofem
#endif // structuralfe2material_h
