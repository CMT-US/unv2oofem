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

/*
 * trustregionsolver.h
 *
 *  Created on: Feb 24, 2017
 *      Author: svennine
 */

#ifndef LINESEARCHSOLVER_H_
#define LINESEARCHSOLVER_H_


#define _IFT_TrustRegionSolver_Name "trustregionsolver"


#include "nrsolver.h"

namespace oofem {

/**
 *
 * Trust-region algorithm.
 * Implementation based on the "Basic trust-region algorithm"
 * on p.116 in Conn et al. (2000).
 *
 */
class OOFEM_EXPORT TrustRegionSolver : public NRSolver {
public:
	TrustRegionSolver(Domain * d, EngngModel * m);
	virtual ~TrustRegionSolver();

    virtual NM_Status solve(SparseMtrx &k, FloatArray &R, FloatArray *R0,
                            FloatArray &X, FloatArray &dX, FloatArray &F,
                            const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                            int &nite, TimeStep *);

    bool checkConvergence(FloatArray &RT, FloatArray &F, FloatArray &rhs, FloatArray &ddX, FloatArray &X,
                          double RRT, const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange, bool printToScreen = true);

protected:

    /// Trust region parameters
    double mEta1, mEta2, mGamma1, mGamma2;

    /// Trust region size
    double mTrustRegionSize;

};

} /* namespace oofem */

#endif /* LINESEARCHSOLVER_H_ */
