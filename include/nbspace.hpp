// Name:     nbspace.hpp
// Author:   A. M. Kazachkov
// Date:     2019-Feb-20
//-----------------------------------------------------------------------------

#pragma once

#include <vector>

// COIN-OR
#include <OsiSolverInterface.hpp>

// VPC
#include "VPCParameters.hpp"

void setCompNBCoorPoint(CoinPackedVector& vec, double& objViolation,
    const VPCParameters& params, const OsiSolverInterface* const tmpSolver,
    const OsiSolverInterface* const origSolver,
    const std::vector<int>& nonBasicVarIndex,
    const std::vector<double>& nonBasicReducedCost, const int deletedVar = -1);
void setCompNBCoorRay(CoinPackedVector& vec, const double* ray, double& objViolation, double& scale,
    const VPCParameters& params, const OsiSolverInterface* const tmpSolver,
    const OsiSolverInterface* const origSolver,
    const std::vector<int>& rowOfOrigNBVar,
    const std::vector<int>& nonBasicVarIndex,
    const std::vector<double>& nonBasicReducedCost,
    const int tmpNBVar, const int deletedVar = -1,
    const bool rayNeedsCalculation = false);
void setCompNBCoor(CoinPackedVector& vec, double& objViolation,
    const VPCParameters& params, const double* const currColValue,
    const double* const currSlackValue,
    const OsiSolverInterface* const origSolver,
    const std::vector<int>& nonBasicVarIndex,
    const std::vector<double>& nonBasicReducedCost, const int deletedVar = -1);
