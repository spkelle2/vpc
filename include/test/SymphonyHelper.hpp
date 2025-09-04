/**
 * @file SymphonyHelper.hpp
 * @author Shannon Kelley
 * @date 2025-Aug-04
 */
#pragma once

class OsiSolverInterface;
class OsiCuts;

struct BBInfo;
namespace VPCParametersNamespace {
  struct VPCParameters;
}

#ifdef USE_SYMPHONY

#include "OsiSymSolverInterface.hpp"
#include <memory> // shared_ptr

std::shared_ptr<CoinWarmStart> getWarmStartShared(OsiSymSolverInterface& model);

// set requested parameters for Symphony
void setStrategyForBBTestSymphony(const VPCParametersNamespace::VPCParameters& params,
                                  const int strategy, OsiSymSolverInterface& model);

// helper function to solve a MILP with Symphony
std::shared_ptr<CoinWarmStart> doBranchAndBoundWithSymphony(
    const VPCParametersNamespace::VPCParameters& params, int strategy,
    const OsiSolverInterface* const si, BBInfo& info,
    const OsiCuts* cuts = nullptr, const CoinWarmStart* ws = nullptr,
    const OsiSolverInterface* const si_init = nullptr,
    const OsiCuts* cuts_init = nullptr);
#endif /* USE_SYMPHONY */
