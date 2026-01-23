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

typedef struct NODE_TIMES node_times;

std::shared_ptr<CoinWarmStart> getWarmStartShared(OsiSymSolverInterface& model);

// set requested parameters for Symphony
void setStrategyForBBTestSymphony(const VPCParametersNamespace::VPCParameters& params,
                                  const int strategy, OsiSymSolverInterface& model);

// helper function to solve a MILP with Symphony
void doBranchAndBoundWithSymphony(
    const VPCParametersNamespace::VPCParameters& params, int strategy,
    const OsiSolverInterface* const si, BBInfo& info, const OsiCuts* cuts,
    std::shared_ptr<OsiSymSolverInterface>& parametric_model, node_times* times=nullptr);
#endif /* USE_SYMPHONY */
