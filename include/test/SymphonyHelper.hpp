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

// set requested parameters for Symphony
void setStrategyForBBTestSymphony(const VPCParametersNamespace::VPCParameters& params,
                                  const int strategy, OsiSymSolverInterface& model,
                                  const double best_bound, int seed = -1);

// helper function to solve a MILP with Symphony
void doBranchAndBoundWithSymphony(const VPCParametersNamespace::VPCParameters& params,
                                  int strategy, const OsiSolverInterface* const solver,
                                  BBInfo& info, const OsiCuts* cuts = nullptr,
                                  const double best_bound = std::numeric_limits<double>::max());
#endif /* USE_SYMPHONY */
