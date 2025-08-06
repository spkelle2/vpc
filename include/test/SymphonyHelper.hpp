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
// helper function to solve a MILP with Symphony
void doBranchAndBoundWithSymphony(const VPCParametersNamespace::VPCParameters& params,
                                  int strategy, const OsiSolverInterface* const solver,
                                  BBInfo& info, const double best_bound);
#endif /* USE_SYMPHONY */
