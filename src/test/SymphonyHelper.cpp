/**
 * @file SymphonyHelper.cpp
 * @author Shannon Kelley
 * @date 2025-Aug-04
 */

#include "SymphonyHelper.hpp"

// Project files
#include "BBHelper.hpp"
#include "CutHelper.hpp" // applyCuts
#include "SolverHelper.hpp"
#include "VPCParameters.hpp"
using namespace VPCParametersNamespace;

#ifdef USE_SYMPHONY

// Symphony modules
#include "OsiSymSolverInterface.hpp"

void doBranchAndBoundWithSymphony(const VPCParameters& params, int strategy,
                                  const OsiSolverInterface* const solver,
                                  BBInfo& info, const double best_bound) {

  // Copy the OsiSolverInterface into a SYMPHONY OSI solver
  OsiSymSolverInterface symphony_si;
  std::string f_name;
  createTmpFileCopy(params, solver, f_name);
  symphony_si.readMps(f_name.c_str());

  // remove temporary files from createTmpFileCopy
  std::string f_name_no_ext = f_name.substr(0, f_name.size() - 4);
  std::string f_name_gz = f_name + ".gz";
  remove(f_name.c_str());
  remove(f_name_gz.c_str());
  remove(f_name_no_ext.c_str());

  // Optionally: set parameters here
  // symphony_si.setSymParam("max_nodes", 100);  // example

  // Solve the problem
  symphony_si.initialSolve(); // if you want LP relaxation first
  symphony_si.branchAndBound();

  // Report solution
  if (symphony_si.isProvenOptimal()) {
    std::cout << "SYMPHONY found an optimal solution.\n";
    std::cout << "Objective value: " << symphony_si.getObjValue() << std::endl;
  } else {
    std::cout << "Solver stopped without finding optimal solution." << std::endl;
  }
}

#endif /* USE_SYMPHONY */
