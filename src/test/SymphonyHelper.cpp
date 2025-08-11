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


// set requested parameters for Symphony
void setStrategyForBBTestSymphony(const VPCParameters& params,
                                  const int strategy,
                                  OsiSymSolverInterface& model,
                                  const double best_bound,
                                  int seed /* = -1 */) {

  if (seed < 0) seed = params.get(intParam::RANDOM_SEED);
  if (seed >= 0) model.setSymParam("random_seed", seed);

  // ---- Always-set basics ----
  model.setSymParam("time_limit", params.get(doubleParam::BB_TIMELIMIT));
  model.setSymParam("gap_limit", .01);  // .01% gap limit
  model.setSymParam("verbosity", params.get(VERBOSITY));

  // (new) enable dual warm-starts
  // todo turn on logging and log_interval to save tree and warm_start to read it back in with tree_log_file_name and cut_log_file_name
  // todo save pruned nodes with keep_description_of_pruned in pruned_node_file_name
  // todo save cuts with cp_warm_start in cp_warm_start_file_name

  // ---- Strategy-controlled toggles ----
  if (strategy > 0) {
    if (use_bb_option(strategy, BB_Strategy_Options::user_cuts)) {
      // print warning if the user provides cuts
        std::cerr << "User provided cuts, but this feature has not yet been "
                     "integrated with Symphony. Ignoring them and proceeding "
                     "with branch-and-bound without warm-started cut pool." << std::endl;
    }

    if (use_bb_option(strategy, BB_Strategy_Options::all_cuts_off)) {
      model.setSymParam("generate_cgl_cuts", false);
    }

    if (use_bb_option(strategy, BB_Strategy_Options::presolve_off)) {
      model.setSymParam("prep_level", 0);
    }

    if (use_bb_option(strategy, BB_Strategy_Options::heuristics_off)) {
      model.setSymParam("do_primal_heuristic", false);
    }

    if (use_bb_option(strategy, BB_Strategy_Options::use_best_bound)) {
      if (!isInfinity(std::abs(best_bound))) {
        //  - "upper_bound": prune nodes with obj >= UB (cutoff)
        // todo check to make sure this works
        model.setSymParam("upper_bound", best_bound * (1 + 1e-4));
      }
      // Check if user provides mip start or solution file
      std::string solfile = params.get(stringParam::SOLFILE);
      std::string ext1 = "_gurobi.sol.gz";
      std::string ext2 = "_gurobi.sol";
      std::string ext3 = "_gurobi.mst.gz";
      std::string ext4 = "_gurobi.mst";
      bool user_provides_start = false;
      user_provides_start |= (solfile.size() > ext1.size()) && (solfile.compare(solfile.size() - ext1.size(), ext1.size(), ext1) == 0);
      user_provides_start |= (solfile.size() > ext2.size()) && (solfile.compare(solfile.size() - ext2.size(), ext2.size(), ext2) == 0);
      user_provides_start |= (solfile.size() > ext3.size()) && (solfile.compare(solfile.size() - ext3.size(), ext3.size(), ext3) == 0);
      user_provides_start |= (solfile.size() > ext4.size()) && (solfile.compare(solfile.size() - ext4.size(), ext4.size(), ext4) == 0);
      if (user_provides_start) {
        // print warning if the user provides a solution file
        std::cerr << "User provided solution file" << solfile <<
          ". This is not supported with SYMPHONY. Ignoring it and proceeding "
          "with branch-and-bound without a MIP start." << std::endl;
      }
    }
  }

  // Enable full strong Strong branching
  if (use_bb_option(std::abs(strategy), BB_Strategy_Options::strong_branching_on)) {
    model.setSymParam("max_presolve_iter", 1e6);  // high enough to not limit strong branching iterations
    model.setSymParam("limit_strong_branching_time", false);
  }
}

void doBranchAndBoundWithSymphony(const VPCParameters& params, int strategy,
                                  const OsiSolverInterface* const solver,
                                  BBInfo& info, const double best_bound) {

  // Copy the OsiSolverInterface into a SYMPHONY OSI solver
  OsiSymSolverInterface model;
  std::string f_name;
  createTmpFileCopy(params, solver, f_name);
  model.readMps(f_name.c_str());

  // remove temporary files from createTmpFileCopy
  std::string f_name_no_ext = f_name.substr(0, f_name.size() - 4);
  std::string f_name_gz = f_name + ".gz";
  remove(f_name.c_str());
  remove(f_name_gz.c_str());
  remove(f_name_no_ext.c_str());

  // set parameters
  setStrategyForBBTestSymphony(params, strategy, model, best_bound);

  // todo: enable call back for data collection

  // Solve the problem
  model.initialSolve();
  model.branchAndBound();

  // Report solution
  if (model.isProvenOptimal()) {
    std::cout << "SYMPHONY found an optimal solution.\n";
    std::cout << "Objective value: " << model.getObjValue() << std::endl;
  } else {
    std::cout << "Solver stopped without finding optimal solution." << std::endl;
  }
}

#endif /* USE_SYMPHONY */
