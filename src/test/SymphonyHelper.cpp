/**
 * @file SymphonyHelper.cpp
 * @author Shannon Kelley
 * @date 2025-Aug-04
 */

#ifdef USE_SYMPHONY

//// Symphony C headers FIRST
//extern "C" {
//
#include <stdio.h>   // defines FILE
//
//  // 1) Pull in TM first so we get the 4-arg prototypes.
#include "sym_tm.h"
//
//  // 2) Mask legacy prototypes ONLY for sym_master.h.
//  #define read_node  read_node__masked__do_not_use
//  #define write_node write_node__masked__do_not_use
//
//  // 3) Now include sym_master.h; any legacy 2-arg decls will be renamed,
//  //    avoiding a signature conflict with the already-seen 4-arg versions.
#include "sym_master.h"
//
//  // 4) Unmask so the rest of your TU sees normal names again (we're not calling them here).
//  #undef read_node
//  #undef write_node
//}

#endif

#include "SymphonyHelper.hpp"

// Project files
#include "BBHelper.hpp"
#include "CutHelper.hpp" // applyCuts
#include "SolverHelper.hpp"
#include "VPCParameters.hpp"
using namespace VPCParametersNamespace;

// COIN-OR
#include <CoinTime.hpp>
#include <OsiCuts.hpp>

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
        model.setSymParam("upper_bound", (1 + 1e-6) * best_bound); // set to granularity based on symphony make the same calculation
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
                                  BBInfo& info, const OsiCuts* cuts,
                                  const double best_bound) {

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

  // add user cuts
  if (cuts && cuts->sizeCuts() > 0) {
    model.applyCuts(*cuts);
  }

  // solve model
  model.initialSolve();

  // collect statistics
  sym_environment * env = model.getSymphonyEnvironment();

  // bounds
  info.last_cut_pass = env->tm->stat.root_lb;
  info.bound = env->tm->lb;
  info.obj = env->tm->ub;

  // times todo: get root processing time
  info.time = env->comp_times.readtime + env->comp_times.ub_overhead +
      env->comp_times.ub_heurtime + env->comp_times.lb_overhead +
      env->comp_times.lb_heurtime + env->tm->comp_times.communication +
      env->tm->comp_times.lp + env->tm->comp_times.lp_setup +
      env->tm->comp_times.separation + env->tm->comp_times.fixing +
      env->tm->comp_times.pricing + env->tm->comp_times.strong_branching +
      env->tm->comp_times.cut_pool + env->tm->comp_times.primal_heur;

  // processing steps
  info.nodes = env->tm->stat.analyzed;
  info.iters = env->tm->lp_stat.lp_iter_num;
}

#endif /* USE_SYMPHONY */
