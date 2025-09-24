/**
 * @file SymphonyHelper.cpp
 * @author Shannon Kelley
 * @date 2025-Aug-04
 */

#ifdef USE_SYMPHONY

#include <stdio.h>
#include "sym_tm.h"
#include "sym_master.h"
#include <numeric> // for std::iota, std::inner_product, std::abs

#endif

#include "SymphonyHelper.hpp"

// Project files
#include "BBHelper.hpp"
#include "CutHelper.hpp" // applyCuts
#include "SolverHelper.hpp"
#include "VPCParameters.hpp"
#include "utility.hpp"

using namespace VPCParametersNamespace;

// COIN-OR
#include <CoinTime.hpp>
#include <OsiCuts.hpp>

#ifdef USE_SYMPHONY

// set requested parameters for Symphony
void setStrategyForBBTestSymphony(const VPCParameters& params, const int strategy,
                                  std::shared_ptr<OsiSymSolverInterface> model) {

  // get pointer to sym environment for later use
  sym_environment * env = model->getSymphonyEnvironment();

  // set parameters
  int seed = params.get(intParam::RANDOM_SEED);
  if (seed >= 0) model->setSymParam("random_seed", seed);

  // ---- Always-set basics ----
  model->setSymParam("time_limit", params.get(doubleParam::BB_TIMELIMIT));
  model->setSymParam("gap_limit", .01);  // .01% gap limit
  model->setSymParam("verbosity", params.get(VERBOSITY));
  model->setSymParam("keep_warm_start", true); // always keep the warm start tree

  // ---- Strategy-controlled toggles ----
  if (strategy > 0) {
    if (use_bb_option(strategy, BB_Strategy_Options::all_cuts_off)) {
      model->setSymParam("generate_cgl_cuts", false);
    }

    if (use_bb_option(strategy, BB_Strategy_Options::presolve_off)) {
      model->setSymParam("prep_level", 0);
    }

    if (use_bb_option(strategy, BB_Strategy_Options::heuristics_off)) {
      model->setSymParam("do_primal_heuristic", false);
    }
  }

  // Enable full strong Strong branching
  if (use_bb_option(std::abs(strategy), BB_Strategy_Options::strong_branching_on)) {
    model->setSymParam("max_presolve_iter", 1e6);  // high enough to not limit strong branching iterations
    model->setSymParam("limit_strong_branching_time", false);
  }

  for (int i = 0; i < env->mip->change_num; i++) {
    int change_type = env->mip->change_type[i];
    if (change_type == RHS_CHANGED) {
      // have to turn off cut generation under RHS changes. this should not be necessary
      // in our case if Symphony cleaned out all cuts generated from previous solve
      model->setSymParam("generate_cgl_cuts", false);

      // check to make sure we don't have a conflicting strategy
      verify(strategy <= 0 || use_bb_option(strategy, BB_Strategy_Options::all_cuts_off),
             "Symphony requires cuts turned off for disjunctive warm-starts under RHS changes");

    } else if (change_type == OBJ_COEFF_CHANGED) {
      // have to turn off reduced cost fixing to reuse disjunctions under objective changes
      model->setSymParam("do_reduced_cost_fixing", false);
    }
  }
  model->setSymParam("generate_cgl_cuts", false);
  model->setSymParam("do_reduced_cost_fixing", false);
}

// todo update to use SolverInterface instead of OsiSolverInterface
void doBranchAndBoundWithSymphony(
    const VPCParameters& params, int strategy, const OsiSolverInterface* const si,
    BBInfo& info, const OsiCuts* cuts, std::shared_ptr<OsiSymSolverInterface>& parametric_model) {

  // create a Symphony model for the input problem si
  std::shared_ptr<OsiSymSolverInterface> input_model = std::make_shared<OsiSymSolverInterface>();
  std::string f_name;
  createTmpFileCopy(params, si, f_name);
  input_model->readMps(f_name.c_str());
  if (cuts && cuts->sizeCuts() > 0) {
    input_model->applyCuts(*cuts);
  }

  // if provided a parametric model, modify it to match input_model
  if (parametric_model){
    // verify that both solvers have the same constraint coefficient matrix and dimensions
    // same constraint matrix isn't technically required, but Symphony expects it
    verify(sameCoefficientMatrix(input_model.get(), parametric_model.get()),
           "initial solver must have same constraint matrix (and dimensions) as main solver");
    for (int i = 0; i < input_model->getNumCols(); i++) {
      verify(input_model->isInteger(i) == parametric_model->isInteger(i),
             "initial solver must have same integer variables as main solver");
    }
    
    // modify parametric_model to match input_model
    // need to do it this way so symphony doesn't get skeeved and ignore the warm start
    for (int j = 0; j < input_model->getNumCols(); j++) {
      if (!isVal(input_model->getColLower()[j], parametric_model->getColLower()[j])){
        parametric_model->setColLower(j, input_model->getColLower()[j]);
      }
      if (!isVal(input_model->getColUpper()[j], parametric_model->getColUpper()[j])){
        parametric_model->setColUpper(j, input_model->getColUpper()[j]);
      }
      if (!isVal(input_model->getObjCoefficients()[j], parametric_model->getObjCoefficients()[j])){
        parametric_model->setObjCoeff(j, input_model->getObjCoefficients()[j]);
      }
    }
    for (int i = 0; i < input_model->getNumRows(); i++) {
      if (!isVal(input_model->getRowLower()[i], parametric_model->getRowLower()[i])){
        parametric_model->setRowLower(i, input_model->getRowLower()[i]);
      }
      if (!isVal(input_model->getRowUpper()[i], parametric_model->getRowUpper()[i])){
        parametric_model->setRowUpper(i, input_model->getRowUpper()[i]);
      }
    }
  } else {
    // otherwise this will be our base instance for the parametric model, so just copy it over
    parametric_model = input_model;
  }

  sym_environment * env = parametric_model->getSymphonyEnvironment();

  // set strategy parameters
  setStrategyForBBTestSymphony(params, strategy, parametric_model);

  // solve the branch-and-bound tree
  parametric_model->resolve();

  // bounds
  info.last_cut_pass = env->tm->stat.root_lb;
  info.bound = std::min(env->tm->lb, env->tm->ub); // sometimes symphony returns bad bounds
  info.obj = env->tm->ub; // ub should always match the best known integer solution

  // times
  info.time = env->comp_times.readtime + env->comp_times.ub_overhead +
      env->comp_times.ub_heurtime + env->comp_times.lb_overhead +
      env->comp_times.lb_heurtime + env->tm->comp_times.communication +
      env->tm->comp_times.lp + env->tm->comp_times.lp_setup +
      env->tm->comp_times.separation + env->tm->comp_times.fixing +
      env->tm->comp_times.pricing + env->tm->comp_times.strong_branching +
      env->tm->comp_times.cut_pool + env->tm->comp_times.primal_heur +
      info.root_time;  // symphony resets the tree manager on the resolve so add back in previous time

  // processing steps
  info.nodes = env->tm->stat.analyzed;
  info.iters = env->tm->lp_stat.lp_iter_num;

  // remove temporary files from createTmpFileCopy
  std::string f_name_no_ext = f_name.substr(0, f_name.size() - 4);
  std::string f_name_gz = f_name + ".gz";
  remove(f_name.c_str());
  remove(f_name_gz.c_str());
  remove(f_name_no_ext.c_str());
}

#endif /* USE_SYMPHONY */
