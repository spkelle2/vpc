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


using namespace VPCParametersNamespace;

// COIN-OR
#include <CoinTime.hpp>
#include <OsiCuts.hpp>

#ifdef USE_SYMPHONY

// set requested parameters for Symphony
void setStrategyForBBTestSymphony(const VPCParameters& params, const int strategy,
                                  OsiSymSolverInterface& model) {

  // get pointer to sym environment for later use
  sym_environment * env = model.getSymphonyEnvironment();

  // set parameters
  int seed = params.get(intParam::RANDOM_SEED);
  if (seed >= 0) model.setSymParam("random_seed", seed);

  // ---- Always-set basics ----
  model.setSymParam("time_limit", params.get(doubleParam::BB_TIMELIMIT));
  model.setSymParam("gap_limit", .01);  // .01% gap limit
  model.setSymParam("verbosity", params.get(VERBOSITY));
  model.setSymParam("keep_warm_start", true); // always keep the warm start tree

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
  }

  // Enable full strong Strong branching
  if (use_bb_option(std::abs(strategy), BB_Strategy_Options::strong_branching_on)) {
    model.setSymParam("max_presolve_iter", 1e6);  // high enough to not limit strong branching iterations
    model.setSymParam("limit_strong_branching_time", false);
  }

  for (int i = 0; i < env->mip->change_num; i++) {
    int change_type = env->mip->change_type[i];
    if (change_type == RHS_CHANGED) {
      // have to turn off cut generation under RHS changes. this should not be necessary
      // in our case if Symphony cleaned out all cuts generated from previous solve
      model.setSymParam("generate_cgl_cuts", false);

      // check to make sure we don't have a conflicting strategy
      verify(strategy <= 0 || use_bb_option(strategy, BB_Strategy_Options::all_cuts_off),
             "Symphony requires cuts turned off for disjunctive warm-starts under RHS changes");

    } else if (change_type == OBJ_COEFF_CHANGED) {
      // have to turn off reduced cost fixing to reuse disjunctions under objective changes
      model.setSymParam("do_reduced_cost_fixing", false);
    }
  }
}

// get a shared_ptr to a CoinWarmStart from the raw pointer in OsiSymSolverInterface
std::shared_ptr<CoinWarmStart> getWarmStartShared(OsiSymSolverInterface& model) {
    CoinWarmStart* ws_raw = model.getWarmStart();
    // set a custom deleter to safely delete the derived class
    return std::shared_ptr<CoinWarmStart>(ws_raw, [](CoinWarmStart* ptr) {
        delete ptr; // safely deletes derived class when shared_ptr is done
    });
}

// todo update to use SolverInterface instead of OsiSolverInterface
std::shared_ptr<CoinWarmStart> doBranchAndBoundWithSymphony(
    const VPCParameters& params, int strategy, const OsiSolverInterface* const si,
    BBInfo& info, const OsiCuts* cuts, const CoinWarmStart* ws,
    const OsiSolverInterface* const si_init, const OsiCuts* cuts_init) {

  // do some sanity checks on the initial solver
  verify((si_init == nullptr) == (ws == nullptr),
         "warm start and initial solver must be both provided or both null");
  if (ws) {
    // verify that si and si_init have the same constraint coefficient matrix and dimensions
    // same constraint matrix isn't technically required, but Symphony expects it
    verify(sameCoefficientMatrix(si, si_init),
           "initial solver must have same constraint matrix (and dimensions) as main solver");
    for (int i = 0; i < si->getNumCols(); i++) {
      verify(si_init->isInteger(i) == si->isInteger(i),
             "initial solver must have same integer variables as main solver");
    }
  }

  // copy the warm start so we don't corrupt the original in case it gets reused
  CoinWarmStart* ws_tmp = nullptr;
  if (ws){
    ws_tmp = dynamic_cast<CoinWarmStart*>(ws->clone());
  }

  // create solution and pool containers in case needed later
  // this needs to happen before we create the model so they still exist when
  // model destructor cleans them up
  int * xind = (int*)malloc(si->getNumCols() * sizeof(int));
  double * xval = (double*)malloc(si->getNumCols() * sizeof(double));
  sp_solution* sol = (sp_solution*)malloc(sizeof(sp_solution));
  sp_solution** solutions = (sp_solution**)malloc(sizeof(sp_solution*));
  sp_desc* pool = (sp_desc*)malloc(sizeof(sp_desc));

  // create models and get pointers to their environments
  OsiSymSolverInterface model;
  sym_environment * env = model.getSymphonyEnvironment();
  OsiSymSolverInterface model_init;

  // Copy the OsiSolverInterface into a/both SYMPHONY OSI solver/s and add cuts
  std::string f_name;
  if (si_init) {
    // if we have si_init, base model off of it and then manually modify to match si
    createTmpFileCopy(params, si_init, f_name);
    model_init.readMps(f_name.c_str());
    if (cuts_init && cuts_init->sizeCuts() > 0) {
      model_init.applyCuts(*cuts_init);
    }
  } else {
    // otherwise just base model off of si from the get go
    createTmpFileCopy(params, si, f_name);
  }
  model.readMps(f_name.c_str());
  if (cuts && cuts->sizeCuts() > 0) {
    model.applyCuts(*cuts);
  }

  if (si_init){
    // this would be a good place to check that the warm start certifies optimality for the initial problem
    // however, symphony does not seem to be able to do this properly
    // so we will just assume the user knows what they are doing

    // modify model to match its SolverInterface
    // need to do it this way so symphony doesn't get skeeved and ignore the warm start
    for (int j = 0; j < si->getNumCols(); j++) {
      if (!isVal(si->getColLower()[j], si_init->getColLower()[j])){
        model.setColLower(j, si->getColLower()[j]);
      }
      if (!isVal(si->getColUpper()[j], si_init->getColUpper()[j])){
        model.setColUpper(j, si->getColUpper()[j]);
      }
      if (!isVal(si->getObjCoefficients()[j], si_init->getObjCoefficients()[j])){
        model.setObjCoeff(j, si->getObjCoefficients()[j]);
      }
    }
    for (int i = 0; i < si->getNumRows(); i++) {
      if (!isVal(si->getRowLower()[i], si_init->getRowLower()[i])){
        model.setRowLower(i, si->getRowLower()[i]);
      }
      if (!isVal(si->getRowUpper()[i], si_init->getRowUpper()[i])){
        model.setRowUpper(i, si->getRowUpper()[i]);
      }
    }
  }

  // set strategy parameters
  setStrategyForBBTestSymphony(params, strategy, model);

  // set primal warm start if requested
  if ((strategy > 0) && use_bb_option(strategy, BB_Strategy_Options::use_best_bound)) {

    // Check if user provides a solution file
    std::string solfile = params.get(stringParam::SOLFILE);
    verify((solfile.size() > 4) && (solfile.compare(solfile.size() - 4, 4, ".sol") == 0),
           "VPC requires a .sol file to primal warm-start Symphony");

    // read in the solution file
    std::vector<double> vals;
    getSolFromFile(solfile.c_str(), vals);

    // check that the solution is valid
    verify(vals.size() == model.getNumCols(), "solution has wrong dimension");
    verify(isFeasible(model, vals, false), "solution is not feasible");

    // copy it over to C-style arrays for Symphony
    for (int i = 0; i < vals.size(); i++) {
      xind[i] = i;
      xval[i] = vals[i];
    }

    // get the objective value of the solution
    double obj_value = std::inner_product(vals.begin(), vals.end(), model.getObjCoefficients(), 0.0);

    // put the solution into a Symphony solution structure
    sol->objval = obj_value;
    sol->xlength = model.getNumCols();
    sol->xind = xind;
    sol->xval = xval;
    sol->node_index = 0;
    sol->node_level = 0;

    // create a solution pool from solution
    pool->max_solutions = 1;  // only keep best solution in the pool - we don't really care
    pool->num_solutions = 1;
    pool->total_num_sols_found = 1;
    solutions[0] = sol;
    pool->solutions = solutions;

    // set the primal warm start
    env->sp = pool;
  }

  // set dual warm start if requested
  model.setWarmStart(ws_tmp);

  // solve root node
  model.setSymParam("node_limit", 1);  // getLowerBound
  if (ws_tmp){
    model.resolve();
  } else {
    model.initialSolve();
  }

  // collect root statistics
  info.root_time = env->comp_times.readtime + env->comp_times.ub_overhead +
      env->comp_times.ub_heurtime + env->comp_times.lb_overhead +
      env->comp_times.lb_heurtime + env->tm->comp_times.communication +
      env->tm->comp_times.lp + env->tm->comp_times.lp_setup +
      env->tm->comp_times.separation + env->tm->comp_times.fixing +
      env->tm->comp_times.pricing + env->tm->comp_times.strong_branching +
      env->tm->comp_times.cut_pool + env->tm->comp_times.primal_heur;
  info.last_cut_pass = env->tm->lb;

  // solve the branch-and-bound tree
  model.setSymParam("node_limit", -1);
  model.resolve();

  // bounds
  info.bound = env->tm->lb;
  info.obj = env->tm->ub;

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

  // return the warm start in case we'd like to use it again
  std::shared_ptr<CoinWarmStart> ws_new = getWarmStartShared(model);

  // remove temporary files from createTmpFileCopy
  std::string f_name_no_ext = f_name.substr(0, f_name.size() - 4);
  std::string f_name_gz = f_name + ".gz";
  remove(f_name.c_str());
  remove(f_name_gz.c_str());
  remove(f_name_no_ext.c_str());

  return ws_new;
}

#endif /* USE_SYMPHONY */
