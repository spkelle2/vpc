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
#include <cmath> // for std::isfinite

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
  if (seed >= 0) {
    model->setSymParam("random_seed", seed);
  }

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

// no getting around providing OsiSymSolverInterface, passing just a warm start does not work
void doBranchAndBoundWithSymphony(
    const VPCParameters& params, int strategy, const OsiSolverInterface* const si,
    BBInfo& info, const OsiCuts* cuts, std::shared_ptr<OsiSymSolverInterface>& parametric_model) {

  // create solution and pool containers in case needed later
  // this needs to happen before we create the model so they still exist when
  // model destructor cleans them up
  int * xind = (int*)malloc(si->getNumCols() * sizeof(int));
  double * xval = (double*)malloc(si->getNumCols() * sizeof(double));
  sp_solution* sol = (sp_solution*)malloc(sizeof(sp_solution));
  sp_solution** solutions = (sp_solution**)malloc(sizeof(sp_solution*));
  sp_desc* pool = (sp_desc*)malloc(sizeof(sp_desc));
  double obj_value = std::numeric_limits<double>::max();
  bool provide_sol = (strategy > 0) && use_bb_option(strategy, BB_Strategy_Options::use_best_bound);
  bool provide_parametric = (parametric_model != nullptr);

  // create a Symphony model for the input problem si
  std::shared_ptr<OsiSymSolverInterface> input_model = std::make_shared<OsiSymSolverInterface>();
  std::shared_ptr<OsiSymSolverInterface> root_model = std::make_shared<OsiSymSolverInterface>();
  std::string f_name;
  createTmpFileCopy(params, si, f_name);
  input_model->readMps(f_name.c_str());
  root_model->readMps(f_name.c_str());

  // input cuts in a standardized <= sense if provided
  OsiCuts cuts_std;
  if (cuts && cuts->sizeCuts() > 0) {
    // Flip sense of cuts to be <= so that Symphony doesn't complain about changes
    // Symphony internally stores all constraints as <=
    for (int i = 0; i < cuts->sizeRowCuts(); i++) {
      const OsiRowCut *oldCut = cuts->rowCutPtr(i);

      double lb = oldCut->lb();
      double ub = oldCut->ub();

      // If it has a finite lower bound, flip to <= form
      bool flipped = false;
      if (lb > -COIN_DBL_MAX) {
        ub = -lb;               // new upper bound
        lb = -COIN_DBL_MAX;     // no lower bound
        flipped = true;
      }

      // Copy row vector
      const CoinPackedVector &oldRow = oldCut->row();
      const int *indices   = oldRow.getIndices();
      const double *values = oldRow.getElements();
      int n = oldRow.getNumElements();

      CoinPackedVector newRow;
      for (int j = 0; j < n; ++j) {
        double val = values[j];
        if (flipped) {
          val = -val;  // flip coefficients
        }
        newRow.insert(indices[j], val);
      }

      // Construct new cut
      OsiRowCut newCut;
      newCut.setLb(lb);
      newCut.setUb(ub);
      newCut.setRow(newRow);

      cuts_std.insert(newCut);
    }

    input_model->applyCuts(cuts_std);
    root_model->applyCuts(cuts_std);
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

    // clear out old stats
    sym_environment * env = parametric_model->getSymphonyEnvironment();
    env->warm_start->lp_stat = lp_stat_desc();
    env->warm_start->stat = problem_stat();
    env->warm_start->comp_times = node_times();
  } else {
    // otherwise this will be our base instance for the parametric model, so just copy it over
    parametric_model = input_model;
  }

  // get pointer to sym environments for later use
  sym_environment * env = parametric_model->getSymphonyEnvironment();
  sym_environment * root_env = root_model->getSymphonyEnvironment();

  // set strategy parameters
  setStrategyForBBTestSymphony(params, strategy, parametric_model);
  setStrategyForBBTestSymphony(params, strategy, root_model);

  // copy over warm start from the parametric model to the root model
  if (env->warm_start){
    // remove any existing solution pool to make it fair vs cut generation only
    FREE(env->sp);
    env->warm_start->best_sol = lp_sol();

    // copy over the warm start structure
    root_env->warm_start = create_copy_warm_start(env->warm_start);
    root_env->warm_start->force_resolve_tree = true;  // resolve the tree to get the correct bound
    root_env->mip = create_copy_mip_desc(env->mip);
  }

  // resolve each node in the tree to get the bound
  root_model->setSymParam("node_limit", 1);
  root_model->resolve();
  info.last_cut_pass = root_env->tm->lb;

  // process just the first node (without resolving the tree to get the correct time info) to start
  parametric_model->setSymParam("node_limit", 1);
  parametric_model->resolve();

  // capture first node specific stats
  info.root_time = env->comp_times.readtime + env->comp_times.ub_overhead +
    env->comp_times.ub_heurtime + env->comp_times.lb_overhead +
    env->comp_times.lb_heurtime + env->tm->comp_times.communication +
    env->tm->comp_times.lp + env->tm->comp_times.lp_setup +
    env->tm->comp_times.separation + env->tm->comp_times.fixing +
    env->tm->comp_times.pricing + env->tm->comp_times.strong_branching +
    env->tm->comp_times.cut_pool + env->tm->comp_times.primal_heur;
  info.root_iters =  env->tm->lp_stat.lp_iter_num;
  info.root_passes = env->tm->stat.analyzed;  // represents nodes processed so far

  // set primal warm start if requested
  if (provide_sol) {

    // Check if user provides a solution file
    std::string solfile = params.get(stringParam::SOLFILE);
    verify((solfile.size() > 4) && (solfile.compare(solfile.size() - 4, 4, ".sol") == 0),
           "VPC requires a .sol file to primal warm-start Symphony");

    // read in the solution file
    std::shared_ptr<double> cached_obj = std::make_shared<double>(std::numeric_limits<double>::quiet_NaN());
    std::vector<double> vals;
    getSolFromFile(solfile.c_str(), vals, cached_obj.get());

    // check that the solution is valid
    verify(vals.size() == parametric_model->getNumCols(), "solution has wrong dimension");
    verify(isFeasible(*parametric_model.get(), vals, false), "solution is not feasible");

    // copy it over to C-style arrays for Symphony
    for (int i = 0; i < vals.size(); i++) {
      xind[i] = i;
      xval[i] = vals[i];
    }

    // get the objective value of the solution
    obj_value = std::inner_product(vals.begin(), vals.end(), parametric_model->getObjCoefficients(), 0.0);
    if (cached_obj && std::isfinite(*cached_obj)) {
      // sometimes when symphony warm-starts it doesn't seem to compute the objective value
      // correctly, so just use the cached value if available
      verify(std::abs(obj_value - *cached_obj) / std::abs(*cached_obj) <= 1e-5,
             "provided solution objective value does not match computed objective value");
    }

    // put the solution into a Symphony solution structure
    sol->objval = obj_value;
    sol->xlength = parametric_model->getNumCols();
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

    // add the solution pool to the symphony environment
    env->sp = pool;
  }

  // now solve the branch-and-bound tree to optimality
  parametric_model->setSymParam("node_limit", params.get(BB_NODE_LIMIT));
  parametric_model->resolve();

  // bounds
  info.bound = env->tm->lb;
  // if we don't have a primal bound to report, use the provided bound
  info.obj = env->tm->ub != 0 ? env->tm->ub : obj_value;

  // nodes - they're cumulative across resolves so no special handling
  info.nodes = env->tm->stat.analyzed;

  // iterations - also cumulative across reesolves
  info.iters = env->tm->lp_stat.lp_iter_num;

  // total time
  info.time = env->comp_times.readtime + env->comp_times.ub_overhead +
    env->comp_times.ub_heurtime + env->comp_times.lb_overhead +
    env->comp_times.lb_heurtime + env->tm->comp_times.communication +
    env->tm->comp_times.lp + env->tm->comp_times.lp_setup +
    env->tm->comp_times.separation + env->tm->comp_times.fixing +
    env->tm->comp_times.pricing + env->tm->comp_times.strong_branching +
    env->tm->comp_times.cut_pool + env->tm->comp_times.primal_heur;

  // remove temporary files from createTmpFileCopy
  std::string f_name_no_ext = f_name.substr(0, f_name.size() - 4);
  std::string f_name_gz = f_name + ".gz";
  remove(f_name.c_str());
  remove(f_name_gz.c_str());
  remove(f_name_no_ext.c_str());
}

#endif /* USE_SYMPHONY */
