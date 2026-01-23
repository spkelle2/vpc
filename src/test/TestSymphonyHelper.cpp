/**
 * @file TestSymphonyHelper.cpp
 * @author Shannon Kelley
 * @date 2025-08-27
 */

#define CATCH_CONFIG_MAIN

// standard library
#include <cstdlib> // abs
#include <memory>
#include <vector> // vector

// unit test library
#include "catch.hpp"

// coin-or modules
#include <CoinWarmStart.hpp>
#include "OsiCuts.hpp" // OsiCuts

// project modules
#include "VPCParameters.hpp"
#include "CglVPC.hpp"
#include "BBHelper.hpp"

using namespace VPCParametersNamespace;

#ifdef USE_SYMPHONY

#include "SymphonyHelper.hpp" // doBranchAndBoundWithSymphony
#include "sym_tm.h"
#include "sym_master.h"


// --------------------- test current behavior remains -------------------------
TEST_CASE("Test doBranchAndBoundWithSymphony", "[SymphonyHelper::doBranchAndBoundWithSymphony]") {

  // parameters
  VPCParametersNamespace::VPCParameters vpc_params;
  vpc_params.set(VPCParametersNamespace::DISJ_TERMS, 64);
  vpc_params.set(VPCParametersNamespace::MODE, 0);  // partial BB tree
  vpc_params.set(VPCParametersNamespace::PARTIAL_BB_KEEP_PRUNED_NODES, 1);
  vpc_params.set(VPCParametersNamespace::SOLFILE, "../test/bm23.sol");
  vpc_params.set(DISJUNCTION_SOLVER, "SYMPHONY");

  // set parameters to use provided bound and skip heuristics
  vpc_params.set(BB_STRATEGY, get_bb_option_value({
      BB_Strategy_Options::user_cuts, // to allow VPCs and data collection
      BB_Strategy_Options::presolve_off, // instances will be presolved already
      BB_Strategy_Options::heuristics_on,  // already providing bound
//      BB_Strategy_Options::use_best_bound,  // use provided solution
      BB_Strategy_Options::all_cuts_off // don't use any cuts other than VPCs
  }));

  SECTION( "Test first instance solve" ) {

    // solver
    OsiClpSolverInterface si;
    SolverInterface* solver;
    si.readMps("../test/bm23.mps");
    si.initialSolve();
    solver = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(&si));

    // make vpcs
    OsiCuts vpcs;
    CglVPC gen = CglVPC(vpc_params);
    gen.generateCuts(si, vpcs);

    // check when we use symphony to generate vpcs we see reasonable improvements
    si.applyCuts(vpcs);
    si.resolve();
    REQUIRE(si.getObjValue() > 25);  // should improve a lot over LP bound of 20.57
    REQUIRE(si.getObjValue() < 30);  // but not too much

    // solve with symphony
    BBInfo info;
    node_times times;
    std::shared_ptr<OsiSymSolverInterface> parametric_model = std::shared_ptr<OsiSymSolverInterface>();
    doBranchAndBoundWithSymphony(vpc_params, vpc_params.get(BB_STRATEGY), solver,
                                 info, &vpcs, parametric_model, &times);

    // check that we're optimal
    REQUIRE(info.obj == info.bound);
    REQUIRE(info.bound == 34);
    // should have done some branching with a few LP iterations each
    REQUIRE(0 < info.nodes);
    REQUIRE(info.nodes < info.iters);
    REQUIRE(0 < info.time);
    // dual bound monotonically improves
    REQUIRE(si.getObjValue() < info.bound);

    // times that should be nonzero
    REQUIRE(times.lp > 0.01);
    REQUIRE(times.strong_branching > 0.01);
    REQUIRE(times.primal_heur > 0.01);
  }

  SECTION( "Test objective perturbed warm-start" ) {

    // small changes should result in warm-start being effective

    // solver
    OsiClpSolverInterface si;
    SolverInterface* solver;
    si.readMps("../test/bm23.mps");
    si.initialSolve();
    solver = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(&si));

    // solve initial instance with symphony to get warm start
    vpc_params.set(VPCParametersNamespace::BB_NODE_LIMIT, 64);  // partial BB tree
    BBInfo info_initial;
    std::shared_ptr<OsiSymSolverInterface> parametric_model = std::shared_ptr<OsiSymSolverInterface>();
    doBranchAndBoundWithSymphony(vpc_params, vpc_params.get(BB_STRATEGY), solver,
                                 info_initial, nullptr, parametric_model);
    vpc_params.set(VPCParametersNamespace::BB_NODE_LIMIT, -1);  // partial BB tree

    // create a solver with a perturbed objective to force a different, but previously found solution
    OsiSolverInterface * si_ptb = si.clone();
    si_ptb->setObjCoeff(0, 4.);
    SolverInterface* solver_ptb;
    si_ptb->initialSolve();
    solver_ptb = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(si_ptb));

    // solve with warm start
    BBInfo info_ws;
    doBranchAndBoundWithSymphony(vpc_params, vpc_params.get(BB_STRATEGY), solver_ptb,
                                 info_ws, nullptr, parametric_model);

    // solve with no warm start
    BBInfo info;
    std::shared_ptr<OsiSymSolverInterface> dummy_model = std::shared_ptr<OsiSymSolverInterface>();
    doBranchAndBoundWithSymphony(vpc_params, vpc_params.get(BB_STRATEGY), solver_ptb,
                                 info, nullptr, dummy_model);

    // check monotonicity for warm start solve
    // time increases monotonically
    REQUIRE(0 < info_ws.time);
    // nodes increase monotonically
    REQUIRE(0 < info_ws.nodes);
    // iterations increase monotonically
    REQUIRE(0 < info_ws.iters);

    // this small of warm start should improve iterations but not time
    REQUIRE(info_ws.iters < info.iters);
    REQUIRE(info_ws.time < info.time);

    // final bounds should be the same
    REQUIRE(info_ws.bound == info_ws.obj);
    REQUIRE(info_ws.bound == info.bound);
    REQUIRE(info_ws.obj == info.obj);
  }

  SECTION( "Test rhs-perturbed warm-start" ){

    // small changes should result in warm-start being effective

    // solver
    OsiClpSolverInterface si;
    SolverInterface* solver;
    si.readMps("../test/bm23.mps");
    si.initialSolve();
    solver = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(&si));

    // solve initial instance with symphony to get warm start
    BBInfo info_initial;
    std::shared_ptr<OsiSymSolverInterface> parametric_model = std::shared_ptr<OsiSymSolverInterface>();
    doBranchAndBoundWithSymphony(vpc_params, vpc_params.get(BB_STRATEGY), solver,
                                 info_initial, nullptr, parametric_model);

    // create a solver with a perturbed RHS but same solution
    OsiSolverInterface * si_ptb = si.clone();
    si_ptb->setRowUpper(0, 60);
    si_ptb->initialSolve();
    SolverInterface* solver_ptb;
    solver_ptb = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(si_ptb));

    // solve with warm start
    BBInfo info_ws;
    doBranchAndBoundWithSymphony(vpc_params, vpc_params.get(BB_STRATEGY), solver_ptb,
                                 info_ws, nullptr, parametric_model);

    // solve with no warm start
    BBInfo info;
    std::shared_ptr<OsiSymSolverInterface> dummy_model = std::shared_ptr<OsiSymSolverInterface>();
    doBranchAndBoundWithSymphony(vpc_params, vpc_params.get(BB_STRATEGY), solver_ptb,
                                 info, nullptr, dummy_model);

    // check monotonicity for warm start solve
    // time increases monotonically
    REQUIRE(0 < info_ws.time);
    // nodes increase monotonically
    REQUIRE(0 < info_ws.nodes);
    // iterations increase monotonically
    REQUIRE(0 < info_ws.iters);

    // warm start should improve performance
    REQUIRE(info_ws.iters < info.iters);
    REQUIRE(info_ws.time < info.time);

    // final bounds should be the same
    REQUIRE(info_ws.bound == info_ws.obj);
    REQUIRE(info_ws.bound == info.bound);
    REQUIRE(info_ws.obj == info.obj);
  }

  SECTION( "Test warm-start lower bound for objective changes" ) {

    // update parameters to ditch provided bound since we perturb the objective significantly
//    vpc_params.set(BB_STRATEGY, get_bb_option_value({
//        BB_Strategy_Options::user_cuts, // to allow VPCs and data collection
//        BB_Strategy_Options::presolve_off, // instances will be presolved already
//        BB_Strategy_Options::heuristics_off,  // already providing bound
//        BB_Strategy_Options::all_cuts_off // don't use any cuts other than VPCs
//    }));

    // solver
    OsiClpSolverInterface si;
    SolverInterface* solver;
    si.readMps("../test/bm23.mps");
    si.initialSolve();
    solver = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(&si));

    // solve initial instance with symphony to get warm start
    BBInfo info_initial;
    std::shared_ptr<OsiSymSolverInterface> parametric_model = std::shared_ptr<OsiSymSolverInterface>();
    doBranchAndBoundWithSymphony(vpc_params, vpc_params.get(BB_STRATEGY), solver,
                                 info_initial, nullptr, parametric_model);

    // create a solver with a perturbed objective to force a different, but previously found solution
    OsiSolverInterface * si_ptb = si.clone();
    for (int i = 0; i < si_ptb->getNumCols(); i++){
      si_ptb->setObjCoeff(i, -1 * si_ptb->getObjCoefficients()[i]);
    }
    SolverInterface* solver_ptb;
    si_ptb->initialSolve();
    solver_ptb = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(si_ptb));

    // solve with warm start
    BBInfo info_ws;
    doBranchAndBoundWithSymphony(vpc_params, vpc_params.get(BB_STRATEGY), solver_ptb,
                                 info_ws, nullptr, parametric_model);

    // check that we're optimal
    REQUIRE(info_ws.obj == info_ws.bound);
    REQUIRE(info_ws.bound == -93);

    // dual bound and time monotonically improve
    REQUIRE(si_ptb->getObjValue() < info_ws.bound);

    // check monotonicity for warm start solve
    // time increases monotonically
    REQUIRE(0 < info_ws.time);
    // nodes increase monotonically
    REQUIRE(0 < info_ws.nodes);
    // iterations increase monotonically
    REQUIRE(0 < info_ws.iters);
  }

  SECTION( "Test warm-start lower bound for RHS changes" ) {

    vpc_params.set(VPCParametersNamespace::SOLFILE, "../test/bm23_rhs.sol");

    // solver
    OsiClpSolverInterface si;
    SolverInterface* solver;
    si.readMps("../test/bm23.mps");
    si.initialSolve();
    solver = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(&si));

    // solve initial instance with symphony to get warm start
    BBInfo info_initial;
    std::shared_ptr<OsiSymSolverInterface> parametric_model = std::shared_ptr<OsiSymSolverInterface>();
    doBranchAndBoundWithSymphony(vpc_params, vpc_params.get(BB_STRATEGY), solver,
                                 info_initial, nullptr, parametric_model);

    // create a solver with a perturbed objective to force a different, but previously found solution
    OsiSolverInterface * si_ptb = si.clone();
    std::vector<int> constraint_idxs = {4, 8, 15};
    for (int i = 0; i < si_ptb->getNumRows(); i++){
      if (std::find(constraint_idxs.begin(), constraint_idxs.end(), i) != constraint_idxs.end()){
        si_ptb->setRowUpper(i, si_ptb->getRowUpper()[i] - 4);
      } else {
        si_ptb->setRowUpper(i, si_ptb->getRowUpper()[i] - 7);
      }
    }
    SolverInterface* solver_ptb;
    si_ptb->initialSolve(); // 43.41 LP objective
    solver_ptb = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(si_ptb));

    // solve with warm start - disjunctive dual bound was 32.59 for initial problem
    BBInfo info_ws;
    doBranchAndBoundWithSymphony(vpc_params, vpc_params.get(BB_STRATEGY), solver_ptb,
                                 info_ws, nullptr, parametric_model);

    // check that we're optimal
    REQUIRE(info_ws.obj == info_ws.bound);
    REQUIRE(info_ws.bound == 63);

    // dual bound and time monotonically improve
    // main test here is just finding a gap at last cut pass (i.e. after root node)
    REQUIRE(si_ptb->getObjValue() < info_ws.bound);
    // check monotonicity for warm start solve
    // time increases monotonically
    REQUIRE(0 < info_ws.time);
    // nodes increase monotonically
    REQUIRE(0 < info_ws.nodes);
    // iterations increase monotonically
    REQUIRE(0 < info_ws.iters);
  }

}

#endif // USE_SYMPHONY

