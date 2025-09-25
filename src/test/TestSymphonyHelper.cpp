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


// --------------------- test current behavior remains -------------------------
TEST_CASE("Test doBranchAndBoundWithSymphony", "[SymphonyHelper::doBranchAndBoundWithSymphony]") {

  // parameters
  VPCParametersNamespace::VPCParameters vpc_params;
  vpc_params.set(VPCParametersNamespace::DISJ_TERMS, 64);
  vpc_params.set(VPCParametersNamespace::MODE, 0);  // partial BB tree
  vpc_params.set(VPCParametersNamespace::PARTIAL_BB_KEEP_PRUNED_NODES, 1);
  vpc_params.set(VPCParametersNamespace::SOLFILE, "../test/bm23.sol");

  // set parameters to use provided bound and skip heuristics
  vpc_params.set(BB_STRATEGY, get_bb_option_value({
      BB_Strategy_Options::user_cuts, // to allow VPCs and data collection
      BB_Strategy_Options::presolve_off, // instances will be presolved already
      BB_Strategy_Options::heuristics_off,  // already providing bound
      BB_Strategy_Options::use_best_bound,  // use provided solution
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
    solver->applyCuts(vpcs);  // cuts have to be added to model else Symphony will quit

    // solve with symphony
    BBInfo info;
    std::shared_ptr<OsiSymSolverInterface> parametric_model = std::shared_ptr<OsiSymSolverInterface>();
    doBranchAndBoundWithSymphony(vpc_params, vpc_params.get(BB_STRATEGY), solver,
                                 info, &vpcs, parametric_model);


    // check that we're optimal
    REQUIRE(info.obj == info.bound);
    REQUIRE(info.bound == 34);
    // should have done some branching with a few LP iterations each
    REQUIRE(0 < info.nodes);
    REQUIRE(info.nodes < info.iters);
    // dual bound monotonically improves
    REQUIRE(si.getObjValue() < info.last_cut_pass);
    REQUIRE(info.last_cut_pass < info.bound);
    // time increases monotonically
    REQUIRE(0 < info.time);
  }

  SECTION( "Test objective perturbed warm-start" ) {

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

    // LP iterations, and time should be less with warm start
    REQUIRE(info_ws.iters < info.iters);
    REQUIRE(info_ws.time < info.time);

    // final bounds should be the same
    REQUIRE(info_ws.bound == info_ws.obj);
    REQUIRE(info_ws.bound == info.bound);
    REQUIRE(info_ws.obj == info.obj);

  }

  SECTION( "Test rhs-perturbed warm-start" ){

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

    // LP iterations, and time should be less with warm start
    REQUIRE(info_ws.iters < info.iters);
    REQUIRE(info_ws.time < info.time);

    // final bounds should be the same
    REQUIRE(info_ws.bound == info_ws.obj);
    REQUIRE(info_ws.bound == info.bound);
    REQUIRE(info_ws.obj == info.obj);
  }

}

#endif // USE_SYMPHONY

