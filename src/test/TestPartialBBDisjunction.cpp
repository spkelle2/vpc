/**
 * @file TestVPCEventHandler.cpp
 * @author Sean Kelley
 * @date 2023-09-20
 */

#define CATCH_CONFIG_MAIN

// standard library
#include <cstdlib> // abs
#include <memory>
#include <vector> // vector

// unit test library
#include "catch.hpp"

// coin-or modules
#include <CoinWarmStartBasis.hpp>
#include "OsiClpSolverInterface.hpp" // OsiClpSolverInterface
#include "OsiCuts.hpp" // OsiCuts

// project modules
#include "CglVPC.hpp" // CglVPC
#include "Disjunction.hpp" // DisjExitReason
#include "PartialBBDisjunction.hpp" // PartialBBDisjunction
#include "SolverInterface.hpp" // SolverInterface

bool sameBasis(CoinWarmStart* basis1, CoinWarmStart* basis2) {
  // Check for null pointers
  if (basis1 == nullptr || basis2 == nullptr) {
    return basis1 == basis2;  // Both should be null to be considered equal
  }

  // Attempt to cast to CoinWarmStartBasis*
  auto* basis1Cast = dynamic_cast<CoinWarmStartBasis*>(basis1);
  auto* basis2Cast = dynamic_cast<CoinWarmStartBasis*>(basis2);

  // If either cast fails, the bases are not comparable
  if (basis1Cast == nullptr || basis2Cast == nullptr) {
    return false;
  }

  // Check structural and artificial counts
  if (basis1Cast->getNumStructural() != basis2Cast->getNumStructural() ||
      basis1Cast->getNumArtificial() != basis2Cast->getNumArtificial()) {
    return false;
  }

  // Compare structural statuses
  for (int i = 0; i < basis1Cast->getNumStructural(); i++) {
    if (basis1Cast->getStructStatus(i) != basis2Cast->getStructStatus(i)) {
      return false;
    }
  }

  // Compare artificial statuses
  for (int i = 0; i < basis1Cast->getNumArtificial(); i++) {
    if (basis1Cast->getArtifStatus(i) != basis2Cast->getArtifStatus(i)) {
      return false;
    }
  }

  return true;  // All checks passed
}

CoinWarmStart* createWarmStartAllArtificialActive(int numStructural, int numArtificial) {
  // Create a CoinWarmStartBasis instance
  CoinWarmStartBasis* warmStart = new CoinWarmStartBasis();

  // Resize to accommodate the structural and artificial variables
  warmStart->setSize(numStructural, numArtificial);

  // Set all structural variables to their bounds
  // Typically, at bounds means "at lower bound" or "at upper bound"
  for (int i = 0; i < numStructural; i++) {
    warmStart->setStructStatus(i, CoinWarmStartBasis::atUpperBound);  // or atLowerBound
  }

  // Set all artificial variables as "basic"
  for (int i = 0; i < numArtificial; i++) {
    warmStart->setArtifStatus(i, CoinWarmStartBasis::basic);
  }

  return warmStart;  // Return the warm start basis
}

// --------------------- test current behavior remains -------------------------
TEST_CASE("Test saveInformation", "[VPCEventHandler::saveInformation]") {

  // parameters
  VPCParametersNamespace::VPCParameters vpc_params;
  vpc_params.set(VPCParametersNamespace::DISJ_TERMS, 4);
  vpc_params.set(VPCParametersNamespace::MODE, 0);  // partial BB tree
  vpc_params.set(VPCParametersNamespace::PARTIAL_BB_KEEP_PRUNED_NODES, 1);

  // solver
  OsiClpSolverInterface si;
  SolverInterface* solver;
  si.readMps("../test/bm23.mps");
  si.initialSolve();
  std::vector<double> sol(si.getColSolution(), si.getColSolution() + si.getNumCols());
  solver = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(&si));

  SECTION( "Test without strong branching or pruned terms" ) {

    // make disjunction
    PartialBBDisjunction disj = PartialBBDisjunction(vpc_params);
    disj.prepareDisjunction(solver);

    // check to make sure warm start works
    for (int term_idx = 0; term_idx < 4; term_idx++){
      OsiSolverInterface* termSolver;
      disj.getSolverForTerm(termSolver, term_idx, solver, true, .001, NULL, false);
      // see how many iterations it takes to resolve from scratch
      CoinWarmStart* origin = createWarmStartAllArtificialActive(
          solver->getNumCols(), solver->getNumRows());
      termSolver->setWarmStart(origin);
      termSolver->resolve();
      REQUIRE(termSolver->getIterationCount() > 0);
      // see how many iterations it takes to resolve with the warm start
      termSolver->setWarmStart(disj.terms[term_idx].basis);
      termSolver->resolve();
      REQUIRE(termSolver->getIterationCount() == 0);

      OsiSolverInterface* termSolverExt;
      disj.getSolverForTerm(termSolverExt, term_idx, solver, true, .001, NULL, false);
      // see how many iterations it takes to resolve from scratch
      CoinWarmStart* origin_ext = createWarmStartAllArtificialActive(
          termSolverExt->getNumCols(), termSolverExt->getNumRows());
      termSolverExt->setWarmStart(origin_ext);
      termSolverExt->resolve();
      REQUIRE(termSolverExt->getIterationCount() > 0);
      // see how many iterations it takes to resolve with the warm start
      termSolverExt->setWarmStart(disj.terms[term_idx].basis_extended);
      termSolverExt->resolve();
      REQUIRE(termSolverExt->getIterationCount() == 0);
    }

    // perturb the problem
    for (int col_idx=0; col_idx < si.getNumCols(); col_idx++){
      bool perturb = col_idx == 25 || col_idx == 0 || col_idx == 13 || col_idx == 22;
      si.setObjCoeff(
          col_idx, perturb ? si.getObjCoefficients()[col_idx] * 5 :
          si.getObjCoefficients()[col_idx]);
    }

    // update the constraint bounds
    for (int row_idx=0; row_idx < si.getNumRows(); row_idx++){
      double b = si.getRowUpper()[row_idx];
      double looser_b = b > 0 ? b * 1.5 : b / 1.5;
      bool perturb = row_idx == 2 || row_idx == 7 || row_idx == 12 || row_idx == 17;
      si.setRowUpper(row_idx, perturb ? looser_b: b);
    }
    si.resolve();

    std::vector<std::unique_ptr<OsiSolverInterface>> term_solvers;
    PartialBBDisjunction param_disj = disj.parameterize(&si, &term_solvers);

    // check that we have the right metadata
    REQUIRE(param_disj.num_terms == 4);
    REQUIRE(param_disj.terms.size() == 4);
    REQUIRE(param_disj.common_changed_var.size() == 0);
    REQUIRE(param_disj.common_changed_bound.size() == 0);
    REQUIRE(param_disj.common_changed_value.size() == 0);
    REQUIRE(param_disj.common_ineqs.size() == 0);
    REQUIRE(param_disj.integer_sol.size() == 0);
    REQUIRE(param_disj.integer_obj > 1e300);
    REQUIRE(param_disj.root_obj == si.getObjValue());
    REQUIRE(isVal(param_disj.best_obj, 19.33, .1));
    REQUIRE(isVal(param_disj.worst_obj, 49.17, .1));

    for (int term_idx = 0; term_idx < 4; term_idx++){

      // these attributes should be the same
      REQUIRE(param_disj.terms[term_idx].type == disj.terms[term_idx].type);
      REQUIRE(param_disj.terms[term_idx].changed_var == disj.terms[term_idx].changed_var);
      REQUIRE(param_disj.terms[term_idx].changed_bound == disj.terms[term_idx].changed_bound);
      REQUIRE(param_disj.terms[term_idx].changed_value == disj.terms[term_idx].changed_value);
      REQUIRE(param_disj.terms[term_idx].is_feasible);
      // we don't ever get solver for term with disj constraints as tightened bounds so warm start stays same
      REQUIRE(sameBasis(param_disj.terms[term_idx].basis, disj.terms[term_idx].basis));

      // these attributes should be different
      REQUIRE(param_disj.terms[term_idx].obj != disj.terms[term_idx].obj);
      // warm start could be different but what's important is that it's a different object
      REQUIRE(param_disj.terms[term_idx].basis_extended != disj.terms[term_idx].basis_extended);

      // make sure the warm start works
      OsiSolverInterface* termSolverExt;
      param_disj.getSolverForTerm(termSolverExt, term_idx, solver, false, .001, NULL, true);

      // see how many iterations it takes to resolve from scratch
      CoinWarmStart* origin_ext = createWarmStartAllArtificialActive(
          termSolverExt->getNumCols(), termSolverExt->getNumRows());
      termSolverExt->setWarmStart(origin_ext);
      termSolverExt->resolve();
      REQUIRE(termSolverExt->getIterationCount() > 0);

      // see how many iterations it takes to resolve with the warm start
      termSolverExt->setWarmStart(param_disj.terms[term_idx].basis_extended);
      termSolverExt->resolve();
      REQUIRE(termSolverExt->getIterationCount() == 0);

      // check the cached solver matches the one used to create the disjunctive term
      OsiSolverInterface* cached_term_solver = term_solvers[term_idx].get();
      REQUIRE(sameBasis(cached_term_solver->getWarmStart(), termSolverExt->getWarmStart()));
    }
  }
}


