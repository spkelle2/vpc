/**
 * @file PartialBBDisjunction.cpp
 * @author A. M. Kazachkov
 * @date 2018-02-22
 */
#include "PartialBBDisjunction.hpp"

// std lib
#include <deque>

// COIN-OR
#ifdef USE_CBC
#include <CbcModel.hpp>

// Cbc
#include <CbcSolver.hpp>
#include <CbcTree.hpp>

// Variable selection
#include <CbcBranchDefaultDecision.hpp>
#include <CbcBranchStrongDecision.hpp>
#include <OsiChooseStrongCustom.hpp>
#include <CbcBranchDynamic.hpp>

// Node selection
#include <CbcCompareDefault.hpp>
#include <CbcCompareBFS.hpp>
#include <CbcCompareDepth.hpp>
#include <CbcCompareEstimate.hpp>
#include <CbcCompareObjective.hpp>

// General strategy
#include <CbcStrategy.hpp>
#endif // USE_CBC

// Project files
#include "CglVPC.hpp" // get timer information too
#include "SolverHelper.hpp"
#include "utility.hpp"
#include "VPCEventHandler.hpp" // currently requires Cbc
#include "BBHelper.hpp"
#include "SymphonyHelper.hpp"

// Symphony
#include "symphony.h"
#include "sym_types.h"

#ifdef TRACE
#include "vpc_debug.hpp"
#endif

using namespace VPCParametersNamespace;

//<---***************** PUBLIC  **********************-->

PartialBBDisjunction::PartialBBDisjunction(const VPCParametersNamespace::VPCParameters& param) {
  initialize(NULL, &param);
} /* param constructor */

PartialBBDisjunction::PartialBBDisjunction(const PartialBBDisjunction& source, const VPCParametersNamespace::VPCParameters& param) {
  initialize(&source, &param);
} /* copy & param constructor */

PartialBBDisjunction::PartialBBDisjunction() {
  initialize(NULL, NULL);
} /* default constructor */

PartialBBDisjunction::PartialBBDisjunction(const PartialBBDisjunction& source) {
  initialize(&source, NULL);
} /* copy constructor */

PartialBBDisjunction::~PartialBBDisjunction() {
} /* destructor */

PartialBBDisjunction& PartialBBDisjunction::operator=(const PartialBBDisjunction& source) {
  if (this != &source) {
    initialize(&source, NULL);
  }
  return *this;
} /* assignment operator */

PartialBBDisjunction* PartialBBDisjunction::clone() const {
  return new PartialBBDisjunction(*this);
} /* clone */

/** @details Set up the disjunction class as new (except the timer pointer, and do not reset params) */
void PartialBBDisjunction::setupAsNew() {
  VPCDisjunction::setupAsNew();
  this->root.initialize();
  this->data.num_nodes_on_tree = 0;
  this->data.num_partial_bb_nodes = 0;
  this->data.num_pruned_nodes = 0;
  this->data.min_node_depth = std::numeric_limits<int>::max();
  this->data.max_node_depth = 0;
  this->data.num_fixed_vars = 0;
  //    this->data.stats.resize(0);
  //    this->data.pruned_stats.resize(0);
  //    this->data.node_id.resize(0);
} /* setupAsNew */

/**
 * @details This will throw away all the information from the old disjunction, except it will not reset the timer
 * setupAsNew() is called, then a new Cbc object is created using a clone of \pi si
 * @returns a DisjExitReason for whether we were suceessful at generating the disjunction
 */
DisjExitReason PartialBBDisjunction::prepareDisjunction(const OsiSolverInterface* const si) {
  // Reset things in case we are reusing the class for some reason
  setupAsNew();
//  if (!timer) {
//    error_msg(errorstring, "Timer is not set.\n");
//    writeErrorToLog(errorstring, params.logfile);
//    verify(false, "VPC tried to exit with error code 1");
//  }
  if (this->params.get(intParam::DISJ_TERMS) < 2) {
    return DisjExitReason::NO_DISJUNCTION_EXIT;
  }

  if (this->params.get(DISJUNCTION_SOLVER) == "SYMPHONY") {
    generatePartialBBTreeSymphony(this, si);
    // we raise failures if when we encounter them above
    return DisjExitReason::SUCCESS_EXIT;
  }

  // Set up solver
  SolverInterface* BBSolver;
  try {
    BBSolver = dynamic_cast<SolverInterface*>(si->clone());
  } catch (std::exception& e) {
    error_msg(errorstring,
        "Unable to clone solver into desired SolverInterface.\n");
    writeErrorToLog(errorstring, params.logfile);
    verify(false, "VPC tried to exit with error code 1");
  }

  { // DEBUG DEBUG DEBUG
#ifdef BB_NEW_SOLVER
    try {
      if (BBSolver) { delete BBSolver; }
      BBSolver = new SolverInterface;
      setLPSolverParameters(BBSolver, params.get(VERBOSITY), params.get(TIMELIMIT));
      BBSolver->readMps(params.get(stringParam::FILENAME).c_str());
      if (BBSolver->getObjSense() < 1e-3) {
        BBSolver->setObjSense(1.0);
        const double* obj = BBSolver->getObjCoefficients();
        for (int col = 0; col < BBSolver->getNumCols(); col++) {
          BBSolver->setObjCoeff(col, -1. * obj[col]);
        }
        double objOffset = 0.;
        BBSolver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
        if (objOffset != 0.) {
          BBSolver->setDblParam(OsiDblParam::OsiObjOffset, -1. * objOffset);
        }
      }
      BBSolver->initialSolve();
    } catch (std::exception& e) {
      error_msg(errorstring,
          "Unable to create branch-and-bound solver into desired SolverInterface.\n");
      writeErrorToLog(errorstring, params.logfile);
      verify(false, "VPC tried to exit with error code 1");
    }
#endif
  } // DEBUG DEBUG DEBUG
  
  // Setup LP for use in partial tree generation
  setLPSolverParameters(BBSolver, params.get(VERBOSITY), std::numeric_limits<double>::max());
  BBSolver->setHintParam(OsiDoPresolveInInitial, false);
  BBSolver->setHintParam(OsiDoPresolveInResolve, false);
  BBSolver->setIntParam(OsiMaxNumIterationHotStart, std::numeric_limits<int>::max());
#ifdef USE_CLP
  try {
    dynamic_cast<OsiClpSolverInterface*>(BBSolver)->setSpecialOptions(16); // use standard strong branching rather than clp's
    // Do not switch from dual to primal, or something to this effect;
    // This allows infeasible branches to be fixed during strong branching
    dynamic_cast<OsiClpSolverInterface*>(BBSolver)->getModelPtr()->setMoreSpecialOptions(
        dynamic_cast<OsiClpSolverInterface*>(BBSolver)->getModelPtr()->moreSpecialOptions() + 256);
  } catch (std::exception& e) {
    std::cerr << "Unable to cast solver as OsiClpSolverInterface." << std::endl;
    verify(false, "VPC tried to exit with error code 1");
  }
#endif

  printf("## Generating partial branch-and-bound tree. ##\n");
  CbcModel* cbc_model = new CbcModel; // for obtaining the disjunction
  cbc_model->swapSolver(BBSolver);
  cbc_model->setModelOwnsSolver(true); // solver will be deleted with cbc object
  setIPSolverParameters(cbc_model, params.get(VERBOSITY));

  if (timer)
    timer->start_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);
  int num_strong = params.get(intParam::PARTIAL_BB_NUM_STRONG);
  if (num_strong == -1) {
    num_strong = si->getNumCols();
  }
  if (num_strong == -2) {
    num_strong = static_cast<int>(std::ceil(std::sqrt(si->getNumCols())));
  }
  const int num_before_trusted = std::numeric_limits<int>::max(); // 10;
  generatePartialBBTree(this, cbc_model, si, params.get(intParam::DISJ_TERMS),
      num_strong, num_before_trusted);
  printf("PartialBBDisjunction::prepareDisjunction: Finished generating partial branch-and-bound tree");
  if (timer) {
    timer->end_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);
    printf(" (%.3f seconds)", timer->get_time(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]));
  }
  printf("\n");

  VPCEventHandler* eventHandler = NULL;
  try {
    eventHandler = dynamic_cast<VPCEventHandler*>(cbc_model->getEventHandler());
  } catch (std::exception& e) {
  }
  if (!eventHandler) {
    error_msg(errstr, "Could not get event handler.\n");
    writeErrorToLog(errstr, params.logfile);
    verify(false, "VPC tried to exit with error code 1");
  }

#ifdef TRACE
  std::vector<NodeStatistics> stats = eventHandler->getStatsVector();
  printNodeStatistics(stats, false);
  if (eventHandler->getPrunedStatsVector().size() > 0) {
    printf("\n");
    printNodeStatistics(eventHandler->getPrunedStatsVector(), false);
  }
#endif
  
#ifdef TRACE
  const int TEMP_VAL = params.get(intParam::TEMP);
  if (use_temp_option(std::abs(TEMP_VAL), TempOptions::GEN_TIKZ_STRING)) {
    std::string dir, instname, ext;
    parseFilename(dir, instname, ext, params.get(FILENAME), params.logfile);
    const std::string filename_stub = dir + "/" + instname + "-round" + std::to_string(this->num_rounds);
    generateTikzTreeString(eventHandler, params, params.get(intParam::PARTIAL_BB_STRATEGY), si->getObjValue(), filename_stub);
    if (use_temp_option(TEMP_VAL, TempOptions::GEN_TIKZ_STRING_AND_RETURN)) {
      // Free
      if (BBSolver && !cbc_model->modelOwnsSolver()) { delete BBSolver; }
      if (cbc_model) { delete cbc_model; }
      return DisjExitReason::SUCCESS_EXIT;
    }
    if (use_temp_option(TEMP_VAL, TempOptions::GEN_TIKZ_STRING_AND_EXIT)) {
      // Free
      if (BBSolver && !cbc_model->modelOwnsSolver()) { delete BBSolver; }
      if (cbc_model) { delete cbc_model; }
      verify(false, "VPC tried to exit with error code 1"); // this should only be used during debug and will have memory leaks
    }
  }
#endif

  // If branch-and-bound finished (was not stopped by a user event), check why and exit
  if (cbc_model->status() == 0 || cbc_model->status() == 1
      || cbc_model->status() == 2 || this->num_terms <= 1) {
    if (cbc_model->isProvenOptimal()) {
      // Set strong branching lb and ub to both be this value
      this->best_obj = cbc_model->getObjValue();
      this->worst_obj = cbc_model->getObjValue();
      this->integer_obj = cbc_model->getObjValue();
      const double* sol = cbc_model->getColSolution();
      this->integer_sol.assign(sol, sol + si->getNumCols());
      if (BBSolver && !cbc_model->modelOwnsSolver()) { delete BBSolver; }
      if (cbc_model) { delete cbc_model; }
      return DisjExitReason::OPTIMAL_SOLUTION_FOUND_EXIT;
    }
    else if (this->num_terms <= 1) {
      warning_msg(warnstr,
          "Giving up on getting cuts from the partial branch-and-bound tree (too few terms). Model status is %d.\n",
          cbc_model->status());
      if (BBSolver && !cbc_model->modelOwnsSolver()) { delete BBSolver; }
      if (cbc_model) { delete cbc_model; }
      return DisjExitReason::TOO_FEW_TERMS_EXIT;
    }
    else {
      warning_msg(warnstr,
          "Giving up on getting cuts from the partial branch-and-bound tree (bad status). Model status is %d.\n",
          cbc_model->status());
      if (BBSolver && !cbc_model->modelOwnsSolver()) { delete BBSolver; }
      if (cbc_model) { delete cbc_model; }
      return DisjExitReason::UNKNOWN;
    }
  } // exit out early if cbc_model status is 0 or insufficiently many disjunctive terms

  // Make sure that the right number of terms has been saved - skip if we're saving
  // pruned terms because we don't know how many terms to expect from fixes between branching
  if (!params.get(intParam::PARTIAL_BB_KEEP_PRUNED_NODES)){
    if ((num_terms != eventHandler->getNumLeafNodes()) || (num_terms !=
        static_cast<int>(terms.size() + eventHandler->isIntegerSolutionFound()))){
      error_msg(errstr,
                "Number of terms does not match: num terms = %d, num leaf nodes = %d, num bases = %d, found_integer_sol = %d\n",
                num_terms, eventHandler->getNumLeafNodes(), static_cast<int>(terms.size()),
                eventHandler->isIntegerSolutionFound());
      writeErrorToLog(errstr, params.logfile);
      verify(false, "VPC tried to exit with error code 1");
    }
  } else {
    // exit early if we have a solution because we don't know how to handle it
    if (eventHandler->owner->integer_sol.size() > 0){
      printf("PartialBBDisjunction::prepareDisjunction: Handling for integer solution"
             "in partial tree not developed for full tree. Discarding tree.\n");
      return DisjExitReason::NO_DISJUNCTION_EXIT;
    }
  }

  if (BBSolver && !cbc_model->modelOwnsSolver()) { delete BBSolver; }
  if (cbc_model) { delete cbc_model; }
  return DisjExitReason::SUCCESS_EXIT;
} /* prepareDisjunction */

/**
 * @details Create a new disjunction that parameterizes the current with the given solver.
 * This updates in each term the optimal basis and objective value, as well as
 * the dual bounds and integer feasible solution in the disjunction.
 *
 * @param solver to update the disjunction with
 */
PartialBBDisjunction PartialBBDisjunction::parameterize(const OsiSolverInterface* const solver,
                                                        std::vector<std::unique_ptr<OsiSolverInterface>>* term_solvers) const {

  verify(this->common_changed_var.size() == 0 && this->common_changed_bound.size() == 0
         && this->common_changed_value.size() == 0 && this->common_ineqs.size() == 0,
         "Cannot parameterize a disjunction that has common terms or inequalities.");

  // If no vector is provided for the term solvers, use a temp local one
  std::vector<std::unique_ptr<OsiSolverInterface>> local_solvers;
  if (!term_solvers) {
    term_solvers = &local_solvers;
  }

  // get a copy of the solver
  SolverInterface* si = dynamic_cast<SolverInterface*>(solver->clone());
  ensureMinimizationObjective(si); // minimize to keep meaning of disjunctive terms consistent
  si->resolve();
  verify(checkSolverOptimality(si, true), "solver must be feasible");

  // create an empty disjunction
  PartialBBDisjunction disj = PartialBBDisjunction(this->params);

  // reset number of infeasible terms
  disj.num_infeasible_terms = 0;

  // update the root LP relaxation
  disj.root_obj = si->getObjValue();

  // update the integer feasible solution if possible
  if (this->integer_sol.size() > 0 && isFeasible(*si, this->integer_sol)){
    disj.integer_sol = this->integer_sol;
    double obj = 0;
    for (int i = 0; i < si->getNumCols(); i++) {
      obj += si->getObjCoefficients()[i] * this->integer_sol[i];
    }
    disj.integer_obj = obj;
  }

  // parameterize each term
  for (int term_idx = 0; term_idx < this->num_terms; term_idx++){

    // get the solver
    OsiSolverInterface* termSolver;
    this->getSolverForTerm(termSolver, term_idx, si, false, .001, NULL, true);

    // get the term
    DisjunctiveTerm term = this->terms[term_idx];

    // update the necessary parts of the term
    term.is_feasible = checkSolverOptimality(termSolver, true);
    term.obj = term.is_feasible ? termSolver->getObjValue() : std::numeric_limits<double>::max();

    // update the necessary disjunction metadata
    disj.updateObjValue(term.obj);
    disj.terms.push_back(term);
    disj.num_infeasible_terms += !term.is_feasible;
    disj.num_terms++;

    // cache the solver
    term_solvers->emplace_back(std::unique_ptr<OsiSolverInterface>(termSolver->clone()));
  }

  // sanity check - egregious errors should be avoided by not counting objectives from infeasible terms
  if (!lessThanVal(si->getObjValue(), disj.best_obj, -1e-7)) {
    std::cout << "Warning: disjunctive dual bound should not be less than LP"
                 "relaxation objective value. Sometimes this happens due to numerical issues." << std::endl;
  }

  // return the parameterized disjunction
  return disj;
}

//<--***************** PROTECTED **********************-->

/**
 * @details If \p source != NULL, clear #root then point #root to \p source->root, and similarly for #data
 * Otherwise, call RootTerm::initialize() and PartialBBDisjunctionData::initialize()
 */
void PartialBBDisjunction::initialize(const PartialBBDisjunction* const source,
    const VPCParametersNamespace::VPCParameters* const params) {
  VPCDisjunction::initialize(source, params);
  if (source != NULL) {
    this->root.clear();
    this->root = source->root;
    this->data = source->data;
  } else {
    this->root.initialize();
    this->data.initialize();
  }
} /* initialize */

/// @details \p choose_strategy can be positive or negative
void getCbcStrategies(
    /// [in] parameters to read settings from
    const VPCParametersNamespace::VPCParameters& params,
    /// [out] ones digit of intParam::PARTIAL_BB_STRATEGY; given a tree, which node to pick next
    int& compare_strategy,
    /// [out] tens digit of intParam::PARTIAL_BB_STRATEGY; given a branching variable, which direction to choose
    int& branch_strategy,
    /// [out] hundreds digit of intParam::PARTIAL_BB_STRATEGY; which branching variable to pick
    int& choose_strategy) {
  // What is the partial strategy?
  const int strategy = std::abs(params.get(intParam::PARTIAL_BB_STRATEGY));
  const int sign = (params.get(intParam::PARTIAL_BB_STRATEGY) < 0) ? -1 : 1;
  compare_strategy = strategy % 10; // ones digit
  branch_strategy = (strategy % 100 - compare_strategy) / 10; // tens digit
  choose_strategy = sign * (strategy % 1000 - compare_strategy - 10 * branch_strategy) / 100; // hundreds digit
} /* getCbcStrategies */

#ifdef USE_CBC
/**
 * @details Set maximum number of hot start iterations to 100,
 * set the CbcBranchDecision, CbcCompareBase, CbcEventHandler, and CbcModel::setChooseMethod().
 * Also disable presolve and cuts, set maximum time equal to \p max_time,
 * and set maximum number of strong branching candidates and number of nodes explored before switching to pseudocosts
 */
void setCbcParametersForPartialBB(
    /// [in,out] CbcModel we are deciding settings for
    CbcModel* const cbc_model,
    /// [in] how much information to print, but (in Cbc 2.10 and on) we need to set CbcModel::setSecsPrintFrequency() to -1 to reach the treeStatus event correctly
    const int verbosity,
    /// [in] given a tree, which node to pick next
    const int compare_strategy,
    /// [in] given a branching variable, which direction to choose
    const int branch_strategy,
    /// [in] which branching variable to pick
    const int choose_strategy,
    /// [in] maximum number of strong branching candidates
    const int numStrong,
    /// [in] number of explored nodes before switching to pseudocosts; 0 disables dynamic strong branching
    const int numBeforeTrusted,
    /// [in] time limit for Cbc to generate a tree
    const double max_time,
    /// [in,out] custom event handler (\sa VPCEventHandler); will not be changed here, but not deleted either
    CbcEventHandler* eventHandler) {
  setIPSolverParameters(cbc_model, verbosity);
#ifdef CBC_VERSION_210PLUS
  cbc_model->setSecsPrintFrequency(-1); // to reach tree status correctly
#endif
  cbc_model->solver()->setIntParam(OsiMaxNumIterationHotStart, 100);

  // Branching decision (tens digit)
  // Given a branching variable, which direction to choose?
  // (The choice of branching variable is through OsiChooseVariable)
  // 0: default, 1: dynamic, 2: strong, 3: none
  CbcBranchDecision* branch;
  if (branch_strategy == 1) {
    branch = new CbcBranchDynamicDecision();
  } else if (branch_strategy == 2) {
    branch = new CbcBranchStrongDecision();
  } else if (branch_strategy == 3) {
    branch = NULL;
  } else {
    branch = new CbcBranchDefaultDecision();
  }

  // Set comparison for nodes (ones digit)
  // Given a tree, which node to pick next?
  // 0: default: 1: bfs, 2: depth, 3: estimate, 4: objective, 5: objective_reverse
  CbcCompareBase* compare;
  if (compare_strategy == 1) {
    compare = new CbcCompareBFS();
  } else if (compare_strategy == 2) {
    compare = new CbcCompareDepth();
  } else if (compare_strategy == 3) {
    compare = new CbcCompareEstimate();
  } else if (compare_strategy == 4) {
    compare = new CbcCompareObjective();
  } else {
    compare = new CbcCompareDefault();
  }

  cbc_model->setTypePresolve(0);
  cbc_model->setMaximumSeconds(max_time);
  cbc_model->setMaximumCutPassesAtRoot(0);
  cbc_model->setMaximumCutPasses(0);
  cbc_model->setWhenCuts(0);
  if (numStrong >= 0) {
    // Maximum number of strong branching candidates to consider each time
    cbc_model->setNumberStrong(numStrong);
  }
  if (numBeforeTrusted >= 0) {
    // # before switching to pseudocosts, I think; 0 disables dynamic strong branching, doesn't work well
    cbc_model->setNumberBeforeTrust(numBeforeTrusted);
  }

  if (branch) {
    OsiChooseStrongCustom choose;
    if (numStrong >= 0) {
      choose.setNumberStrong(numStrong);
    }
    if (numBeforeTrusted >= 0) {
      choose.setNumberBeforeTrusted(numBeforeTrusted);
    }
    choose.setMethod(choose_strategy);
    branch->setChooseMethod(choose);
    // From CbcModel::convertToDynamic, we see that the branching decision may be ignored without a choose method
    if ((branch->whichMethod()&1) == 0 && !branch->chooseMethod()) {
      OsiChooseStrong choose;
      choose.setNumberStrong(5);
      choose.setNumberBeforeTrusted(0);
      branch->setChooseMethod(choose);
    }
    cbc_model->setBranchingMethod(*branch);
  } // check that branch is not NULL

  if (eventHandler) {
    cbc_model->passInEventHandler(eventHandler);
  }

  cbc_model->setNodeComparison(compare);

  if (branch) {
    delete branch;
  }
  if (compare) {
    delete compare;
  }
} /* setCbcParametersForPartialBB */

//<!--***********************************************************-->
void generatePartialBBTree(PartialBBDisjunction* const owner, CbcModel* cbc_model,
    const OsiSolverInterface* const solver, const int max_leaf_nodes,
    const int num_strong, const int num_before_trusted) {
  //  const double partial_timelimit = 100 * max_leaf_nodes * solver->getNumCols()
  //      * GlobalVariables::timeStats.get_time(GlobalConstants::INIT_SOLVE_TIME);
  const double partial_timelimit = owner->params.get(PARTIAL_BB_TIMELIMIT); // will be checked manually by the eventHandler

  // Set up options
  VPCEventHandler* eventHandler = new VPCEventHandler(owner, max_leaf_nodes, partial_timelimit);
  eventHandler->setOriginalSolver(solver);

  // This sets branching decision, event handling, etc.
  int compare_strategy, branch_strategy, choose_strategy;
  getCbcStrategies(owner->params, compare_strategy, branch_strategy, choose_strategy);
  setCbcParametersForPartialBB(cbc_model, owner->params.get(intParam::VERBOSITY),
      compare_strategy, branch_strategy, choose_strategy, num_strong,
      num_before_trusted, std::numeric_limits<double>::max(), eventHandler);

#ifdef TRACE
  cbc_model->branchAndBound(3);
#else
  cbc_model->branchAndBound(0);
#endif

  // Free
  // When eventHandler is passed, Cbc currently clones it and does not delete the original
  if (eventHandler) { // in case this behavior gets changed in the future
    delete eventHandler;
  }
} /* generatePartialBBTree */
#endif // USE_CBC

void generatePartialBBTreeSymphony(
    PartialBBDisjunction* const owner,
    const OsiSolverInterface* const solver) {

  BBInfo info;
  std::shared_ptr <OsiSymSolverInterface> parametric_model = std::shared_ptr<OsiSymSolverInterface>();
  int node_limit = owner->params.get(intParam::BB_NODE_LIMIT);

  // generate a partial BB tree with DISJ_TERMS leaf nodes
  owner->params.set(intParam::BB_NODE_LIMIT, owner->params.get(intParam::DISJ_TERMS));
  doBranchAndBoundWithSymphony(
      owner->params, owner->params.get(VPCParametersNamespace::BB_STRATEGY),
      solver, info, nullptr, parametric_model);
  owner->params.set(intParam::BB_NODE_LIMIT, node_limit);

  // get its root node
  sym_environment * env = parametric_model->getSymphonyEnvironment();
  warm_start_desc* ws = sym_get_warm_start(env, true);
  bc_node* root = ws->rootnode;

  // set up the disjunction object
  owner->root_obj = root->lower_bound;

  // bfs on nodes. track branching decisions to reach each node ordered from root to leaf
  std::deque<bc_node*> nodes{root};
  std::deque<std::vector<int>> branch_paths{{}};
  while (!nodes.empty()){

    // get current node
    bc_node* current = nodes.front();  // current node
    nodes.pop_front();
    // and path to it (represented by the child index each node in the path was for its parent)
    std::vector<int> branch_path = branch_paths.front();
    branch_paths.pop_front();

    if (current->bobj.child_num > 0) {
      // internal node - add children and their indices to bfs queues
      // indices used later to traverse from root to leaf
      for (int i = 0; i < current->bobj.child_num; i++){
        nodes.push_back(current->children[i]);
        branch_paths.push_back(branch_path);
        branch_paths.back().push_back(i);
      }
    } else {
      // leaf node - create disjunctive term
      DisjunctiveTerm term;
      term.type = "leaf";

      // traverse path from root to current node to get branching decisions
      verify(branch_path.size() == current->bc_level,
             "branch path size should equal node depth");
      bc_node* path_node = root;
      for (int level = 0; level < branch_path.size(); level++){
        int child_idx = branch_path[level];

        // L is <=, thus upper bound -> term.changed_bound is 1
        // G is >=, thus lower bound -> term.changed_bound is 0
        // source: Disjunction.hpp::DisjunctiveTerm
        char sense = path_node->bobj.sense[child_idx];
        verify(sense == 'L' || sense == 'G', "VPC only supports <= and >= branching");
        term.changed_var.push_back(path_node->bobj.position);
        term.changed_bound.push_back(sense == 'L' ? 1 : 0);
        term.changed_value.push_back(path_node->bobj.rhs[child_idx]);

        if (level < branch_path.size() - 1) {
          // move to next node in path
          path_node = path_node->children[child_idx];
        } else {
          verify(path_node->children[child_idx] == current,
                 "current node should be end of path");
        }
      }

      // add term to disjunction
      owner->terms.push_back(term);
      owner->num_terms++;
    }
  }

  // call parameterize with existing solver to set feasibility, objective values, and warm starts
  std::vector<std::unique_ptr<OsiSolverInterface>> term_solvers;
  *owner = owner->parameterize(solver, &term_solvers);
  for (int i = 0; i < owner->num_terms; i++) {
    // save the basis for each term - currently in added row format b/c pdc expects this format
    // not bothering with tightened bound format since we don't have it handy here and VPC can figure it out
    owner->terms[i].basis_extended = dynamic_cast<CoinWarmStartBasis*>(term_solvers[i]->getWarmStart());
  }

  // verify that the disjunction represents the leaves of a full binary tree
  owner->isFullBinaryTree();

} /* generatePartialBBTree */

/**
 * @details raise an assert if the disjunction does not represent the leaves of a full binary tree.
 * This function assumes that the branching decisions in each disjunctive term
 * are sorted in the order they occurred.
 *
 * @return void
 */
void PartialBBDisjunction::isFullBinaryTree(){

  std::vector<DisjunctiveTerm> terms_copy = terms;

  // check that term does not contain term2
  // we can get away with just checking for ancestral relationship because it
  // isn't possible to end up with the same branching decisions just in different order
  for (const DisjunctiveTerm& term : terms_copy){
    for (const DisjunctiveTerm& term2 : terms_copy){

      // term can only contain term2 if term2 is at least as deep in the tree
      if (&term != &term2 && term.changed_var.size() <= term2.changed_var.size()){

        // make sure that term is not an ancestor of term2 by checking they branch
        // on different variables or at least in different directions
        verify(!std::equal(term2.changed_var.begin(), term2.changed_var.begin() +
                                                      term.changed_var.size(), term.changed_var.begin()) ||
               !std::equal(term2.changed_bound.begin(), term2.changed_bound.begin() +
                                                        term.changed_bound.size(), term.changed_bound.begin()),
               "the LP relaxation of term contains the LP relaxation of term2");
      }
    }
  }

  // get the maximum depth
  int max_depth = 0;
  for (DisjunctiveTerm term : terms_copy){
    max_depth = term.changed_var.size() > max_depth ? term.changed_var.size() : max_depth;
  }

  // check each leaf has a sibling
  for (int depth = max_depth; depth > 0; depth--){

    // get the terms at this depth
    std::vector<DisjunctiveTerm> depth_terms;
    for (DisjunctiveTerm term : terms_copy){
      if (term.changed_var.size() == depth){
        depth_terms.push_back(term);
      }
    }

    // keep a running list of paired terms
    std::set<const DisjunctiveTerm*> paired_terms;

    // find a sibling for each term
    for (const DisjunctiveTerm& term : depth_terms){

      // Check if the term was found to be another's sibling earlier
      if (paired_terms.find(&term) != paired_terms.end()){
        continue;
      }

      // check term against all other terms at this depth
      for (const DisjunctiveTerm& term2 : depth_terms){
        std::vector<int> differing_idx;

        // we can only be siblings if we share the same variables that were
        // branched on and we're not the same term
        if (term.changed_var == term2.changed_var && &term != &term2){

          // record branching directions we differ on
          for (int i = 0; i < depth; i++) {
            if (term.changed_bound[i] != term2.changed_bound[i]) {
              differing_idx.push_back(i);
            }
          }

          // if we differ by only the last branching decision, we're siblings
          if (differing_idx.size() == 1 && differing_idx[0] == depth - 1){
            // check that we have the reciprocal branching decision
            double expected_val = term.changed_value[differing_idx[0]] +
                                  (term.changed_bound[differing_idx[0]] == 0 ? -1 : 1);
            verify(term2.changed_value[differing_idx[0]] == expected_val,
                   "term2 branch value doesnt meet expectation");

            // record siblings
            paired_terms.insert(&term);
            paired_terms.insert(&term2);

            // create a parent node to leave in the next level above
            DisjunctiveTerm parent_term = term;
            parent_term.changed_var.erase(parent_term.changed_var.begin() + differing_idx[0]);
            parent_term.changed_bound.erase(parent_term.changed_bound.begin() + differing_idx[0]);
            parent_term.changed_value.erase(parent_term.changed_value.begin() + differing_idx[0]);
            parent_term.type = "parent";
            terms_copy.push_back(parent_term);
            break;
          }
        }
      }

      // if we didn't find a sibling, the disjunction is not complete
      if (paired_terms.find(&term) == paired_terms.end()){
        verify(false, "Disjunction does not represent a full binary tree.");
      }

    } // find a sibling for each term
  }
} /* isFullBinaryTree */
