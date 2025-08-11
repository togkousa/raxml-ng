#ifndef RAXML_OPTIMIZER_H_
#define RAXML_OPTIMIZER_H_

#include "TreeInfo.hpp"
#include "Checkpoint.hpp"

#include <stdio.h>
#include "adaptive/StoppingCriterion.hpp"

class Optimizer
{
public:
  Optimizer (const Options& opts);
  virtual
  ~Optimizer ();

  double optimize_model(TreeInfo& treeinfo, double lh_epsilon, bool testing_sites = false);
  double optimize_model(TreeInfo& treeinfo, bool testing_sites = false) { return optimize_model(treeinfo, _lh_epsilon, testing_sites); };
  
  double evaluate_testing_sites(TreeInfo& training_treeinfo, 
                                TreeInfo* testing_treeinfo,
                                double br_len_epsilon,
                                double mod_opt_epsilon,
                                bool substitute = true);

  // optimization routines
  double optimize_topology(TreeInfo& treeinfo, TreeInfo* treeinfo_testing, CheckpointManager& cm, PartitionedMSA& parted_msa);
  double optimize_topology_adaptive(TreeInfo& treeinfo, TreeInfo* treeinfo_testing, CheckpointManager& cm, PartitionedMSA& parted_msa);
  double optimize_topology_modified(TreeInfo& treeinfo, TreeInfo* treeinfo_testing, CheckpointManager& cm, PartitionedMSA& parted_msa);
  
  double evaluate(TreeInfo& treeinfo, CheckpointManager& cm, PartitionedMSA& parted_msa);
  void nni(TreeInfo& treeinfo, nni_round_params& nni_params, double& loglh);

  void set_stopping_criterion(StoppingCriterion* _criterion) { criterion = _criterion; } // have to fix this

private:
  double _lh_epsilon;
  double _lh_epsilon_brlen_triplet;
  int _spr_radius;
  double _spr_cutoff;

  // nni params
  double _nni_epsilon;
  double _nni_tolerance;

  // stoping criteria
  int _stopping_criterion;
  bool _modified_version;
  StoppingCriterion *criterion;

  // cross validation
  bool _use_holdout_es;
  int _convergence_iters;
  corax_unode_t *_best_holdout_es_solution;
  double _best_holdout_es_loglh;

  // functions for adaptive mode
  int fast_spr_radius_adaptive(double difficulty);
  bool call_modified_version() {return _modified_version; }

  // houldout es fucntion
  double revert_to_best_solution(TreeInfo& training_treeinfo, 
                                TreeInfo* testing_treeinfo,
                                double br_len_epsilon,
                                double mod_opt_epsilon,
                                double &training_loglh,
                                double &test_loglh);
  
  void set_best_holdout_es_solution(const corax_unode_t* root, double testing_loglh);
};

#endif /* RAXML_OPTIMIZER_H_ */
