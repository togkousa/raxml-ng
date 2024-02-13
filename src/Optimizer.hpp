#ifndef RAXML_OPTIMIZER_H_
#define RAXML_OPTIMIZER_H_

#include "TreeInfo.hpp"
#include "Checkpoint.hpp"

#include <stdio.h>


class Optimizer
{
public:
  Optimizer (const Options& opts);
  virtual
  ~Optimizer ();

  double optimize_model(TreeInfo& treeinfo, double lh_epsilon);
  double optimize_model(TreeInfo& treeinfo) { return optimize_model(treeinfo, _lh_epsilon); };
  
  // optimization routines
  double optimize_topology(TreeInfo& treeinfo, CheckpointManager& cm);
  double optimize_topology_noise(TreeInfo& treeinfo, CheckpointManager& cm);
  double optimize_topology_adaptive(TreeInfo& treeinfo, CheckpointManager& cm, double difficulty);
  
  double evaluate(TreeInfo& treeinfo, CheckpointManager& cm);
  void nni(TreeInfo& treeinfo, nni_round_params& nni_params, double& loglh);

  bool noise_quantification_mode() { return ( _msa_error_rate || _sampling_noise); }

  std::string _msa_error_file;

private:
  double _lh_epsilon;
  double _lh_epsilon_brlen_triplet;
  int _spr_radius;
  double _spr_cutoff;

  // nni params
  double _nni_epsilon;
  double _nni_tolerance;

  // msa error-rate
  double _msa_error_rate;
  bool _msa_error_randomized;
  bool _msa_error_blo;
  bool _sampling_noise;
  int _sampling_noise_mode;
  unsigned int _seed;
  bool _count_spr_moves;

  // functions for adaptive mode
  double convergence_rate(CheckpointManager& cm, double test_loglh);
  bool first_search_done(CheckpointManager& cm);
  bool converged(CheckpointManager& cm, double test_loglh, double epsilon);
  int adaptive_slow_spr_radius(double difficulty);
  
  double sampling_noise_epsilon(corax_treeinfo_t* treeinfo, 
                              double** persite_lnl, 
                              double loglh, 
                              int sampling_noise_mode);

  double **initialize_persite_vector(TreeInfo& treeinfo);
  void free_persite_vector(TreeInfo& treeinfo, double ** persite_lnl);

  double kh_like(TreeInfo& treeinfo, double** persite_lnl, double** persite_lnl_new);

  // void delete_save_to_file(std::string filename, int counter, TreeInfo& treeinfo, double ** persite_lnl);

};

#endif /* RAXML_OPTIMIZER_H_ */
