#include "Optimizer.hpp"
#include "topology/RFDistCalculator.hpp"
#include "adaptive/MSAErrorHandler.hpp"
#include "adaptive/KHtester.hpp"
#include <memory>

#define _USE_MATH_DEFINES 
#include <cmath>


/* 
#include <thread>         // std::thread
#include <mutex>          // std::mutex
 */
using namespace std;

Optimizer::Optimizer (const Options &opts) :
    _msa_error_file(opts.msa_error_file),
    _lh_epsilon(opts.lh_epsilon), _lh_epsilon_brlen_triplet(opts.lh_epsilon_brlen_triplet), 
    _spr_radius(opts.spr_radius), _spr_cutoff(opts.spr_cutoff), 
    _nni_epsilon(opts.nni_epsilon), _nni_tolerance(opts.nni_tolerance), 
    _msa_error_rate(opts.msa_error_rate),
    _msa_error_randomized(opts.msa_error_randomized), _msa_error_blo(opts.msa_error_blo),
    _sampling_noise(opts.sampling_noise), _noise_rell(opts.noise_rell) , _seed(opts.random_seed)
{
} 

Optimizer::~Optimizer ()
{
  // TODO Auto-generated destructor stub
}

double Optimizer::optimize_model(TreeInfo& treeinfo, double lh_epsilon)
{
  double new_loglh = treeinfo.loglh();

//  if (!params_to_optimize)
//    return new_loglh;

  int iter_num = 0;
  double cur_loglh;
  do
  {
    cur_loglh = new_loglh;

    treeinfo.optimize_params_all(lh_epsilon);

    new_loglh = treeinfo.loglh();

//      printf("old: %f, new: %f\n", cur_loglh, new_loglh);

    iter_num++;
    LOG_DEBUG << "Iteration " << iter_num <<  ": logLH = " << new_loglh << endl;
  }
  while (new_loglh - cur_loglh > lh_epsilon);

  return new_loglh;
}

void Optimizer::nni(TreeInfo& treeinfo, nni_round_params& nni_params, double& loglh){
  // nni round
  LOG_PROGRESS(loglh) << "NNI round tolerance = " <<  nni_params.tolerance << ", epsilon = " << nni_params.lh_epsilon << endl;
  loglh = treeinfo.nni_round(nni_params);
  
}

double Optimizer::optimize_topology(TreeInfo& treeinfo, CheckpointManager& cm)
{
  const double fast_modopt_eps = 10.;
  const double interim_modopt_eps = 3.;
  const double final_modopt_eps = 0.1;

  SearchState local_search_state = cm.search_state();
  auto& search_state = ParallelContext::group_master_thread() ? cm.search_state() : local_search_state;
  ParallelContext::barrier();

  /* set references such that we can work directly with checkpoint values */
  double &loglh = search_state.loglh;
  int& iter = search_state.iteration;
  spr_round_params& spr_params = search_state.spr_params;
  int& best_fast_radius = search_state.fast_spr_radius;
  spr_params.lh_epsilon_brlen_full = _lh_epsilon;
  spr_params.lh_epsilon_brlen_triplet = _lh_epsilon_brlen_triplet;

  // msa-error-rate configuration
  std::shared_ptr<MSAErrorHandler> msa_error_handler;
  double **persite_lnl;

  if(_msa_error_rate){

    assert(treeinfo.kh_test() == false);
    assert(_sampling_noise == false);

    msa_error_handler = make_shared<MSAErrorHandler>(&treeinfo.pll_treeinfo(),
                                                      _msa_error_rate,
                                                      _seed,
                                                      _msa_error_file);
  }

  CheckpointStep resume_step = search_state.step;

  /* Compute initial LH of the starting tree */
  loglh = treeinfo.loglh();

  auto do_step = [&search_state,resume_step](CheckpointStep step) -> bool
      {
        if (step >= resume_step)
        {
          search_state.step = step;
          return true;
        }
        else
          return false;;
      };

  if (do_step(CheckpointStep::brlenOpt))
  {
    cm.update_and_write(treeinfo);
    LOG_PROGRESS(loglh) << "Initial branch length optimization" << endl;
    loglh = treeinfo.optimize_branches(fast_modopt_eps, 1);
  }

  /* Initial fast model optimization */
  if (do_step(CheckpointStep::modOpt1))
  {
    cm.update_and_write(treeinfo);
    LOG_PROGRESS(loglh) << "Model parameter optimization (eps = " << fast_modopt_eps << ")" << endl;
    loglh = optimize_model(treeinfo, fast_modopt_eps);

    /* start spr rounds from the beginning */
    iter = 0;
  }

  if(_sampling_noise){

    LOG_PROGRESS(loglh) << "Noise Sampling Using the " << (_noise_rell ? "RELL" : "NO-RELL") << " apporach." << endl;
    assert(_msa_error_rate == 0);
    unsigned int partition_count = treeinfo.pll_treeinfo().partition_count;
    
    persite_lnl = new double*[partition_count];
    for (unsigned int part = 0; part<partition_count; part++)
      persite_lnl[part] = new double[treeinfo.pll_treeinfo().partitions[part]->sites];
    
    corax_treeinfo_t* tmp_treeinfo = &treeinfo.pll_treeinfo_unconst();
    loglh = corax_treeinfo_compute_loglh_persite(tmp_treeinfo, 1, 0, persite_lnl);
    _lh_epsilon = sampling_noise_epsilon(tmp_treeinfo, persite_lnl, loglh, _noise_rell);
    
    LOG_PROGRESS(loglh) << (_noise_rell ? "RELL" : "NO-RELL") << " apporach. Epsilon = " << _lh_epsilon << endl;
  }

  // Do it once here !!!
  unsigned int dist_size = 1000;
  if(_msa_error_rate)
    msa_error_handler->msa_error_dist(treeinfo, dist_size, loglh, true, fast_modopt_eps, _msa_error_blo);
  
  // do SPRs
  const int radius_limit = min(22, (int) treeinfo.pll_treeinfo().tip_count - 3 );
  const int radius_step = 5;

  if (_spr_radius > 0)
    best_fast_radius = _spr_radius;
  else
  {
    /* auto detect best radius for fast SPRs */

    if (do_step(CheckpointStep::radiusDetect))
    {
      if (iter == 0)
      {
        spr_params.thorough = 0;
        spr_params.radius_min = 1;
        best_fast_radius = spr_params.radius_max = 5;
        spr_params.ntopol_keep = 0;
        spr_params.subtree_cutoff = 0.;
      }

      double best_loglh = loglh;
      double epsilon;

      while (spr_params.radius_min < radius_limit)
      {
        epsilon = _msa_error_rate ? msa_error_handler->draw_proportionately_from_distribution(_msa_error_randomized) : _lh_epsilon ;

        cm.update_and_write(treeinfo);

        ++iter;
        LOG_PROGRESS(best_loglh) << "AUTODETECT spr round " << iter << " (radius: " <<
            spr_params.radius_max << "), epsilon = " << epsilon << endl;
        loglh = treeinfo.spr_round(spr_params);

        if (loglh - best_loglh > epsilon)
        {
          /* LH improved, try to increase the radius */
          best_fast_radius = spr_params.radius_max;
          spr_params.radius_min += radius_step;
          spr_params.radius_max += radius_step;
          best_loglh = loglh;
        }
        else
          break;
      }
    }
  }

  LOG_PROGRESS(loglh) << "SPR radius for FAST iterations: " << best_fast_radius << " (" <<
                 (_spr_radius > 0 ? "user-specified" : "autodetect") << ")" << endl;

  if (do_step(CheckpointStep::modOpt2))
  {
    cm.update_and_write(treeinfo);

    /* optimize model parameters a bit more thoroughly */
    LOG_PROGRESS(loglh) << "Model parameter optimization (eps = " <<
                                                            interim_modopt_eps << ")" << endl;
    loglh = optimize_model(treeinfo, interim_modopt_eps);

    /* reset iteration counter for fast SPRs */
    iter = 0;

    /* initialize search params */
    spr_params.thorough = 0;
    spr_params.radius_min = 1;
    spr_params.radius_max = best_fast_radius;
    spr_params.ntopol_keep = 20;
    spr_params.subtree_cutoff = _spr_cutoff;
    spr_params.reset_cutoff_info(loglh);
  }

  double old_loglh;

  if (do_step(CheckpointStep::fastSPR))
  {
    double epsilon;
    do
    { 
      epsilon = _msa_error_rate ? msa_error_handler->draw_proportionately_from_distribution(_msa_error_randomized) : _lh_epsilon ;
      
      cm.update_and_write(treeinfo);
      ++iter;

      old_loglh = loglh;
      LOG_PROGRESS(old_loglh) << (spr_params.thorough ? "SLOW" : "FAST") <<
          " spr round " << iter << " (radius: " << spr_params.radius_max << "), epsilon = " << epsilon << endl;
      loglh = treeinfo.spr_round(spr_params);

      /* optimize ALL branches */
      loglh = treeinfo.optimize_branches(_lh_epsilon, 1);

    }
    while (loglh - old_loglh > epsilon);
  }

  if (do_step(CheckpointStep::modOpt3))
  {
    cm.update_and_write(treeinfo);
    LOG_PROGRESS(loglh) << "Model parameter optimization (eps = " << 1.0 << ")" << endl;
    loglh = optimize_model(treeinfo, 1.0);

    /* init slow SPRs */
    spr_params.thorough = 1;
    spr_params.radius_min = 1;
    spr_params.radius_max = 10;
    iter = 0;
  }
  
  bool impr = true;
  if (do_step(CheckpointStep::slowSPR))
  {
    double epsilon;

    do
    {
      epsilon = _msa_error_rate ? msa_error_handler->draw_proportionately_from_distribution(_msa_error_randomized) : _lh_epsilon ;

      cm.update_and_write(treeinfo);
      ++iter;
      old_loglh = loglh;

      LOG_PROGRESS(old_loglh) << (spr_params.thorough ? "SLOW" : "FAST") <<
          " spr round " << iter << " (radius: " << spr_params.radius_max << "), epsilon = " << epsilon << endl;
      loglh = treeinfo.spr_round(spr_params);
      loglh = treeinfo.optimize_branches(_lh_epsilon, 1);

      impr = (loglh - old_loglh > epsilon);
      /* if (impr)
      {
        spr_params.radius_min = 1;
        spr_params.radius_max = radius_step;
      }
      else
      {
        spr_params.radius_min = spr_params.radius_max + 1;
        spr_params.radius_max += radius_step;
      } */
    }
    while (impr);
  }

  /* Final thorough model optimization */
  if (do_step(CheckpointStep::modOpt4))
  {
    cm.update_and_write(treeinfo);
    LOG_PROGRESS(loglh) << "Model parameter optimization (eps = " << final_modopt_eps << ")" << endl;
    loglh = optimize_model(treeinfo, final_modopt_eps);
  }

  if (do_step(CheckpointStep::finish))
    cm.update_and_write(treeinfo);
  
  if(_msa_error_rate && _msa_error_file.length() > 0)
    msa_error_handler->msa_error_dist(treeinfo, dist_size, loglh, false, fast_modopt_eps, _msa_error_blo);

  if(_sampling_noise){
    
    unsigned int partition_count = treeinfo.pll_treeinfo().partition_count;
    for (unsigned int part = 0; part<partition_count; part++)
      delete[] persite_lnl[part];

    delete[] persite_lnl;

  }

  return loglh;
}

double Optimizer::optimize_topology_adaptive(TreeInfo& treeinfo, CheckpointManager& cm, double difficulty)
{
  // TODO: connect the command line arguments for nni-epsilon and nni-tolerance with nni_params.lh_epsilon and 
  // nni_params.tolerance
  const double fast_modopt_eps = 10.;
  const double interim_modopt_eps = 3.;
  const double final_modopt_eps = 0.1;

  // to store all the intermediate trees
  TreeList intermediate_trees;

  SearchState local_search_state = cm.search_state();
  auto& search_state = ParallelContext::group_master_thread() ? cm.search_state() : local_search_state;
  ParallelContext::barrier();

  /* set references such that we can work directly with checkpoint values */
  double &loglh = search_state.loglh;
  
  // spr round - basics
  spr_round_params& spr_params = search_state.spr_params;
  int slow_spr_radius = adaptive_slow_spr_radius(difficulty); // slow spr radius is determined by difficulty
  slow_spr_radius = min(slow_spr_radius, (int) treeinfo.pll_treeinfo().tip_count - 3 );
  spr_params.lh_epsilon_brlen_full = _lh_epsilon;
  spr_params.lh_epsilon_brlen_triplet = _lh_epsilon_brlen_triplet;
  int iter = 0;

  // nni round - basics
  nni_round_params& nni_params = search_state.nni_params;
  nni_params.tolerance = _nni_tolerance;
  nni_params.lh_epsilon = _nni_epsilon;

  // rf distance calculator
  std::unique_ptr<RFDistCalculator> rfDist(new RFDistCalculator());

  /* Compute initial LH of the starting tree */
  loglh = treeinfo.loglh();

  /* Initial branch length optimization */
  cm.update_and_write(treeinfo);
  LOG_PROGRESS(loglh) << "Initial branch length optimization" << endl;
  loglh = treeinfo.optimize_branches(fast_modopt_eps, 1);
  
  /* Initial fast model optimization */
  cm.update_and_write(treeinfo);
  LOG_PROGRESS(loglh) << "Model parameter optimization (eps = " << fast_modopt_eps << ")" << endl;
  loglh = optimize_model(treeinfo, fast_modopt_eps);
  
  /* push back the initial tree */
  intermediate_trees.emplace_back(treeinfo.tree());
  rfDist->set_tree_list(intermediate_trees);

  // If the dataset is "easy" or "difficult", start with an NNI round
  if(difficulty < 0.3 || difficulty > 0.7){
    
    cm.update_and_write(treeinfo);
    nni(treeinfo, nni_params, loglh);

    // + model parameter optimization
    LOG_PROGRESS(loglh) << "Model parameter optimization (eps = " << interim_modopt_eps << ")" << endl;
    loglh = optimize_model(treeinfo, fast_modopt_eps);

    intermediate_trees.emplace_back(treeinfo.tree());

  }

  // do SPRs
  const int radius_limit = min(22, (int) treeinfo.pll_treeinfo().tip_count - 3 );
  const int radius_step = 5; // HAVE TO CHANGE THAT

  // setting up fast SPR parameters
  spr_params.thorough = 0;
  spr_params.radius_min = 1;
  spr_params.radius_max = 5;
  spr_params.ntopol_keep = 0;
  spr_params.subtree_cutoff = 0.;

  size_t rf_distance = 1;
  bool skip_intermid_model_opt = false;
  bool impr = true;
  double old_loglh;

  if(converged(cm, loglh, 0.01)) skip_intermid_model_opt = true;

  /* Fast SPR-NNI rounds */
  while( !converged(cm, loglh, 0.01) && rf_distance != 0 && impr) {

    ++iter;
    old_loglh = loglh;
    cm.update_and_write(treeinfo);
    
    // spr round
    LOG_PROGRESS(loglh) << "SPR round " << iter << " (radius: " <<
        spr_params.radius_max << ")" << endl;
    loglh = treeinfo.spr_round(spr_params);
    
    // nni round
    nni(treeinfo, nni_params, loglh);

    if(intermediate_trees.size() == 1){
      
      intermediate_trees.emplace_back(treeinfo.tree());
      rfDist->recalculate_rf();
      rf_distance = rfDist->rf(0,1);
        
    } else {

      intermediate_trees.emplace_back(treeinfo.tree());
      intermediate_trees.erase(intermediate_trees.begin());
      assert(intermediate_trees.size() == 2);

      rfDist->recalculate_rf();
      rf_distance = rfDist->rf(0,1);
    
    }

    if (spr_params.radius_min + radius_step < radius_limit)
    {
      spr_params.radius_min += radius_step;
      spr_params.radius_max += radius_step;
    }

    impr = (loglh - old_loglh > _lh_epsilon);

  }
  
  
  if (!skip_intermid_model_opt)
  {

    cm.update_and_write(treeinfo);

    /* optimize model parameters a bit more thoroughly */
    LOG_PROGRESS(loglh) << "Model parameter optimization (eps = " <<
                                                            interim_modopt_eps << ")" << endl;
    loglh = optimize_model(treeinfo, interim_modopt_eps);

  }

  // slow spr setup
  iter = 0;
  spr_params.thorough = 1;
  spr_params.radius_min = 1;
  spr_params.radius_max = slow_spr_radius;
  spr_params.ntopol_keep = 20;
  spr_params.subtree_cutoff = _spr_cutoff;
  spr_params.reset_cutoff_info(loglh);

  do
  {
    cm.update_and_write(treeinfo);
    ++iter;
    old_loglh = loglh;

    LOG_PROGRESS(old_loglh) << (spr_params.thorough ? "SLOW" : "FAST") <<
        " spr round " << iter << " (radius: " << spr_params.radius_max << ")" << endl;
    loglh = treeinfo.spr_round(spr_params);
    nni(treeinfo, nni_params, loglh);
    
    /* optimize ALL branches */
    loglh = treeinfo.optimize_branches(_lh_epsilon, 1);

    impr = (loglh - old_loglh > _lh_epsilon);
    
  } while (impr);
  
  cm.update_and_write(treeinfo);
  LOG_PROGRESS(loglh) << "Model parameter optimization (eps = " << final_modopt_eps << ")" << endl;
  loglh = optimize_model(treeinfo, final_modopt_eps);

  cm.update_and_write(treeinfo);

  return loglh;
}

double Optimizer::evaluate(TreeInfo& treeinfo, CheckpointManager& cm)
{
  const double fast_modopt_eps = 10.;

  SearchState local_search_state = cm.search_state();
  auto& search_state = ParallelContext::group_master_thread() ? cm.search_state() : local_search_state;
  ParallelContext::barrier();

  double &loglh = search_state.loglh;

  /* Compute initial LH of the starting tree */
  loglh = treeinfo.loglh();

  CheckpointStep resume_step = search_state.step;
  auto do_step = [&search_state,resume_step](CheckpointStep step) -> bool
      {
        if (step >= resume_step)
        {
          search_state.step = step;
          return true;
        }
        else
          return false;;
      };

  if (do_step(CheckpointStep::brlenOpt))
  {
    cm.update_and_write(treeinfo);
    LOG_PROGRESS(loglh) << "Initial branch length optimization" << endl;
    loglh = treeinfo.optimize_branches(fast_modopt_eps, 1);
  }

  /* Model optimization */
  if (do_step(CheckpointStep::modOpt1))
  {
    cm.update_and_write(treeinfo);
    LOG_PROGRESS(loglh) << "Model parameter optimization (eps = " << _lh_epsilon << ")" << endl;
    loglh = optimize_model(treeinfo);
  }

  if (do_step(CheckpointStep::finish))
    cm.update_and_write(treeinfo);

  return loglh;
}

double Optimizer::convergence_rate(CheckpointManager& cm, double test_loglh){
  
  double retval;

  if (cm.best_loglh() == -INFINITY){
    retval = -INFINITY;
  } else {
    if(test_loglh > cm.best_loglh()) retval = 0;
    else {
      retval = (cm.best_loglh() - test_loglh) / fabs(cm.best_loglh());
    }
  }

  return retval;
}

bool Optimizer::first_search_done(CheckpointManager& cm){
  return (cm.best_loglh() != -INFINITY);
}

bool Optimizer::converged(CheckpointManager& cm, double test_loglh, double epsilon){
  if(first_search_done(cm)){
    return convergence_rate(cm,test_loglh) < epsilon ? true : false;
  } 
  else return false;
}

int Optimizer::adaptive_slow_spr_radius(double difficulty){
  if (difficulty <= 0.5){
    return (int) (50*difficulty + 5);
  } else {
    return (int) (-50*difficulty + 55);
  }
}

double Optimizer::sampling_noise_epsilon(corax_treeinfo_t* treeinfo, double** persite_lnl, double loglh, bool rell){

  double epsilon = 0;
  unsigned int total_sites = 0;

  unsigned int partition_count = treeinfo->partition_count;
  unsigned int part, i;

  for(part = 0; part<partition_count; part++) total_sites += treeinfo->partitions[part]->sites;

  if(rell) {

    srand(_seed);
    int pIndex, sIndex, exp;
    int rell_size = 1000;
    double * rell_dist = new double[rell_size];
    double * logl_diff_dist = new double[rell_size];

    for(exp = 0; exp < rell_size; exp++){
      
      double rell_loglh = 0;
      for (i = 0; i<total_sites; i++){
        
        if(partition_count > 1){
          pIndex = rand() % partition_count;
        } else {
          pIndex = 0;
        }

        sIndex =  rand() % treeinfo->partitions[pIndex]->sites;
        rell_loglh += persite_lnl[pIndex][sIndex];
      }

      rell_dist[exp] = rell_loglh;
    }

    for(exp = 0; exp < rell_size; exp++){
      int ind1 = rand() % rell_size;
      int ind2 = rand() % rell_size;

      logl_diff_dist[exp] = fabs(rell_dist[ind1] - rell_dist[ind2]);
    }

    sort(logl_diff_dist, logl_diff_dist + rell_size);

    int epsilon_index = (int) (0.95*rell_size);
    epsilon = logl_diff_dist[epsilon_index];

    delete[] logl_diff_dist;
    delete[] rell_dist;

    //cout << "RELL epsilon " << epsilon << endl;

  } else {

    double stdev = 0;
    double mean = loglh / total_sites;

    for (part = 0; part < partition_count; part++)
      for(i = 0; i < treeinfo->partitions[part]->sites; i++)
        stdev += pow(persite_lnl[part][i] - mean, 2);

    stdev = sqrt(stdev / total_sites);
    
    epsilon = 1.645 * M_SQRT2 * sqrt(total_sites) * stdev;
    //cout << "NO-RELL epsilon " << epsilon << endl;
  }

  return epsilon;
}