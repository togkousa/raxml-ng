#include "MSAErrorHandler.hpp"
#include "../ParallelContext.hpp"
#include <fstream>
#include <cstdio>

MSAErrorHandler::MSAErrorHandler(const corax_treeinfo_t* treeinfo,
                                double _msa_error_rate, 
                                unsigned int _seed,
                                std::string _outfile){
    
    msa_error_rate = _msa_error_rate;
    seed = _seed;
    partition_count = treeinfo->partition_count;
    
    mutations = 0;
    max_loglh = 0;

    n_rates = new unsigned int[partition_count];
    n_states = new unsigned int[partition_count];
    n_states_padded = new unsigned int[partition_count];
    n_sites = new unsigned int[partition_count];

    for(unsigned int p = 0; p<partition_count; p++){
        set_partition_info(p,
                        treeinfo->partitions[p]->rate_cats,
                        treeinfo->partitions[p]->states,
                        treeinfo->partitions[p]->states_padded,
                        treeinfo->partitions[p]->sites);
    }
    
    mutations_info = new std::vector<int *>;
    original_clvs = new std::vector<double *>;

    tmp_brlens = new double*[partition_count];
    original_brlens = new double*[partition_count];

    init_random_seed();

    epsilon = 0;
    outfile_initial = _outfile.length() > 0 ? _outfile + "_initial.csv" : "" ;
    outfile_final = _outfile.length() > 0 ? _outfile + "_final.csv" : "" ;
}

MSAErrorHandler::~MSAErrorHandler(){

    clear_vectors();
    delete[] n_rates;
    delete[] n_states;
    delete[] n_states_padded;
    delete[] n_sites;

    if(delta_loglh_dist) delete[] delta_loglh_dist;
    if(new_loglh_dist) delete[] new_loglh_dist;
    if(errors) delete[] errors;

    if(max_loglh){
        for(unsigned int p = 0; p<partition_count; p++){
            delete[] tmp_brlens[p];
            delete[] original_brlens[p];
        }
    }
    
    delete[] tmp_brlens;
    delete[] original_brlens;
}

void MSAErrorHandler::clear_vectors(){
    
    size_t vec_size;

    // clear mutations_info
    if (mutations_info){
        vec_size = mutations_info->size();
        for (size_t i = 0; i<vec_size; i++) 
            if (mutations_info->at(i)) delete[] mutations_info->at(i);

        mutations_info->clear();
        delete mutations_info;
    }

    // clear original_clv
    if(original_clvs){
        vec_size = original_clvs->size();
        for (size_t i = 0; i<vec_size; i++) 
            if (original_clvs->at(i)) delete[] original_clvs->at(i);
        
        original_clvs->clear();
        delete original_clvs;
    }
}


void MSAErrorHandler::set_partition_info(unsigned int part_index,
                                        unsigned int rates,
                                        unsigned int states,
                                        unsigned int states_padded,
                                        unsigned int sites)
{
    n_rates[part_index] = rates;
    n_states[part_index] = states;
    n_states_padded[part_index] = states_padded;
    n_sites[part_index] = sites;
}

void MSAErrorHandler::init_random_seed(){
    /* if(ParallelContext::group_master_thread()) srand(seed);
    ParallelContext::barrier(); */
    srand(seed);
}

void MSAErrorHandler::apply_single_mutation(double *clv,
                                            unsigned int states,
                                            unsigned int states_padded,
                                            unsigned int rates)
{
    
    bool is_gap;
    unsigned int i, random_integer;
    for(i = 0; i < states; i++)
        if (clv[i] == 0) break;
    
    is_gap = i == states ? true : false;
    if(is_gap) {
        return; // we don't apply mutations to a gap
    }
    
    /* 
    cout << "-------------------------------------------" << endl;
    for(i = 0; i < states_padded*rates; i++)
        cout << "Original clv[" << i << "] = " << clv[i] << endl;
    */

    double *original_clv = new double[states_padded];
    memcpy(original_clv, clv, states_padded*sizeof(double));
    original_clvs->push_back(original_clv);
    
    memset(clv, 0, states_padded*rates*sizeof(double));
    
    do{
        random_integer = rand() % states;
    } while (original_clv[random_integer] == 1);
    

    clv[random_integer] = 1;

    // naive - I can just set 1 to the other positions
    for(i = 1; i<rates; ++i){
        memcpy(clv + i*states_padded, clv, states_padded*sizeof(double));
        //clv += states_padded;
    }
    
    /* 
    cout << "-------------------------------------------" << endl;
    for(i = 0; i < states_padded*rates; i++)
        cout << "New clv[" << i << "] = " << clv[i] << endl;
    
    getchar();
    */

    mutations++;
}

void MSAErrorHandler::apply_random_mutations(const corax_treeinfo_t* treeinfo){

    // declarations
    unsigned int clv_index, rates, states, states_padded, sites, effective_size, i, p;
    double rand_number, _mut;
    corax_partition *partition;

    unsigned int tip_count = treeinfo->tip_count;
    unsigned int partition_count = treeinfo->partition_count;

    for (p = 0; p < partition_count; ++p){
        
        partition = treeinfo->partitions[p];
        rates = n_rates[p];
        states = n_states[p];
        states_padded = n_states_padded[p];
        sites = n_sites[p];
        effective_size = states_padded * rates;

        for (i = 0; i < tip_count; ++i){
            
            clv_index = treeinfo->tree->nodes[i]->clv_index;
            double *clv = partition->clv[clv_index];
            
            for (unsigned int s = 0; s < sites; s++){
                
                _mut = mutations;
                rand_number = ((double) rand() / (RAND_MAX));
                
                if (rand_number < msa_error_rate){ // then insert random mutation
                    
                    // setup mutations info
                    apply_single_mutation(clv, states, states_padded, rates);
                    
                    if (_mut < mutations){
                        int *mut_info = new int[3];
                        mut_info[0] = (int) p;
                        mut_info[1] = (int) clv_index;
                        mut_info[2] = (int) s;
                        mutations_info->push_back(mut_info);
                    }
                }
                clv += effective_size;
                
            }
            
        }
    }
    
    /* if mutation size == 0, apply single mutation */
    while (mutations == 0){
        
        _mut = mutations;
        p = partition_count > 1 ? ( rand() % partition_count ) : 0;
        partition = treeinfo->partitions[p];

        rates = n_rates[p];
        states = n_states[p];
        states_padded = n_states_padded[p];
        sites = n_sites[p];
        effective_size = states_padded * rates;

        unsigned int tip = rand() % tip_count;
        clv_index = treeinfo->tree->nodes[tip]->clv_index;
        unsigned int s = rand() % sites;

        double *clv = partition->clv[clv_index];
        clv += effective_size * s;
        apply_single_mutation(clv, states, states_padded, rates);

        if (_mut < mutations){
            int *mut_info = new int[3];
            mut_info[0] = (int) p;
            mut_info[1] = (int) clv_index;
            mut_info[2] = (int) s;
            mutations_info->push_back(mut_info);
        }
        //cout << "Mutations info size = " << mutations_info->size() << endl;
    }

}

void MSAErrorHandler::reverse_single_mutation(double* clv,
                                            double* original_clv,
                                            unsigned int states_padded,
                                            unsigned int rates)
{

    memcpy(clv, original_clv, states_padded*sizeof(double));

    for(unsigned int i = 1; i<rates; ++i){
        memcpy(clv + i*states_padded, clv, states_padded*sizeof(double));
    }

    mutations--;
}


void MSAErrorHandler::reverse_mutations(const corax_treeinfo_t* treeinfo){

    unsigned int rates, states_padded, effective_size;
    corax_partition *partition;

    assert(mutations == mutations_info->size());

    for (unsigned int i = 0; i < mutations_info->size(); i++){
        int p = mutations_info->at(i)[0];
        int clv_index = mutations_info->at(i)[1];
        int s = mutations_info->at(i)[2];

        partition = treeinfo->partitions[p];
        rates = n_rates[p];
        states_padded = n_states_padded[p];
        effective_size = states_padded * rates;

        double *clv = partition->clv[clv_index];
        clv += s * effective_size;
        double *original_clv = original_clvs->at(i);

        reverse_single_mutation(clv, original_clv, states_padded, rates);

        delete[] original_clvs->at(i);
        delete[] mutations_info->at(i);
    }

    original_clvs->clear();
    mutations_info->clear();

    assert(mutations == 0);

}

void MSAErrorHandler::msa_error_dist(TreeInfo& treeinfo, 
                                    unsigned int _dist_size, 
                                    double init_loglh,
                                    bool initial)
{
    
    double new_loglh;
    dist_size = _dist_size;
    unsigned int errors_size = 10*dist_size;

    // if(delta_loglh_dist) delete[] delta_loglh_dist;

    delta_loglh_dist = new double[dist_size];
    new_loglh_dist = new double[dist_size];
    errors = new double[10*dist_size];

    for (unsigned int experiment = 0; experiment < dist_size; experiment++){
        
        // if(ParallelContext::group_master_thread())
        apply_random_mutations(&treeinfo.pll_treeinfo());
        // ParallelContext::barrier();

        new_loglh = treeinfo.loglh();
        
        double delta_loglh = new_loglh - init_loglh;
        delta_loglh_dist[experiment] = delta_loglh;
        new_loglh_dist[experiment] = new_loglh;

        epsilon += delta_loglh / dist_size;

        // if(ParallelContext::group_master_thread())
        // cout << "Experiment = " << experiment <<", Mutations = " << mutations <<", New loglh = " << delta_loglh << endl;
        
        //getchar();

        // if(ParallelContext::group_master_thread()) 
        // if(ParallelContext::group_master_thread())
        reverse_mutations(&treeinfo.pll_treeinfo());
        
        //ParallelContext::barrier();
        new_loglh = treeinfo.loglh();
        assert(fabs(init_loglh - new_loglh)<1e-4);
        
    }

    assert(delta_loglh_dist != nullptr);
    
    if(outfile_initial.length() > 0) 
        write_dist_to_file(initial ? outfile_initial : outfile_final);
    
    // if(ParallelContext::group_master_thread())
    
    //std::sort(delta_loglh_dist, delta_loglh_dist + dist_size);
    //max_loglh = init_loglh + delta_loglh_dist[dist_size - 1];
    // cout << "Ln-1 " << init_loglh << " and max_loglh " << max_loglh << endl;

    
    // ParallelContext::global_barrier();
    
    /* for (unsigned int rep = 0; rep < errors_size; rep++){

        unsigned int index_1 = 0, index_2 = 0;
        
        while(index_1 == index_2){
            index_1 = rand() % dist_size;
            index_2 = rand() % dist_size;
        }

        errors[rep] = delta_loglh_dist[index_1] - delta_loglh_dist[index_2]; 
    } */

    //int epsilon_index = (int) (0.95*errors_size);
    //epsilon = fabs(epsilon);
    
    return;
}

void MSAErrorHandler::write_dist_to_file(std::string outfile){

    const char* filename = outfile.c_str();
    remove(filename);

    ofstream myfile (filename);
    if (myfile.is_open())
    {
        for(unsigned int count = 0; count < dist_size; count ++){
            myfile << delta_loglh_dist[count] << "\n" ;
        }
        myfile.close();
    }
}



void MSAErrorHandler::store_brlens(corax_treeinfo_t* treeinfo, bool preultimate){

    for(unsigned int p = 0; p<partition_count; p++)
        store_brlens_partition(treeinfo, p, preultimate);

}


void MSAErrorHandler::store_brlens_partition(corax_treeinfo_t* treeinfo, unsigned int partition_id, bool preultimate){

    // cout << "Store with preultimate " << preultimate << endl;

    double **store_matrix = preultimate ? tmp_brlens : original_brlens;
    corax_partition_t* partition = treeinfo->partitions[partition_id];

    unsigned int num_branches = 2*treeinfo->tip_count - 3;
    store_matrix[partition_id] = new double[num_branches];

    memcpy(store_matrix[partition_id], treeinfo->branch_lengths[partition_id], num_branches*sizeof(double));

}

void MSAErrorHandler::set_brlens(corax_treeinfo_t* treeinfo, bool preultimate){

    if(preultimate) store_brlens(treeinfo, false);

    double **copy_pmatrx = preultimate ? tmp_brlens : original_brlens;

    unsigned int num_branches = 2*treeinfo->tip_count - 3;

    for(unsigned int p = 0; p < partition_count; p++){
        memcpy(treeinfo->branch_lengths[p], copy_pmatrx[p], num_branches*sizeof(double));
    }
} 