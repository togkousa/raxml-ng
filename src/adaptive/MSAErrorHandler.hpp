#ifndef RAXML_ADAPTIVE_MSAERRORHANDLER_HPP_
#define RAXML_ADAPTIVE_MSAERRORHANDLER_HPP_

#include "../PartitionedMSA.hpp"
#include "../MSA.hpp"
#include "../Tree.hpp"

#include <vector>
#include <map>
#include <tuple>
#include <math.h>
#include <memory.h>
#include <random>
#include "../TreeInfo.hpp"

using namespace std;

class MSAErrorHandler{

    private:
        
        double msa_error_rate;
        unsigned int seed;

        unsigned int partition_count;
        unsigned int *n_rates;
        unsigned int *n_states;
        unsigned int *n_states_padded;
        unsigned int *n_sites;

        unsigned int mutations;
        unsigned int dist_size;
        
        std::vector<int*> *mutations_info;
        std::vector<double*> *original_clvs;
        double* delta_loglh_dist;

        std::string outfile;

        void init_random_seed();

        void apply_single_mutation(double *clv,
                                    unsigned int states,
                                    unsigned int states_padded,
                                    unsigned int rates);
        
        void reverse_single_mutation(double* clv,
                                    double* original_clv,
                                    unsigned int states_padded,
                                    unsigned int rates);
        
        void clear_vectors();

        void write_dist_to_file();



    public:

        // constructors
        MSAErrorHandler(const corax_treeinfo_t* treeinfo,
                        double _msa_error_rate, 
                        unsigned int _seed,
                        std::string _outfile);
                        
        ~MSAErrorHandler();

        void apply_random_mutations(const corax_treeinfo_t* treeinfo);
        void reverse_mutations(const corax_treeinfo_t* treeinfo);

        void set_partition_info(unsigned int part_index, 
                                unsigned int rates, 
                                unsigned int states,
                                unsigned int states_padded, 
                                unsigned int sites);


        double* generate_delta_loglh_dist(TreeInfo& treeinfo, unsigned int _dist_size, double init_loglh);

};


#endif /* RAXML_ADAPTIVE_MSAERRORHANDLER_HPP_ */
