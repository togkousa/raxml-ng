#ifndef RAXML_ADAPTIVE_KHTESTER_HPP_
#define RAXML_ADAPTIVE_KHTESTER_HPP_

#include <corax/corax.h>
#include "../TreeInfo.hpp"

using namespace std;

class KHtester{

    private:
        unsigned int partition_count;
        unsigned int* sites;
        unsigned int seed;
        double kh_epsilon;

        void init_random_seed();

    public:
        KHtester(const corax_treeinfo_t* treeinfo, unsigned int _seed, double **lnl1, double **lnl2);

        void init_persite_lnl_pointers(double **lnl1, double **lnl2);
        void delete_persite_lnl_pointers(double **lnl1, double **lnl2);
        
        ~KHtester();

        double kh_test(double **persite_lnl_1, 
                        double **persite_lnl_2,
                        unsigned int nBootstrap);

        
};


#endif /* RAXML_ADAPTIVE_KHTESTER_HPP_ */
