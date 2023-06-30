#include "KHtester.hpp"

KHtester::KHtester(const corax_treeinfo_t* treeinfo, unsigned int _seed, double **lnl1, double **lnl2){

    partition_count = treeinfo->partition_count;

    sites = new unsigned int[partition_count];

    for(unsigned int p = 0; p < partition_count; p++) 
        sites[p] = treeinfo->partitions[p]->sites;

    kh_epsilon = 0;

    seed = _seed;
    init_random_seed();

    init_persite_lnl_pointers(lnl1, lnl2);

}

KHtester::~KHtester(){

    delete[] sites;
}

void KHtester::init_random_seed(){
    if(ParallelContext::group_master_thread()) srand(seed);
    ParallelContext::barrier();
}

void KHtester::init_persite_lnl_pointers(double **lnl1, double **lnl2){

    lnl1 = new double*[partition_count];
    lnl2 = new double*[partition_count];

    for (unsigned int part = 0; part<partition_count; part++){
      lnl1[part] = new double[sites[part]];
      lnl2[part] = new double[sites[part]];
    }
}

void KHtester::delete_persite_lnl_pointers(double **lnl1, double **lnl2){

    for (unsigned int part = 0; part<partition_count; part++){
      delete[] lnl1[part];
      delete[] lnl2[part];
    }

    delete[] lnl1;
    delete[] lnl2;
}

double KHtester::kh_test(double **persite_lnl_1, 
                        double **persite_lnl_2,
                        unsigned int nBootstrap)
{

    unsigned int i,j;

    // calculating likelihoods
    double lnl1=0, lnl2=0;

    for (i=0; i<partition_count; i++){
        for (j=0; j<sites[i]; j++){
            lnl1 += persite_lnl_1[i][j];
            lnl2 += persite_lnl_2[i][j];
        }
    }

    // calculating total len
    unsigned int nSites = 0;
    for (i = 0; i<partition_count; i++) nSites += sites[i];

     // delta Lnl - Improvement
    double Delta_Lnl = (lnl1 - lnl2) >= 0 ? (lnl1 - lnl2) : 0; // to avoid numerical errors
    
    // RELL bootstrap process
    double Delta_Ls[nBootstrap];
    //double mean_Dl = 0;
    int KHsupport=0;

    for (i = 0; i<nBootstrap; i++){
        
        int pIndex, sIndex;
        lnl1=0;
        lnl2=0;

        for(j = 0; j<nSites; j++){

            if(partition_count > 1){
                pIndex = rand() % partition_count;
            } else {
                pIndex = 0;
            }

            sIndex =  rand() % sites[pIndex];
            lnl1 += persite_lnl_1[pIndex][sIndex];
            lnl2 += persite_lnl_2[pIndex][sIndex];
        }

        Delta_Ls[i] = lnl1 - lnl2 - Delta_Lnl;
        if ( Delta_Ls[i] >= Delta_Lnl ) KHsupport++;
        // mean_Dl += Delta_Ls[i] / nBootstrap;
    }

    /* 
     // substracting the mean and counting the number of samples worse that Delta_Lnl
    for(i=0; i<nBootstrap; i++) {
        
        Delta_Ls[i] = Delta_Ls[i] - mean_Dl;
        if ( Delta_Ls[i] >= Delta_Lnl ) KHsupport++;
    } 
    */

    double p_value = (double) KHsupport / nBootstrap ;

    std::sort(Delta_Ls, Delta_Ls + nBootstrap);
    int epsilon_index = (int) (0.95*nBootstrap);
    kh_epsilon = Delta_Ls[epsilon_index] > 0 ? Delta_Ls[epsilon_index] : 0;
    
    return p_value;
}

