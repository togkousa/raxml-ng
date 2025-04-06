#ifndef RAXML_PARTITIONEDMSA_HPP_
#define RAXML_PARTITIONEDMSA_HPP_

#include "PartitionInfo.hpp"

class AutoPartitioner;
class CVPartitionedMSA;

class PartitionedMSA
{
public:
  PartitionedMSA();
  PartitionedMSA(const NameList& taxon_names);

  // copy/move constructors and assignments
  PartitionedMSA (PartitionedMSA&& other) = default;
  PartitionedMSA& operator=(PartitionedMSA&& other);

  // getters
  const std::vector<std::string>& full_msa_sequences() const { return _full_msa.sequences(); }
  const MSA& full_msa() const { return (part_count() == 1) ? _part_list.at(0).msa() : _full_msa; };
  const MSA& part_msa(size_t index) const { return _part_list.at(index).msa(); };
  const PartitionInfo& part_info(size_t index) const { return _part_list.at(index); };
  const Model& model(size_t index) const { return _part_list.at(index).model(); };
  ModelCRefMap models() const;
  const std::vector<PartitionInfo>& part_list() const { return _part_list; };
  std::vector<PartitionInfo>& part_list() { return _part_list; };
  const NameList& taxon_names()  const { return _taxon_names; };
  const NameIdMap& taxon_id_map() const { return _taxon_id_map; }

  size_t full_msa_site(size_t index, size_t site) const;
  const uintVector& site_part_map() const;
  IdPairVector full_to_parted_sitemap() const;

  size_t taxon_count() const { return _taxon_names.size(); };
  size_t part_count() const { return _part_list.size(); };
  size_t total_sites() const;
  size_t total_patterns() const;
  size_t total_length() const;

  size_t total_free_model_params() const;

  /* given in elements (NOT in bytes) */
  size_t taxon_clv_size() const;

  // setters
  void full_msa(MSA&& msa);
  void part_msa(size_t index, MSA&& msa) { _part_list.at(index).msa(std::move(msa)); };
  void part_msa(size_t index, const corax_msa_t * pll_msa)
  {
    _part_list.at(index).msa(MSA(pll_msa));
  };
  void model(size_t index, Model&& m) { _part_list.at(index).model(std::move(m)); };
  void model(size_t index, const Model& m) { _part_list.at(index).model(m); };

  double difficulty_score() const { return _difficulty_score; }
  void difficulty_score(double score) { _difficulty_score = score; }

  // operations
  void init_single_model(DataType data_type, const std::string &model_string);
  void append_part_info(PartitionInfo&& part_info) { _part_list.push_back(std::move(part_info)); };

  template <class... Args>
  void emplace_part_info (Args&&... args)
  {
    _part_list.emplace_back(std::forward<Args>(args)...);
  }

  void split_msa();
  void compress_patterns(bool store_backmap = false);
  
  /* Cross validation apporach */
  void split_msa_cross_validation(const Options& opts,
                                  double split_radio, 
                                  unsigned int seed);

  void set_model_empirical_params();

  CVPartitionedMSA& parted_training_msa() { return *(_parted_training_msa.get()); }
  CVPartitionedMSA& parted_testing_msa() { return *(_parted_testing_msa.get()); }
  
  std::shared_ptr<CVPartitionedMSA> parted_training_msa_shared_ptr() { return _parted_training_msa; }
  std::shared_ptr<CVPartitionedMSA> parted_testing_msa_shared_ptr() { return _parted_testing_msa; }

protected:
  NameList _taxon_names;
  NameIdMap _taxon_id_map;
  mutable uintVector _site_part_map;

  void copy_taxon_names(const NameList& taxon_names);

private:
  std::vector<PartitionInfo> _part_list;
  std::shared_ptr<AutoPartitioner> _auto_part;
  
  MSA _full_msa; // why is this a stack object and not a pointer to heap?

  std::shared_ptr<CVPartitionedMSA> _parted_training_msa;
  std::shared_ptr<CVPartitionedMSA> _parted_testing_msa;
  
  //std::shared_ptr<MSA> _training_msa;
  //std::shared_ptr<MSA> _testing_msa;
  
  
  double _difficulty_score;

  IDVector _training_sites_map;
  IDVector _testing_sites_map;

  uintVector get_site_part_assignment() const;
  void set_taxon_names(const NameList& taxon_names);
};

class AutoPartitioner
{
public:
  void init_from_string(PartitionedMSA& part_msa, DataType data_type, const std::string &model_string);
  void update_partition_ranges(PartitionedMSA& part_msa);

private:
  doubleVector get_column_entropies(const PartitionedMSA& part_msa);
  std::string resolve_auto_range(const doubleVector& col_entropies, size_t part_num, double binw);
};

class CVPartitionedMSA : public PartitionedMSA {

  public:
    CVPartitionedMSA(const NameList& taxon_names, 
                    const uintVector& site_part_map) : PartitionedMSA() 
    {
      copy_taxon_names(taxon_names);
      set_part_map(site_part_map);
      //init_part_info(opts, site_part_map);
    }
    
  void init_part_info(const Options& opts);
  
  private:
    void set_part_map(const uintVector& site_part_map) { _site_part_map = site_part_map; }
    void correct_range_strs();
};

std::ostream& operator<<(std::ostream& stream, const PartitionedMSA& part_msa);

#endif /* RAXML_PARTITIONEDMSA_HPP_ */
