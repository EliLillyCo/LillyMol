#include <algorithm>
#include <iostream>
#include <string_view>

#include "Foundational/accumulator/accumulator.h"

#include "needle.h"

namespace needle {

using std::cerr;

uint32_t Needle::_number_neighbours = 1;
double* Needle::_max_possible_association = nullptr;
int Needle::_accumulate_matching_genes = 0;

// If gene ids that are above this number (if non zero) then those
// gene ids get placed into a hash rather than to the _gene_to_index array.
uint32_t Needle::_large_gene_id_threshold = 100000;

void
set_accumulate_matching_genes(int s) {
  Needle::_accumulate_matching_genes = s;
}

void set_large_gene_id_threshold(uint32_t s) {
  Needle::_large_gene_id_threshold = s;
}

Needle::Needle() {
  _number_genes = 0;
  _highest_gene_number = 0;
  _gene_to_index = nullptr;
  _sign = nullptr;
  _lhs = nullptr;
  _rank = nullptr;
}

Needle::~Needle() {
  if (_gene_to_index != nullptr) {
    delete [] _gene_to_index;
  }

  delete [] _sign;
  delete [] _lhs;
  delete [] _rank;
}

void
Needle::SetMaxPossibleAssociation(uint32_t maxrank) {

  _max_possible_association = new double[maxrank + 1];

  // This should never be used since this is a division.
  _max_possible_association[0] = 0;

  uint64_t sum = 0;
  for (uint32_t i = 1; i <= maxrank; ++i) {
    uint64_t tmp = i * i;
    sum += tmp;
    _max_possible_association[i] = static_cast<double>(sum);
    // cerr << i << " sum " << sum << '\n';
  }

  cerr << "Needle::SetMaxPossibleAssociation:max possible rank " << maxrank << 
          " max_sum " << _max_possible_association[maxrank - 1] << '\n';

  return;
}

static constexpr int kSwitchToSparse = 100000;

int
Needle::Build(const gene_expression::Profile& proto, uint32_t maxrank) {
  _name = proto.name();
  _number_genes = proto.gene_size();
  _lhs = new int[_number_genes];
  _rank = new int[_number_genes];
  _sign = new int[_number_genes];

  _highest_gene_number = 0;
  for (const gene_expression::Gene& g : proto.gene()) {
    if (g.gene_id() > _highest_gene_number) {
      _highest_gene_number = g.gene_id();
    }
  }
  
  // cerr << " cmp lar " << _large_gene_id_threshold << " high " << _highest_gene_number << '\n';
  if (_large_gene_id_threshold > 0 && _highest_gene_number > _large_gene_id_threshold) {
    _gene_to_index = new int[_large_gene_id_threshold + 1];
    _max_gene_id_in_index = _large_gene_id_threshold;
  } else {
    _gene_to_index = new int[_highest_gene_number + 1];
    _max_gene_id_in_index = _highest_gene_number;
  }

  if (_max_gene_id_in_index > 10000000) {
    cerr << "WARNING: highest gene number " << _highest_gene_number << 
         " use the -maxgeneid option to lower memory requirements\n";
  }

  std::fill_n(_gene_to_index, _max_gene_id_in_index + 1, -1);

  for (int i = 0; i < proto.gene_size(); ++i) {
    uint32_t g = proto.gene(i).gene_id();
    if (g <= _max_gene_id_in_index) {
      _gene_to_index[g] = i;
    } else {
      _large_gene_to_index[g] = i;
    }
    if (proto.gene(i).score() < 0.0) {
      _sign[i] = -1;
    } else {
      _sign[i] = 1;
    }
  }

  return 1;
}

int
Needle::GeneToIndex(uint32_t gene_id) const {
  if (gene_id <= _max_gene_id_in_index) {
    return _gene_to_index[gene_id];
  }
  if (gene_id > _highest_gene_number) {
    return -1;
  }

  auto iter = _large_gene_to_index.find(gene_id);
  if (iter == _large_gene_to_index.end()) {
    return -1;
  }

  return iter->second;
}

template <typename CMP>
int CompareScores(float s1, float s2, CMP) {
  return CMP(std::abs(s1), std::abs(s2));
}

// #define CHECK_SORTED

// `proto` is a member of the haystack.
// needle_max_rank and haystack_max_rank are what were specified on the command line.
int
Needle::Compare(const gene_expression::Profile& proto,
                uint32_t needle_max_rank,
                uint32_t haystack_max_rank) {

  uint32_t istop_rhs = proto.gene_size();
  // cerr << "Haystack has " << proto.gene_size() << " genes, cmp " << haystack_max_rank << '\n';
  if (istop_rhs > haystack_max_rank) {
    istop_rhs = haystack_max_rank;
  }

  uint32_t istop_lhs = _number_genes;
  if (istop_lhs > needle_max_rank) {
    istop_lhs = needle_max_rank;
  }

#ifdef CHECK_SORTED
  cerr << "Check " << istop_lhs << " genes, needle_max_rank " << needle_max_rank << " istop_rhs " << istop_rhs << '\n';
#endif

  // Mark all the genes as uninvolved.
  std::fill_n(_lhs, istop_lhs, 0);

  // Make a list of the genes from the RHS that are being used.

  std::vector<uint32_t> rhs_needed;
  rhs_needed.reserve(istop_rhs);

  for (uint32_t i = 0; i < istop_rhs; ++i) {
    const gene_expression::Gene& g = proto.gene(i);
    uint32_t gene_id = g.gene_id();
    int ndx = GeneToIndex(gene_id);
    if (ndx < 0) {
      continue;
    }

    if (static_cast<uint32_t>(ndx) > needle_max_rank) {
      continue;
    }

    rhs_needed.push_back(i);
    // Mark the left hand side of the comparison.
    _lhs[ndx] = 1;
  }

  // No genes in common within any maxrank settings.
  if (rhs_needed.empty()) {
    cerr << "No genes to be compared\n";
    return 1;
  }

  // Assign Needle ranks.
  int rank = 0;
  for (uint32_t i = 0; i < istop_lhs; ++i) {
    if (_lhs[i] == 0) {
      continue;
    }
    _rank[i] = rank;
    ++rank;
  }

  // Loop through all the haystack genes being processed.
  // Fetch the corresponding item from the lhs, get ranks and compute...

  int64_t sum = 0;
  uint32_t in_common = rhs_needed.size();
  for (uint32_t i = 0; i < in_common; ++i) {
    uint32_t rhs = rhs_needed[i];
    const gene_expression::Gene& gene = proto.gene(rhs);
    uint32_t gene_id = gene.gene_id();
    uint32_t lhs = GeneToIndex(gene_id);

    int64_t rhs_signed_rank;
    if (gene.score() < 0.0) {
      rhs_signed_rank = -static_cast<int64_t>(in_common - i);
    } else {
      rhs_signed_rank = static_cast<int64_t>(in_common - i);
    }
    // cerr << "in_common " << in_common << " i " << i << " rhs_signed_rank " << rhs_signed_rank << '\n';
    assert(_rank[lhs] >= 0);

    int64_t lhs_signed_rank;
    if (_sign[lhs] < 0) {
      lhs_signed_rank = - static_cast<int64_t>(in_common - _rank[lhs]);
    } else {
      lhs_signed_rank = static_cast<int64_t>(in_common - _rank[lhs]);
    }
    sum += lhs_signed_rank * rhs_signed_rank;

    // cerr << i << " rhs " << rhs << " gene_id " << gene_id << " lhs " << lhs << " rank " << lhs_signed_rank << " rhs " << rhs_signed_rank << " sum " << sum << '\n';
  }

  if (sum == 0) {
    return 1;
  }

  if (sum > _max_possible_association[in_common]) {
    cerr << " sum " << sum << " in_common " << in_common << " max " << _max_possible_association[in_common] << '\n';
  }
  float score = (static_cast<double>(sum) / _max_possible_association[in_common]);

#ifdef STORE_MATCHING_GENES
  // Convert rhs_needed from indices into gene ids
  for (uint32_t i = 0; i < rhs_needed.size(); ++i) {
    uint32_t j = rhs_needed[i];
    rhs_needed[i] = proto.gene(j).gene_id();
  }
#endif
  return InsertIntoNbrList(proto.name(), score, in_common, rhs_needed);
}

int
Needle::InsertIntoNbrList(const std::string& name, float score, int in_common,
                const std::vector<uint32_t>& rhs_needed) {

  if (rhs_needed.empty()) {
    cerr << "rhs_needed empty!\n";
  }
#ifdef CHECK_SORTED
  cerr << "sum " << sum << " Storing " << score << " current size " << _nbr_list.size() << '\n';
  cerr << _max_possible_association[in_common] << " max possible, " << in_common << " in common\n";
#endif
  assert(std::abs(score) <= 1.0f);

  if (_nbr_list.empty()) {
    // Need to reserve one extra so that the insert_before call below ramains valid
    // if something is added near the end of a list.
    _nbr_list.reserve(_number_neighbours + 1);
    _nbr_list << new Neighbour(name, score, in_common);
#ifdef STORE_MATCHING_GENES
    _nbr_list.back()->genes = std::move(rhs_needed);
#endif
    return 1;
  } 

  // cerr << "Last item on list " << _nbr_list.back()->score  << '\n';

  // Lower score than the least similar we have. If the list is full
  // just return, otherwise append to end.
  if (std::abs(score) <= std::abs(_nbr_list.back()->score)) {
    if (_nbr_list.size() < _number_neighbours) {
      _nbr_list << new Neighbour(name, score, in_common);
#ifdef STORE_MATCHING_GENES
      _nbr_list.back()->genes = std::move(rhs_needed);
#endif
    }
    return 1;
  }

  // If better than first on the queue, insert at beginning.
  if (std::abs(score) >= std::abs(_nbr_list[0]->score)) {
    if (_nbr_list.size() >= _number_neighbours) {
      _nbr_list.pop();
    }

    _nbr_list.insert_at_beginning(new Neighbour(name, score, in_common));
#ifdef STORE_MATCHING_GENES
    _nbr_list.front()->genes = std::move(rhs_needed);
#endif
    return 1;
  }

  // We need to find it.

  auto cmp = [](const Neighbour* n, float score) {
    return std::abs(n->score) >= std::abs(score);
  };
  auto iter = std::lower_bound(_nbr_list.rawdata(), _nbr_list.rawdata() + _nbr_list.size(),
                   score, cmp);
#ifdef CHECK_SORTED
  for (uint32_t i = 0; i < _nbr_list.size(); ++i) {
    cerr << i << ' ' << _nbr_list[i]->score << '\n';
  }
  cerr << " score " << score << " goes at " << (iter - _nbr_list.rawdata()) << '\n';
#endif

  // Very important that the pop operation be done after the insertion, and not before.
#ifdef STORE_MATCHING_GENES
  Neighbour* tmp = new Neighbour(name, score, in_common);
  tmp->genes = std::move(rhs_needed);
  _nbr_list.insert_before(iter - _nbr_list.rawdata(), tmp);
#else
  _nbr_list.insert_before(iter - _nbr_list.rawdata(), new Neighbour(name, score, in_common));
#endif

  if (_nbr_list.size() > _number_neighbours) {
    _nbr_list.pop();
  }

#ifdef CHECK_SORTED
  CheckSorted();
#endif

  return 1;
}

int
Needle::CheckSorted() const {
  for (uint32_t i = 1; i < _nbr_list.size(); ++i) {
    if (std::abs(_nbr_list[i-1]->score) >= std::abs(_nbr_list[i]->score)) {
      continue;
    }
    cerr << "Needle::CheckSorted:out of order at " << i << ' ' << (_nbr_list[i-1]->score - _nbr_list[i]->score) << '\n';
    for (uint32_t j = 0; j < _nbr_list.size(); ++j) {
      cerr << j << ' ' << _nbr_list[j]->score << '\n';
    }
    return 0;
  }

  return 1;
}

int
Needle::WriteNeighbours(uint32_t max_genes_in_common, IWString_and_File_Descriptor& output) {
  static constexpr char kSep = ' ';

#ifdef STORE_MATCHING_GENES
  static constexpr std::string_view kMissing = "-1";
#endif

  // First Needle to write a neighbour list must write the header.
  static bool first_call = true;
  if (first_call) {
    output << "Needle" << kSep << "Haystack" << kSep << "Score" << kSep << "Common";
#ifdef STORE_MATCHING_GENES
    for (uint32_t i = 0; i < max_genes_in_common; ++i) {
      output << kSep << "gid" << i;
    }
#endif
    output << "\n";

    first_call = false;
  }

  for (uint32_t i = 0; i < _nbr_list.size(); ++i) {
    const Neighbour* n = _nbr_list[i];
    output << _name << kSep << n->name << kSep << n->score << kSep << n->genes_in_common;
#ifdef STORE_MATCHING_GENES
    for (uint32_t gene : n->genes) {
      output << kSep << gene;
    }

    for (uint32_t i = n->genes.size(); i < max_genes_in_common; ++i) {
      output << kSep << kMissing;
    }
#endif
    output << '\n';

    output.write_if_buffer_holds_more_than(4096);
  }

  return output.good();
}

void
Needle::UpdateNnbrStatistics(uint64_t& exact_matches_found, Accumulator<double>& acc_sim) const {
  if (_nbr_list.empty()) {
    return;
  }
  if (abs(_nbr_list.front()->score - 1.0f) < 1.0e-06) {
    ++exact_matches_found;
  }

  _nbr_list.each_lambda([&acc_sim] (const Neighbour* nbr) {
    acc_sim.extra(nbr->score);
  });
}

bool
operator==(const Neighbour& n1, const Neighbour& n2) {
#ifdef DEBUG_COMPARE_NEIGHGOURS
  cerr << "cmp " << n1.name << " and " << n2.name << '\n';
#endif
  if (n1.name != n2.name) {
    return false;
  }

#ifdef DEBUG_COMPARE_NEIGHGOURS
  cerr << " cmp score " << n1.score << " and " << n2.score << '\n';
#endif
  if (n1.score == n2.score) {
  } else if (std::abs(n1.score - n2.score) < 1.0e-05) {
  } else {
    return false;
  }

#ifdef DEBUG_COMPARE_NEIGHGOURS
  cerr << " genes_in_common " << n1.genes_in_common << " and " << n2.genes_in_common << '\n';
#endif
  return n1.genes_in_common == n2.genes_in_common;
}

std::ostream& operator<<(std::ostream& output, const Neighbour& nbr) {
  output << "name: " << nbr.name << " score: " << nbr.score << " common: " << nbr.genes_in_common;
  return output;
}

uint32_t
Needle::MaxNumberGenes() const {
  uint32_t rc = 0;
  for (const Neighbour* nbr : _nbr_list) {
    if (nbr->genes.size() > rc) {
      rc = nbr->genes.size();
    }
  }

  return rc;
}

} // namespace needle
