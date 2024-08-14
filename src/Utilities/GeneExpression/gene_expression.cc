#include "gene_expression.h"

namespace gene_expression {

GeneProfile::GeneProfile() {
  _number_genes = 0;
  _rank = nullptr;
  _max_positive_connection_strength = 0;
}

GeneProfile::~GeneProfile() {
  _number_genes = 0;
  delete [] _rank;
}

double
GeneProfile::Association(const GeneProfile& rhs) const {
  assert(_number_genes == rhs._number_genes);

  int64_t sum = 0;

  // could use std::inner_product, but why...
  for (uint32_t i = 0; i < _number_genes; ++i) {
    sum += _rank[i] * rhs._rank[i];
  }

  return static_cast<double>(sum) / _max_positive_connection_strength;
}

int
GeneProfile::Build(const Profile& proto) {
  _id = proto.name();

#ifdef REIMPLEMENT_QQ
  _number_genes = proto.rank_size();

  _rank = new int[_number_genes];

  for (uint32_t i = 0; i < _number_genes; ++i) {
    _rank[i] = proto.rank(i);
  }

  uint64_t sum = 0;

  for (uint32_t i = 0; i < _number_genes; ++i) {
    sum += _rank[i] * _rank[i];
  }

  _max_positive_connection_strength = static_cast<double>(sum);
#endif

  return 1;
}

}  // namespace gene_expression
