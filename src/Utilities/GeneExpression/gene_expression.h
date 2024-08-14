#ifndef UTILITIES_GENEEXPRESSION_GENE_EXPRESSION_H_
#define UTILITIES_GENEEXPRESSION_GENE_EXPRESSION_H_

#include <string>

#include "Utilities/GeneExpression/gene_expression.pb.h"

namespace gene_expression {

class GeneProfile {
  private:
    uint32_t _number_genes;
    int32_t* _rank;

    double _max_positive_connection_strength;

    // Profile _profile;

    std::string _id;

  public:
    GeneProfile();
    ~GeneProfile();

    int Build(const Profile& proto);

    double Association(const GeneProfile& rhs) const;

    const std::string& id() const {
      return _id;
    }
};


}  // namespace gene_expression

#endif // UTILITIES_GENEEXPRESSION_GENE_EXPRESSION_H_
