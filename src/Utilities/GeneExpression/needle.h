#ifndef UTILITIES_GENEEXPRESSION_NEEDLE_H_
#define UTILITIES_GENEEXPRESSION_NEEDLE_H_

#include <iostream>
#include <string>

#include "absl/container/flat_hash_map.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwstring/iwstring.h"

#include "Utilities/GeneExpression/gene_expression.pb.h"

// Set this in order to accumulate with each neighbour the
// list of matched genes.
#define STORE_MATCHING_GENES

namespace needle {

struct Neighbour {
  public:
    std::string name;
    float score;
    // genes_in_common will be the same as genes.size()
    // but we may not always accumulate the matching genes.
    uint32_t genes_in_common;
    // A list of the gene ids shared.
    std::vector<uint32_t> genes;

  public:
};

struct NeedleCompare {
  bool operator()(const Neighbour& n1, const Neighbour& n2) {
    return n1.score < n2.score;
  }
};

class Needle {
  private:
    // The number of genes we read.
    uint32_t _number_genes;
    // We read a bunch of gene ids from our proto. We record the highest number encountered
    // so we can quickly figure out if we have data for a particular gene.
    uint32_t _highest_gene_number;

    //  This is a mapping from gene id
    // to where that gene occurs in our proto. A negative entry means not present.
    int* _gene_to_index;

    // Gene ids larger than this threshold are stored in a hash rather than the
    // _gene_to_index array. That is because gene_id's can be quite large and
    // that consumes a lot of memory.
    static uint32_t _large_gene_id_threshold;

    // The size of the _gene_to_index array can vary, so store the highest
    // gene number that is in that array. Genes above that number will need
    // to go to _large_gene_to_index.
    uint32_t _max_gene_id_in_index;

    // For gene id's above _max_gene_id_in_index;
    absl::flat_hash_map<uint32_t, int> _large_gene_to_index;

    // Observe that we can get very high gene id's (100M) which will cause memory
    // problems if we are using a large number of needles. We can have a sparse
    // index for dealing with large gene id's. TODO:ianwatson figure this out.

    // For each gene, the sign of the response.
    int* _sign;

    // When comparing with a haystack gene expression, this keeps track of
    // which of our genes are to be used in the computation.
    // This is a mapping from an index over the haystack member to
    // our genes.
    int* _lhs;

    // Once we have decided which items are needed, we scan the list and
    // sequentially assign a rank to each.
    int* _rank;

    // Each item will be assigned a signed rank, depending on whether or not the
    // score is positive or negative.
    int* _signed_rank;

    // std::priority_queue<Neighbour, std::vector<Neighbour>, NeedleCompare> _neighbours;

    resizable_array_p<Neighbour> _nbr_list;

    // Every class instance needs to know the number of neighbours to retain.
  public:
    static uint32_t _number_neighbours;
  private:

    // Every class instance needs to know the maximum possible association value for
    // the number of genes found to be in common.
    static double* _max_possible_association;

    // Whether or not we keep track of the matching gene ids during comparisons.
    static int _accumulate_matching_genes;

    IWString _name;

  // Private functions
    int CheckSorted() const;
    int GeneToIndex(uint32_t gene_id) const;
    int InsertIntoNbrList(const std::string& name, float score, int in_common,
                const std::vector<uint32_t>& rhs_needed);

  public:
    Needle();
    ~ Needle();

    int Build(const gene_expression::Profile& proto, uint32_t maxrank);

    void set_name(IWString& s) {
      _name = s;
    }

    uint32_t number_genes() const {
      return _number_genes;
    }

    void SetMaxPossibleAssociation(uint32_t maxrank);

    int Compare(const gene_expression::Profile& proto,
                uint32_t needle_max_rank,
                uint32_t haystack_max_rank);
      
    void UpdateNnbrStatistics(uint64_t& exact_matches_found, Accumulator<double>& acc_sim) const;

    // Scan the _neighbours list and return the max number of genes_in_common
    // among those neighbours.
    uint32_t MaxNumberGenes() const;

    // max_genes_in_common is only used if we are also writing the matching gene id's.
    // It is simple to compute so it is always included, but may not get used.
    int WriteNeighbours(uint32_t max_genes_in_common, IWString_and_File_Descriptor& output);

    int number_neighbours() const {
      return _nbr_list.number_elements();
    }
    const Neighbour* nbr(int ndx) const {
      return _nbr_list[ndx];
    }

    // Now many gene ids are in the overflow hash.
    uint32_t GeneIdsInHash() const {
      return _large_gene_to_index.size();
    }

    friend void set_large_gene_id_threshold(uint32_t);
    friend void set_accumulate_matching_genes(int);

};

// control whether or not the list of matching gene id's is accumulated and
// written. This will have performance implications, since keeping track of
// the matching gene id's will be expensive.
void set_accumulate_matching_genes(int s);

// Set the threshold above which genes id's are put into a hash rather
// than the _gene_to_index array.
void set_large_gene_id_threshold(uint32_t s);

bool operator==(const Neighbour& n1, const Neighbour& n2);

std::ostream& operator<<(std::ostream& output, const Neighbour& nbr);

}  // namespace needle


#endif  // UTILITIES_GENEEXPRESSION_NEEDLE_H_
