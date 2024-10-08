#include <cmath>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Utilities/GeneExpression/gene_expression.pb.h"

#include "needle.h"

namespace {

using needle::Needle;
using needle::Neighbour;

TEST(TestNeedle, TestIdentical) {
  std::string string_proto1 = R"pb(name: "t1",
gene {
  gene_id: 10
  score: 1.0
}
gene {
  gene_id: 12
  score: 0.5
}
)pb";

  gene_expression::Profile proto1;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto1, &proto1));

  // This test is not about maxrank, so use a number larger than the number of genes.
  constexpr uint32_t kLarge = 10;

  Needle needle;
  ASSERT_TRUE(needle.Build(proto1, kLarge));

  needle.SetMaxPossibleAssociation(kLarge);
  needle::Needle::_number_neighbours = kLarge;
  
  needle.Compare(proto1, kLarge, kLarge);

  EXPECT_EQ(needle.number_neighbours(), 1);

  Neighbour expected{"t1", 1.0f, 2};
  EXPECT_EQ(*needle.nbr(0), expected);

  needle::delete_max_possible_association();
}

TEST(TestNeedle, TestIdenticalOppositeSign) {
  std::string string_proto1 = R"pb(name: "t1",
gene {
  gene_id: 10
  score: 1.0
}
gene {
  gene_id: 12
  score: 0.5
}
)pb";

  gene_expression::Profile proto1;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto1, &proto1));

  // This test is not about maxrank, so use a number larger than the number of genes.
  constexpr uint32_t kLarge = 10;

  Needle needle;
  ASSERT_TRUE(needle.Build(proto1, kLarge));

  needle.SetMaxPossibleAssociation(kLarge);
  needle::Needle::_number_neighbours = kLarge;
  
  proto1.mutable_gene(0)->set_score(-1.0);
  proto1.mutable_gene(1)->set_score(-0.5);
  needle.Compare(proto1, kLarge, kLarge);

  EXPECT_EQ(needle.number_neighbours(), 1);

  Neighbour expected{"t1", -1.0f, 2};
  EXPECT_EQ(*needle.nbr(0), expected);

  needle::delete_max_possible_association();
}

// Hmmm, made a mistake here. This struct cannot explore
// mutiple neighbours since there are only 2 string protos.
// We should have needle_string_proto and a vector of
// haystack_string_proto. But discovered that too late
// and did not want to change things.
// The neighbour list handling has been tested in the
// neighbour finding code. These tests really just test
// the score.
struct Cmp {
  std::string string_proto1;
  std::string string_proto2;
  uint32_t max_association;
  int number_neighbours;
  uint32_t needle_maxrank;
  uint32_t haystack_maxrank;
  std::vector<Neighbour> expected;
};

class TestGeneExpression : public testing::TestWithParam<Cmp> {
  protected:
    gene_expression::Profile _proto1;
    gene_expression::Profile _proto2;

    Needle _n1;
    Needle _n2;

  // Protected Functions
    void SetUp();
    void TearDown();
};

void
TestGeneExpression::SetUp() {
  needle::initialise_max_possible_association();
}
void
TestGeneExpression::TearDown() {
  needle::delete_max_possible_association();
}

TEST_P(TestGeneExpression, Test1) {
  const auto params = GetParam();
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.string_proto1, &_proto1));
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.string_proto2, &_proto2));

  ASSERT_TRUE(_n1.Build(_proto1, params.needle_maxrank));
  ASSERT_TRUE(_n2.Build(_proto2, params.haystack_maxrank));

  _n1.SetMaxPossibleAssociation(params.max_association);
  needle::Needle::_number_neighbours = params.number_neighbours;

  _n1.Compare(_proto2, params.needle_maxrank, params.haystack_maxrank);

  EXPECT_EQ(_n1.number_neighbours(), params.expected.size());
  if (_n1.number_neighbours() > 0) {
    for (int i = 0; i < _n1.number_neighbours(); ++i) {
      EXPECT_EQ(*_n1.nbr(i), params.expected[i]);
    }
  }
}
INSTANTIATE_TEST_SUITE_P(TestOrder, TestGeneExpression, testing::Values(
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 1
    score: 1.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 1.0, 1}}
},

  // No genes in common
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 2
    score: 1
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{}
},

  // extra genes in the haystack.
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 4
    score: 1
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 1.0f, 1}}
},

  // extra genes in the haystack, diff sign
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: -1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 2
    score: 0.3
  }
  gene {
    gene_id: 55
    score: 0.3
  }
  gene {
    gene_id: 1
    score: 2.0
  }
  gene {
    gene_id: 4
    score: 3.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", -1.0f, 1}}
},

  // extra genes in the needle
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 8
    score: -1
  }
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 3
    score: -1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 1
    score: 2.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 1.0f, 1}}
},

  // same two genes, just order swapped
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 3
    score: -1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 2
    score: 2.0
  }
  gene {
    gene_id: 1
    score: 3.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 0.8f, 2}}
},

  // same two genes, order swapped and extra needle genes
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 9
    score: 1
  }
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 11
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 13
    score: -1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 2
    score: 2.0
  }
  gene {
    gene_id: 1
    score: 3.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 0.8f, 2}}
},

  // same two genes, order swapped and exra haystack genes
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 3
    score: -1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 53
    score: 2.0
  }
  gene {
    gene_id: 2
    score: 2.0
  }
  gene {
    gene_id: 55
    score: 2.0
  }
  gene {
    gene_id: 1
    score: 3.0
  }
  gene {
    gene_id: 65
    score: 2.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 0.8f, 2}}
},


  // three genes in common same order, all same sign
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 1
    score: 2.0
  }
  gene {
    gene_id: 2
    score: 2.0
  }
  gene {
    gene_id: 3
    score: 2.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 1.0f, 3}}
},

  // three genes in common same order, first pair different signs
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 1
    score: -2.0
  }
  gene {
    gene_id: 2
    score: 2.0
  }
  gene {
    gene_id: 3
    score: 2.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", -4.0f/14.0f, 3}}
},

  // three genes in common same order, second pair different signs
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 1
    score: 2.0
  }
  gene {
    gene_id: 2
    score: -2.0
  }
  gene {
    gene_id: 3
    score: 2.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 6.0f/14.0f, 3}}
},

  // three genes in common same order, third pair different signs
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 1
    score: 2.0
  }
  gene {
    gene_id: 2
    score: 2.0
  }
  gene {
    gene_id: 3
    score: -2.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 12.0f/14.0f, 3}}
},

  // three genes in common same order, third pair different signs. Extra genes
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 21
    score: 1
  }
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 31
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 41
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
  gene {
    gene_id: 51
    score: 1
  }
  gene {
    gene_id: 61
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 91
    score: 2.0
  }
  gene {
    gene_id: 1
    score: 2.0
  }
  gene {
    gene_id: 81
    score: 2.0
  }
  gene {
    gene_id: 2
    score: 2.0
  }
  gene {
    gene_id: 82
    score: 2.0
  }
  gene {
    gene_id: 3
    score: -2.0
  }
  gene {
    gene_id: 83
    score: 2.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 12.0f/14.0f, 3}}
},

  // three genes in common first two out of order
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 2
    score: 2.0
  }
  gene {
    gene_id: 1
    score: 2.0
  }
  gene {
    gene_id: 3
    score: -2.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 11.0f/14.0f, 3}}
},

  // three genes in common, swap 1 and 3
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 3
    score: 2.0
  }
  gene {
    gene_id: 2
    score: 2.0
  }
  gene {
    gene_id: 1
    score: -2.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 4.0f/14.0f, 3}}
},

  // three genes in common, swap 2 and 3
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 1
    score: 2.0
  }
  gene {
    gene_id: 3
    score: 2.0
  }
  gene {
    gene_id: 2
    score: -2.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 9.0f/14.0f, 3}}
},

  // Same as previous, extra genes added
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 101
    score: 1
  }
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 102
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 103
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
  gene {
    gene_id: 104
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 201
    score: 1
  }
  gene {
    gene_id: 1
    score: 2.0
  }
  gene {
    gene_id: 202
    score: 1
  }
  gene {
    gene_id: 3
    score: 2.0
  }
  gene {
    gene_id: 203
    score: 1
  }
  gene {
    gene_id: 2
    score: -2.0
  }
  gene {
    gene_id: 204
    score: 1
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 9.0f/14.0f, 3}}
},

  // Same as previous, adjust haystack_maxrank
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 101
    score: 1
  }
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 102
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 103
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
  gene {
    gene_id: 104
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 201
    score: 1
  }
  gene {
    gene_id: 1
    score: 2.0
  }
  gene {
    gene_id: 202
    score: 1
  }
  gene {
    gene_id: 3
    score: 2.0
  }
  gene {
    gene_id: 203
    score: 1
  }
  gene {
    gene_id: 2
    score: -2.0
  }
  gene {
    gene_id: 204
    score: 1
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 2,  // brings in gene 1
{Neighbour{"t2", 14.0f/14.0f, 1}}
},


  // Same as previous, adjust haystack_maxrank
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 101
    score: 1
  }
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 102
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 103
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
  gene {
    gene_id: 104
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 201
    score: 1
  }
  gene {
    gene_id: 1
    score: 2.0
  }
  gene {
    gene_id: 202
    score: 1
  }
  gene {
    gene_id: 3
    score: 2.0
  }
  gene {
    gene_id: 203
    score: 1
  }
  gene {
    gene_id: 2
    score: -2.0
  }
  gene {
    gene_id: 204
    score: 1
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 4,  // brings in genes 1 and 3
{Neighbour{"t2", 14.0f/14.0f, 2}}
},

  // three genes in common, reverse order
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 3
    score: 2.0
  }
  gene {
    gene_id: 2
    score: 2.0
  }
  gene {
    gene_id: 1
    score: -2.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 4.0f/14.0f, 3}}
},

  // three genes in common, reverse same sign
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 3
    score: 2.0
  }
  gene {
    gene_id: 2
    score: 2.0
  }
  gene {
    gene_id: 1
    score: 2.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", 10.0f/14.0f, 3}}
},

  // three genes in common, reverse order, opposite sign
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 3
    score: -2.0
  }
  gene {
    gene_id: 2
    score: -2.0
  }
  gene {
    gene_id: 1
    score: -2.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", -10.0f/14.0f, 3}}
},

  // three genes in common, reverse order, opposite sign, added genes
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 101
    score: 1
  }
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 102
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 103
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
  gene {
    gene_id: 104
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 201
    score: 1
  }
  gene {
    gene_id: 3
    score: -2.0
  }
  gene {
    gene_id: 202
    score: 1
  }
  gene {
    gene_id: 2
    score: -2.0
  }
  gene {
    gene_id: 203
    score: 1
  }
  gene {
    gene_id: 1
    score: -2.0
  }
  gene {
    gene_id: 204
    score: 1
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", -10.0f/14.0f, 3}}
},


  // three genes in common, same order, opposite sign
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 1
  }
  gene {
    gene_id: 2
    score: 1
  }
  gene {
    gene_id: 3
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 1
    score: -2.0
  }
  gene {
    gene_id: 2
    score: -2.0
  }
  gene {
    gene_id: 3
    score: -2.0
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
10, 10, 10, 10,
{Neighbour{"t2", -14.0f/14.0f, 3}}
},

  // example from Rick Higgs
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: 6.1
  }
  gene {
    gene_id: 2
    score: -4.2
  }
  gene {
    gene_id: 3
    score: 2.5
  }
  gene {
    gene_id: 4
    score: -2.1
  }
  gene {
    gene_id: 5
    score: 1.9
  }
  gene {
    gene_id: 6
    score: 1.8
  }
  gene {
    gene_id: 7
    score: -1.7
  }
  gene {
    gene_id: 8
    score: -1.6
  }
  gene {
    gene_id: 9
    score: 1.4
  }
  gene {
    gene_id: 10
    score: -1.2
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 2
    score: -5.9
  }
  gene {
    gene_id: 1
    score: 5.5
  }
  gene {
    gene_id: 4
    score: -2.8
  }
  gene {
    gene_id: 3
    score: 2.5
  }
  gene {
    gene_id: 5
    score: 2
  }
  gene {
    gene_id: 7
    score: -1.9
  }
  gene {
    gene_id: 6
    score: 1.8
  }
  gene {
    gene_id: 10
    score: -1.5
  }
  gene {
    gene_id: 9
    score: -1.4
  }
  gene {
    gene_id: 8
    score: -1.3
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
20, 20, 20, 20,
{Neighbour{"t2", 370.0f/385.0f, 10}}
},


  // example from Rick Higgs, with extra genes inserted
  Cmp{
R"pb(
  name: "t1",
  gene {
    gene_id: 100
    score: 2.5
  }
  gene {
    gene_id: 101
    score: 2.5
  }
  gene {
    gene_id: 102
    score: 2.5
  }
  gene {
    gene_id: 103
    score: 2.5
  }
  gene {
    gene_id: 1
    score: 6.1
  }
  gene {
    gene_id: 104
    score: 2.5
  }
  gene {
    gene_id: 105
    score: 2.5
  }
  gene {
    gene_id: 106
    score: 2.5
  }
  gene {
    gene_id: 2
    score: -4.2
  }
  gene {
    gene_id: 107
    score: 2.5
  }
  gene {
    gene_id: 3
    score: 2.5
  }
  gene {
    gene_id: 108
    score: 2.5
  }
  gene {
    gene_id: 109
    score: 2.5
  }
  gene {
    gene_id: 110
    score: 2.5
  }
  gene {
    gene_id: 4
    score: -2.1
  }
  gene {
    gene_id: 111
    score: 2.5
  }
  gene {
    gene_id: 112
    score: 2.5
  }
  gene {
    gene_id: 113
    score: 2.5
  }
  gene {
    gene_id: 114
    score: 2.5
  }
  gene {
    gene_id: 5
    score: 1.9
  }
  gene {
    gene_id: 115
    score: 2.5
  }
  gene {
    gene_id: 6
    score: 1.8
  }
  gene {
    gene_id: 7
    score: -1.7
  }
  gene {
    gene_id: 116
    score: 2.5
  }
  gene {
    gene_id: 117
    score: 2.5
  }
  gene {
    gene_id: 8
    score: -1.6
  }
  gene {
    gene_id: 118
    score: 2.5
  }
  gene {
    gene_id: 9
    score: 1.4
  }
  gene {
    gene_id: 119
    score: 2.5
  }
  gene {
    gene_id: 120
    score: 2.5
  }
  gene {
    gene_id: 10
    score: -1.2
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 200
    score: -1.2
  }
  gene {
    gene_id: 201
    score: -1.2
  }
  gene {
    gene_id: 202
    score: -1.2
  }
  gene {
    gene_id: 2
    score: -5.9
  }
  gene {
    gene_id: 1
    score: 5.5
  }
  gene {
    gene_id: 4
    score: -2.8
  }
  gene {
    gene_id: 3
    score: 2.5
  }
  gene {
    gene_id: 203
    score: -1.2
  }
  gene {
    gene_id: 204
    score: -1.2
  }
  gene {
    gene_id: 205
    score: -1.2
  }
  gene {
    gene_id: 5
    score: 2
  }
  gene {
    gene_id: 206
    score: -1.2
  }
  gene {
    gene_id: 207
    score: -1.2
  }
  gene {
    gene_id: 7
    score: -1.9
  }
  gene {
    gene_id: 6
    score: 1.8
  }
  gene {
    gene_id: 10
    score: -1.5
  }
  gene {
    gene_id: 9
    score: -1.4
  }
  gene {
    gene_id: 208
    score: -1.2
  }
  gene {
    gene_id: 8
    score: -1.3
  }
)pb",
// max association, number nbrs,  needle_maxrank, haystack_maxrank
40, 40, 40, 40,
{Neighbour{"t2", 370.0f/385.0f, 10}}
}
));


struct CmpGenes {
  std::string string_proto1;
  std::string string_proto2;
  uint32_t max_association;
  int number_neighbours;
  std::vector<Neighbour> expected;
};

class TestGeneExpressionThreshold : public testing::TestWithParam<std::tuple<int, CmpGenes>> {
  protected:
    gene_expression::Profile _proto1;
    gene_expression::Profile _proto2;

    Needle _n1;
    Needle _n2;
};

TEST_P(TestGeneExpressionThreshold, Test1) {
  auto p = GetParam();
  int large_gene_threshold = std::get<0>(p);
  const CmpGenes& params = std::get<1>(p);

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.string_proto1, &_proto1));
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.string_proto2, &_proto2));

  needle::set_large_gene_id_threshold(large_gene_threshold);

  constexpr int kMaxRank = 1000;
  ASSERT_TRUE(_n1.Build(_proto1, kMaxRank));
  ASSERT_TRUE(_n2.Build(_proto2, kMaxRank));

  // If no threshold, all genes should be in the array.
  if (large_gene_threshold == 0) {
    EXPECT_EQ(_n1.GeneIdsInHash(), 0);
    EXPECT_EQ(_n2.GeneIdsInHash(), 0);
  } else {
    EXPECT_GT(_n1.GeneIdsInHash(), 0);
    EXPECT_GT(_n2.GeneIdsInHash(), 0);
  }

  _n1.SetMaxPossibleAssociation(params.max_association);
  needle::Needle::_number_neighbours = params.number_neighbours;

  _n1.Compare(_proto2, kMaxRank, kMaxRank);

  EXPECT_EQ(_n1.number_neighbours(), params.expected.size());
  if (_n1.number_neighbours() > 0) {
    for (int i = 0; i < _n1.number_neighbours(); ++i) {
      EXPECT_EQ(*_n1.nbr(i), params.expected[i]);
    }
  }
}
INSTANTIATE_TEST_SUITE_P(TestThreshold, TestGeneExpressionThreshold,
  // Values for the large gene id threshold.
  testing::Combine(testing::Range(0, 10), 
  testing::Values(
CmpGenes{
R"pb(
  name: "t1",
  gene {
    gene_id: 1
    score: -8.2
  }
  gene {
    gene_id: 2
    score: 8
  }
  gene {
    gene_id: 3
    score: -7
  } 
  gene {
    gene_id: 4
    score: -6
  }
  gene {
    gene_id: 5
    score: 5
  }
  gene {
    gene_id: 6
    score: -4
  }
  gene {
    gene_id: 9
    score: 3
  }
  gene {
    gene_id: 10
    score: -2
  }
  gene {
    gene_id: 11
    score: 1
  }
)pb",
R"pb(
  name: "t2",
  gene {
    gene_id: 6
    score: -8.2
  }
  gene {
    gene_id: 1
    score: 8
  }
  gene {
    gene_id: 10
    score: -7
  } 
  gene {
    gene_id: 2
    score: -6
  }
  gene {
    gene_id: 3
    score: 5
  }
  gene {
    gene_id: 9
    score: -4
  }
  gene {
    gene_id: 5
    score: 3
  }
  gene {
    gene_id: 7
    score: -2
  }
  gene {
    gene_id: 4
    score: 1
  }
)pb",
// max association, number nbrs
40, 40,
{Neighbour{"t2", -88.0f/204.0f, 8}}
}
)));

}  // namespace
