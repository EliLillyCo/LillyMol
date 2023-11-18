// Tests for geometric constraints

#include "geometric_constraints.h"
#include "substructure.h"
#include "Molecule_Lib/substructure.pb.h"

//#include "googlemock/include/gmock/gmock.h"
//#include "googletest/include/gtest/gtest.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace {
using geometric_constraints::AllowedRange;
using geometric_constraints::DistanceConstraint;
using geometric_constraints::BondAngleConstraint;
using geometric_constraints::TorsionAngleConstraint;

TEST(Test_Distance, TestAtomNumbers) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC"));
  m.setxyz(0, 0.0, 0.0, 0.0);
  m.setxyz(1, 1.0, 0.0, 0.0);

  DistanceConstraint constraint;
  ASSERT_FALSE(constraint.Active());
  ASSERT_TRUE(constraint.IsValid());
  ASSERT_TRUE(constraint.SetAtoms(0, 1));
  ASSERT_FALSE(constraint.Active());
  ASSERT_FALSE(constraint.IsValid());
  EXPECT_FALSE(constraint.Matches(m));
  ASSERT_TRUE(constraint.SetRange(0.9, 1.0));  // Somewhat dangerous float comparison.
  ASSERT_TRUE(constraint.Active());
  ASSERT_TRUE(constraint.IsValid());
  EXPECT_TRUE(constraint.Matches(m));
  m.setxyz(1, 1.1, 0.0, 0.0);
  EXPECT_FALSE(constraint.Matches(m));
}

TEST(Test_Distance, TestEmbedding) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCC"));
  m.setxyz(0, 0.0, 0.0, 0.0);
  m.setxyz(1, 8.0, 0.0, 0.0);
  m.setxyz(2, 9.0, 0.0, 0.0);

  Set_of_Atoms s;
  s << 2;
  s << 1;

  DistanceConstraint constraint;
  ASSERT_TRUE(constraint.SetAtoms(0, 1));
  ASSERT_TRUE(constraint.SetRange(0.0, 1.1));
  EXPECT_TRUE(constraint.Matches(m, s));
  constraint.SetRange(0.5, 0.9);
  EXPECT_FALSE(constraint.Matches(m, s));
  constraint.SetAtoms(1, 0);
  EXPECT_FALSE(constraint.Matches(m, s));
}

TEST(Test_Bond_Angle, TestAtomNumbersDegrees) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCC"));
  m.setxyz(0, 0.0, 0.0, 0.0);
  m.setxyz(1, 1.0, 0.0, 0.0);
  m.setxyz(2, 2.0, 0.0, 0.0);

  BondAngleConstraint constraint;
  constraint.SetAtoms(0, 1, 2);
  constraint.set_range_in_degrees(180.0 - 0.0001, 180.0 + 0.0001);
  EXPECT_TRUE(constraint.Matches(m));
  constraint.SetAtoms(2, 1, 0);
  EXPECT_TRUE(constraint.Matches(m));

  m.setxyz(2, 1.0 - 0.7071067811865475, 0.7071067811865475, 0.0);
  EXPECT_FALSE(constraint.Matches(m));
  constraint.set_range_in_degrees(40.0, 50.0);
  EXPECT_TRUE(constraint.Matches(m));
  constraint.SetAtoms(2, 1, 0);
  EXPECT_TRUE(constraint.Matches(m));

  m.setxyz(2, 1.0 + 0.7071067811865475, 0.7071067811865475, 0.0);
  EXPECT_FALSE(constraint.Matches(m));
  constraint.set_range_in_degrees(135.0 - 0.0001, 135.0 + 0.0001);
  EXPECT_TRUE(constraint.Matches(m));
}

TEST(Test_Bond_Angle, TestEmbeddingRadians) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCCCCCCC"));
  m.setxyz(2, 0.0, 0.0, 0.0);
  m.setxyz(3, 0.0, 1.0, 0.0);
  m.setxyz(7, 1.0, 1.0, 0.0);

  BondAngleConstraint constraint;
  constraint.SetAtoms(0, 1, 2);
  constraint.set_range_in_radians(M_PI / 2 - 0.0001, M_PI + 0.0001);

  Set_of_Atoms s;
  s << 2 << 3 << 7;
  EXPECT_TRUE(constraint.Matches(m));
  m.setxyz(7, 1.0, -1.0, 0.0);
  EXPECT_TRUE(constraint.Matches(m));
}

TEST(Test_Torsion_Angle, TestDegrees) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCCC"));
  m.setxyz(0, -0.5, -0.5, 0.0);
  m.setxyz(1, 0.0, 0.0, 0.0);
  m.setxyz(2, 1.0, 0.0, 0.0);
  m.setxyz(3, 1.5, 0.5, 0.1);

  TorsionAngleConstraint constraint;
  constraint.SetAtoms(0, 1, 2, 3);
  constraint.set_range_in_degrees(168.69 - 0.1, 168.69 + 0.1);
  EXPECT_TRUE(constraint.Matches(m));
  m.setxyz(3, 1.5, 0.5, -0.1);
  EXPECT_TRUE(constraint.Matches(m));
  constraint.set_range_in_degrees(168.69 - 0.01, 168.69 + 0.01);
  EXPECT_TRUE(constraint.Matches(m));
  m.setxyz(3, 0.5, 0.5, -0.1);
  EXPECT_TRUE(constraint.Matches(m));
}

TEST(Test_Range, FromProtoFails) {
  std::string string_proto = R"pb(
    min: 1.0
    max: 0.1
  )pb";
  GeometricConstraints::Range proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto, &proto));
  AllowedRange constraint;
  EXPECT_FALSE(constraint.BuildFromProto(proto));
}

TEST(Test_Torsion_Angle, FromProtoDuplicatedAtoms) {
  std::string string_proto = R"pb(
      range {
         min: -1.5709
         max: -1.5707
      }
      a1: 3
      a2: 6
      a3: 3
      a4: 8
  )pb";

  GeometricConstraints::TorsionAngle proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto, &proto));
  TorsionAngleConstraint constraint;
  ASSERT_FALSE(constraint.BuildFromProto(proto));
}

TEST(Test_Distance, FromProto) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCCCCCC"));
  m.setxyz(3, 0.0, 0.0, 0.0);
  m.setxyz(6, 1.0, 1.0, 1.0);
  std::string string_proto = R"pb(
      range {
         min: 1.732
         max: 1.734
      }
      a1: 3
      a2: 6
  )pb";

  GeometricConstraints::Distance proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto, &proto));
  DistanceConstraint constraint;
  ASSERT_TRUE(constraint.BuildFromProto(proto));
  EXPECT_TRUE(constraint.Matches(m));
}

TEST(Test_Angle, FromProto) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCCCCCCC"));
  m.setxyz(5, 0.0, 0.0, 0.0);
  m.setxyz(6, 1.0, 0.0, 0.0);
  m.setxyz(7, -10.0, -0.1, 0.0);
  std::string string_proto = R"pb(
      range {
         min: 0.1
         max: 1.5
      }
      a1: 5
      a2: 6
      a3: 7
  )pb";

  GeometricConstraints::BondAngle proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto, &proto));
  BondAngleConstraint constraint;
  ASSERT_TRUE(constraint.BuildFromProto(proto));
  EXPECT_TRUE(constraint.Matches(m));
}

TEST(Test_Torsion_Angle, FromProto) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCCCCCCCC"));
  m.setxyz(5, 0.0, -1.0, 0.0);
  m.setxyz(6, 0.0, 0.0, 0.0);
  m.setxyz(7, 1.0, 0.0, 0.0);
  m.setxyz(8, 1.0, 0.0, 1.0);
  std::string string_proto = R"pb(
      range {
         min: 89.0
         max: 91.9
      }
      a1: 5
      a2: 6
      a3: 7
      a4: 8
  )pb";

  GeometricConstraints::TorsionAngle proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto, &proto));
  TorsionAngleConstraint constraint;
  ASSERT_TRUE(constraint.BuildFromProto(proto));
  EXPECT_TRUE(constraint.Matches(m));
}

TEST(TestSet, TestAllFeatures) {
  std::string string_proto = R"pb(
    distance {
      range {
        min: 0.9
        max: 1.1
      }
      a1: 0
      a2: 1
    }
    distance {
      range {
        min: 0.95
        max: 1.05
      }
      a1: 1
      a2: 2
    }
    bond_angle {
      range {
        min: 134.90
        max: 135.10
      }
      a1: 0
      a2: 1
      a3: 2
    }
    torsion_angle {
      range {
        min: 56.3
        max: 56.7
      }
      a1: 0
      a2: 1
      a3: 2
      a4: 3
    }
    number_to_match: 3
  )pb";

  GeometricConstraints::SetOfConstraints proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto, &proto));
  ASSERT_EQ(proto.distance().size(), 2);
  ASSERT_EQ(proto.bond_angle().size(), 1);
  ASSERT_EQ(proto.torsion_angle().size(), 1);
  geometric_constraints::SetOfGeometricConstraints constraints;
  ASSERT_TRUE(constraints.BuildFromProto(proto));

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCCC"));
  m.setxyz(0, -0.7071067811865475, 0.7071067811865475, 0.0);
  m.setxyz(1, 0.0, 0.0, 0.0);
  m.setxyz(2, 1.0, 0.0, 0.0);
  m.setxyz(3, 2.0, 0.2, 0.3);
  EXPECT_TRUE(constraints.Matches(m));
}

TEST(TestDistance, FromSmiles) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C{{0,0,0}}C{{1.5,0,0}}"));
  std::string string_proto = R"pb(
    distance {
      range {
        min: 1.4
        max: 1.6
      }
      a1: 0
      a2: 1
    }
    number_to_match: 1
  )pb";
  GeometricConstraints::SetOfConstraints proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto, &proto));
  ASSERT_EQ(proto.distance().size(), 1);

  geometric_constraints::SetOfGeometricConstraints constraints;
  ASSERT_TRUE(constraints.BuildFromProto(proto));
  EXPECT_GT(constraints.Matches(m), 0);

  Set_of_Atoms embedding;
  embedding << 1 << 0;
  EXPECT_TRUE(constraints.Matches(m, embedding));
}

struct SmilesConstraint {
  IWString smiles;
  std::string proto;
  int expected_result;
};

class TestPositionalConstraintsP: public testing::TestWithParam<SmilesConstraint> {
  protected:
    SubstructureSearch::SubstructureQuery _proto;
    Substructure_Query _query;
    Molecule _mol;
    Substructure_Results _sresults;
};

TEST_P(TestPositionalConstraintsP, Test) {
  const auto& params = GetParam();
  std::cerr << params.proto << '\n';
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  EXPECT_EQ(_query.substructure_search(_mol, _sresults), params.expected_result);
}
INSTANTIATE_TEST_SUITE_P(TestPositionalConstraintsP, TestPositionalConstraintsP, testing::Values(
  SmilesConstraint{"C", 
R"pb(
query {
  smarts: "C"
}
)pb", 1},

  SmilesConstraint{"C{{1,1,1}}", 
R"pb(
query {
  smarts: "C"
}
)pb", 1},

  SmilesConstraint{"C{{1,1,1}}", 
R"pb(
query {
  smarts: "C"
  geometric_constraints {
    positive_spatial_constraint {
      atom: 0
      position {
        x: 0.95
        y: 1
        z: 1.05
      }
      range {
        min: 0.05
        max: 0.2
      }
      match_type: AVERAGE
    }
  }
}
)pb", 1},

  SmilesConstraint{"C{{1,1,1}}", 
R"pb(
query {
  smarts: "C"
  geometric_constraints {
    positive_spatial_constraint {
      atom: 0
      position {
        x: 0.95
        y: 1
        z: 1.05
      }
      range {
        min: 0.05
        max: 0.2
      }
      match_type: ANY
    }
  }
}
)pb", 1},

  SmilesConstraint{"C{{1,1,1}}", 
R"pb(
query {
  smarts: "C"
  geometric_constraints {
    positive_spatial_constraint {
      atom: 0
      position {
        x: 0.95
        y: 1
        z: 1.05
      }
      range {
        min: 0.05
        max: 0.2
      }
      match_type: ALL
    }
  }
}
)pb", 1},

  SmilesConstraint{"C{{1,1,1}}", 
R"pb(
query {
  smarts: "C"
  geometric_constraints {
    positive_spatial_constraint {
      atom: 0
      position {
        x: -0.95
        y: -1
        z: -1.05
      }
      range {
        min: 0.05
        max: 0.2
      }
      match_type: AVERAGE
    }
  }
}
)pb", 0},

  SmilesConstraint{"C{{1,1,1}}", 
R"pb(
query {
  smarts: "C"
  geometric_constraints {
    positive_spatial_constraint {
      atom: 0
      position {
        x: -0.95
        y: -1
        z: -1.05
      }
      range {
        min: 0.05
        max: 0.2
      }
      match_type: ANY
    }
  }
}
)pb", 0},

  SmilesConstraint{"C{{1,1,1}}", 
R"pb(
query {
  smarts: "C"
  geometric_constraints {
    positive_spatial_constraint {
      atom: 0
      position {
        x: -0.95
        y: -1
        z: -1.05
      }
      range {
        min: 0.05
        max: 0.2
      }
      match_type: ALL
    }
  }
}
)pb", 0},

  SmilesConstraint{"C{{1,1,1}}C{{2,2,2}}C{{3,3,3}}", 
R"pb(
query {
  smarts: "CCC"
  geometric_constraints {
    positive_spatial_constraint {
      atom: [0, 1, 2]
      position {
        x: 2.
        y: 2.
        z: 2.
      }
      range {
        max: 0.1
      }
      match_type: AVERAGE
    }
  }
}
)pb", 2},

  SmilesConstraint{"C{{1,1,1}}C{{2,2,2}}C{{3,3,3}}", 
R"pb(
query {
  smarts: "CCC"
  geometric_constraints {
    positive_spatial_constraint {
      atom: [0, 1, 2]
      position {
        x: 2.
        y: 2.
        z: 2.
      }
      range {
        max: 0.1
      }
      match_type: ANY
    }
  }
}
)pb", 2},

  SmilesConstraint{"C{{1,1,1}}C{{2,2,2}}C{{3,3,3}}", 
R"pb(
query {
  smarts: "CCC"
  geometric_constraints {
    positive_spatial_constraint {
      atom: [0, 1, 2]
      position {
        x: 2.
        y: 2.
        z: 2.
      }
      range {
        max: 0.1
      }
      match_type: ALL
    }
  }
}
)pb", 0},

  SmilesConstraint{"C{{1,1,1}}C{{0.99,0.99,0.99}}C{{1.02,1.02,1.02}}", 
R"pb(
query {
  smarts: "CCC"
  geometric_constraints {
    positive_spatial_constraint {
      atom: [0, 1, 2]
      position {
        x: 1.
        y: 1.
        z: 1.
      }
      range {
        max: 0.1
      }
      match_type: ALL
    }
  }
}
)pb", 2},

  SmilesConstraint{"C{{1,1,1}}C{{0.99,0.99,0.99}}C{{1.02,1.02,1.02}}", 
R"pb(
query {
  smarts: "CCC"
  geometric_constraints {
    negative_spatial_constraint {
      atom: [0, 1, 2]
      position {
        x: 1.
        y: 1.
        z: 1.
      }
      range {
        max: 0.1
      }
      match_type: ALL
    }
  }
}
)pb", 0},

  SmilesConstraint{"C{{1,1,1}}C{{0.99,0.99,0.99}}C{{2.02,1.02,1.02}}", 
R"pb(
query {
  smarts: "CCC"
  geometric_constraints {
    negative_spatial_constraint {
      atom: [0, 1, 2]
      position {
        x: 1.
        y: 1.
        z: 1.
      }
      range {
        max: 0.1
      }
      match_type: ALL
    }
  }
}
)pb", 2}

));

}  // namespace
