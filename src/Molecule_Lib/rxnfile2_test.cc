// Tester for functions defined in rxnfile2.cc
#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/rxn_file.h"

namespace {

using std::cerr;
using std::endl;

class TestReactionStringRep : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    const_IWSubstring _rxn_smiles;
    ReactionStringRep _reaction_string_rep;
};

void
TestReactionStringRep::SetUp()
{
  set_unique_smiles_legacy_atom_ordering(1);
  return;
}

TEST_F(TestReactionStringRep, NoReagents)
{
  _rxn_smiles = ">C>C";
  EXPECT_EQ(_reaction_string_rep.Build(_rxn_smiles, '+'), 0);
}

TEST_F(TestReactionStringRep, NoProducts)
{
  _rxn_smiles = "CC>C>";
  EXPECT_EQ(_reaction_string_rep.Build(_rxn_smiles, '+'), 0);
}

TEST_F(TestReactionStringRep, OneEverything)
{
  _rxn_smiles = "R>A>P";
  EXPECT_EQ(_reaction_string_rep.Build(_rxn_smiles, '+'), 1);
  EXPECT_EQ(_reaction_string_rep.number_reagents(), 1);
  EXPECT_EQ(_reaction_string_rep.number_agents(), 1);
  EXPECT_EQ(_reaction_string_rep.number_products(), 1);
}

TEST_F(TestReactionStringRep, MultipleReagents)
{
  _rxn_smiles = "R1+R2>A>P";
  EXPECT_EQ(_reaction_string_rep.Build(_rxn_smiles, '+'), 1);
  EXPECT_EQ(_reaction_string_rep.number_reagents(), 2);
  EXPECT_EQ(_reaction_string_rep.number_agents(), 1);
  EXPECT_EQ(_reaction_string_rep.number_products(), 1);
}

TEST_F(TestReactionStringRep, MultipleAgents)
{
  _rxn_smiles = "R>A1+A2>P";
  EXPECT_EQ(_reaction_string_rep.Build(_rxn_smiles, '+'), 1);
  EXPECT_EQ(_reaction_string_rep.number_reagents(), 1);
  EXPECT_EQ(_reaction_string_rep.number_agents(), 2);
  EXPECT_EQ(_reaction_string_rep.number_products(), 1);
}

TEST_F(TestReactionStringRep, MultipleProducts)
{
  _rxn_smiles = "R>A>P1+P2";
  EXPECT_EQ(_reaction_string_rep.Build(_rxn_smiles, '+'), 1);
  EXPECT_EQ(_reaction_string_rep.number_reagents(), 1);
  EXPECT_EQ(_reaction_string_rep.number_agents(), 1);
  EXPECT_EQ(_reaction_string_rep.number_products(), 2);
}

TEST_F(TestReactionStringRep, DotsAreFragments)
{
  _rxn_smiles = "R1.R2>A1.A2>P1.P2";
  EXPECT_EQ(_reaction_string_rep.Build(_rxn_smiles, '+'), 1);
  EXPECT_EQ(_reaction_string_rep.number_reagents(), 1);
  EXPECT_EQ(_reaction_string_rep.number_agents(), 1);
  EXPECT_EQ(_reaction_string_rep.number_products(), 1);
}

TEST(TestChemAxonFragmentData, EmptyString)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build(""), 1);
}

TEST(TestChemAxonFragmentData, WrongPrefixUnrecognizedAndIgnored)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("F:3.2"), 1);
}

TEST(TestChemAxonFragmentData, EmptySpecification)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("f:"), 0);
}

TEST(TestChemAxonFragmentData, EmptyGroup1)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("f:,"), 0);
}

TEST(TestChemAxonFragmentData, EmptyGroup2)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("f:1,"), 0);
}

TEST(TestChemAxonFragmentData, EmptyGroup3)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("f:,2"), 0);
}

TEST(TestChemAxonFragmentData, GroupTooSmall1)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("f:2"), 0);
}

TEST(TestChemAxonFragmentData, GroupTooSmall2)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("f:1.2.3,4"), 0);
}

TEST(TestChemAxonFragmentData, OneGroup)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("f:1.2.3"), 1);
  const auto& groups = cxfd.connected_group();
  EXPECT_EQ(groups.size(), 1);
  const auto& group0 = *groups[0];
  EXPECT_EQ(group0.size(), 3);
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(i + 1, group0[i]);
  }
}

TEST(TestChemAxonFragmentData, DuplicateFragments)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("f:1.2.3,1.2"), 0);
}

TEST(TestChemAxonFragmentData, InvalidAtom)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("f:1.2.3.Y.2"), 0);
}


TEST(TestChemAxonFragmentData, TwoGroup)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("f:1.2.3,0.4"), 1);
  const auto& groups = cxfd.connected_group();
  EXPECT_EQ(groups.size(), 2);
  const auto& group0 = *groups[0];
  EXPECT_EQ(group0.size(), 3);
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(i + 1, group0[i]);
  }
  const auto & group1 = *groups[1];
  EXPECT_EQ(group1.size(), 2);
  EXPECT_EQ(group1[0], 0);
  EXPECT_EQ(group1[1], 4);
}

// [S:1]([O-:5])([O-:4])(=[O:3])=[O:2].[NH4+:6].[NH4+]>O>[S:1](=[O:3])(=[O:2])([OH:5])[O-:4].[NH4+:6].[S:1]([O-:5])([O-:4])(=[O:3])=[O:2].[NH4+:6].[NH4+:6] |f:0.1.2,4.5,6.7.8|    US03930988              1976            

TEST(Tesof_grouptChemAxonFragmentData, BuildNewSmilesInvalidFragment)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("f:0.1.2,4.5,6.7.8"), 1);
  EXPECT_EQ(cxfd.number_connected_groups(), 3);
  const const_IWSubstring rxn_smiles = "A.B.C.D.E.F>>Q";
  IWString new_rxnsmiles;
  EXPECT_EQ(cxfd.ChangeToPlusForm(rxn_smiles, '+', new_rxnsmiles), 0);
}

TEST(TestChemAxonFragmentData, BuildNewSmiles)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("f:0.1.2,4.5,6.7.8"), 1);
  EXPECT_EQ(cxfd.number_connected_groups(), 3);
  IWString new_rxnsmiles;
  const const_IWSubstring rxn_smiles = "[S:1]([O-:5])([O-:4])(=[O:3])=[O:2].[NH4+:6].[NH4+]>O>[S:1](=[O:3])(=[O:2])([OH:5])[O-:4].[NH4+:6].[S:1]([O-:5])([O-:4])(=[O:3])=[O:2].[NH4+:6].[NH4+:6]";
  EXPECT_TRUE(cxfd.ChangeToPlusForm(rxn_smiles, '+', new_rxnsmiles));
  cerr << "Old is " << rxn_smiles << endl;
  cerr << "New is " << new_rxnsmiles << endl;
  EXPECT_EQ(new_rxnsmiles, "[S:1]([O-:5])([O-:4])(=[O:3])=[O:2].[NH4+:6].[NH4+]>O>[S:1](=[O:3])(=[O:2])([OH:5])[O-:4].[NH4+:6]+[S:1]([O-:5])([O-:4])(=[O:3])=[O:2].[NH4+:6].[NH4+:6]");
}

TEST(TestChemAxonFragmentData, BuildNewSmilesRearange)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("f:0.2.4,1.3.5"), 1);
  EXPECT_EQ(cxfd.number_connected_groups(), 2);
  const const_IWSubstring rxn_smiles = "A.B.C.D.E.F>>Q";
  IWString new_rxnsmiles;
  EXPECT_TRUE(cxfd.ChangeToPlusForm(rxn_smiles, '+', new_rxnsmiles));
  cerr << "Old is " << rxn_smiles << endl;
  cerr << "New is " << new_rxnsmiles << endl;
  EXPECT_EQ(new_rxnsmiles, "A.C.E+B.D.F>>Q");
}

TEST(TestChemAxonFragmentData, FragsOnRhsOfEmptyAgent)
{
  ChemAxonFragmentData cxfd;
  EXPECT_EQ(cxfd.Build("f:2.3,1.4"), 1);
  EXPECT_EQ(cxfd.number_connected_groups(), 2);
  const const_IWSubstring rxn_smiles = "A>>B.C.D.E";
  IWString new_rxnsmiles;
  EXPECT_TRUE(cxfd.ChangeToPlusForm(rxn_smiles, '+', new_rxnsmiles));
  cerr << "Old is " << rxn_smiles << endl;
  cerr << "New is " << new_rxnsmiles << endl;
  EXPECT_EQ(new_rxnsmiles, "A>>B.E+C.D");
}

class TestRxnFile : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    RXN_File _rxn;
};

void
TestRxnFile::SetUp()
{
  _rxn.set_do_automatic_atom_mapping(0);

  return;
}

TEST_F(TestRxnFile, TestRemoveUnchangingFragments1)
{
  const_IWSubstring buffer = "[IH:1].[NH:15]([NH2:16])[C:4]1=[N:10][CH2:9][CH2:8][CH2:7][CH2:6][NH:5]1>>[IH:1].CS[C:4]1=[N:10][CH2:9][CH2:8][CH2:7][CH2:6][NH:5]1.[NH2:15][NH2:16]";

  EXPECT_TRUE(_rxn.build_from_reaction_smiles(buffer, 1));
  EXPECT_EQ(_rxn.number_reagents(), 1);
  const int removed = _rxn.remove_unchanging_fragments();
  EXPECT_EQ(removed, 1);
}

TEST_F(TestRxnFile, TestCommonFragmentIsAComponent)
{
  // Common frag is a fragment on one side, a full component on the other. Not handled here.
  const_IWSubstring buffer = "[CH3:16][S:17](=[O:18])(=[O:19])[OH:31].[Cl:1][C:2]1=[CH:3][CH:4]=[C:5]2[N:6]([CH2:7][CH2:8][CH2:9][OH:10])[C:11](=[O:12])[NH:13][C:14]2=[CH:15]1>[CH3:20][CH2:21][N:22]([CH2:23][CH3:24])[CH2:25][CH3:26]+[Cl:28][CH2:29][Cl:30]>[Cl:1][C:2]1=[CH:3][CH:4]=[C:5]2[N:6]([CH2:7][CH2:8][CH2:9][OH:10])[C:11](=[O:12])[NH:13][C:14]2=[CH:15]1+[CH3:16][S:17](=[O:18])(=[O:19])[Cl:27]";   // 04160836__dup_4216
  EXPECT_TRUE(_rxn.build_from_reaction_smiles(buffer, 1));
  const int removed = _rxn.remove_unchanging_fragments();
  EXPECT_EQ(removed, 0);
}

TEST_F(TestRxnFile, TestAtLeastSomeMappedAtomsCommonBtwReagentsAndProducts)
{
//const_IWSubstring buffer = "[CH3:1][C:2]1[N:3]=[CH:4][C:5]2[C:10]([CH:11]=1)=[C:9]([N+:12]([O-:14])=[O:13])[CH:8]=[CH:7][CH:6]=2.[ClH:15]>>[ClH:15].[CH3:1][C:2]1[N:3]=[CH:4][C:5]2[C:10]([CH:11]=1)=[C:9]([N+:12]([O-:14])=[O:13])[CH:8]=[CH:7][CH:6]=2 |f:2.3| US03930837  1976";
  const_IWSubstring buffer = "CC1N=CC2C(C=1)=C([N+]([O-])=O)C=CC=2.[Cl:15][C:16]1[C:25]2[C:20](=[CH:21][CH:22]=[CH:23][CH:24]=2)[CH:19]=[CH:18][N:17]=1>>[ClH:15].[Cl:15][C:16]1[C:25]2[C:20](=[CH:21][CH:22]=[CH:23][CH:24]=2)[CH:19]=[CH:18][N:17]=1 |f:2.3| US03930837  1976";
  EXPECT_TRUE(_rxn.build_from_reaction_smiles(buffer, 0));
  EXPECT_TRUE(_rxn.at_least_some_mapped_atoms_common_btw_reagents_and_products());
}

TEST_F(TestRxnFile, TestRemoveFragmentsNotParticipating1)
{
  const_IWSubstring buffer = "[Br:1][CH2:2][CH2:3][OH:4].[CH2:5]([S:7](Cl)(=[O:9])=[O:8])[CH3:6].CCOCC>C(N(CC)CC)C>[CH2:5]([S:7]([O:4][CH2:3][CH2:2][Br:1])(=[O:9])=[O:8])[CH3:6]";
  EXPECT_TRUE(_rxn.build_from_reaction_smiles(buffer, 1));
  EXPECT_EQ(_rxn.number_reagents(), 1);
  EXPECT_EQ(_rxn.reagent(0).natoms(), 15);
  EXPECT_EQ(_rxn.remove_fragments_not_participating(), 1);
  EXPECT_EQ(_rxn.reagent(0).number_fragments(), 2);
  EXPECT_EQ(_rxn.reagent(0).natoms(), 10);
}

TEST_F(TestRxnFile, TestRemoveFragmentsNotParticipating2)
{
  const_IWSubstring buffer = "CC1N=CC2C(C=1)=C([N+]([O-])=O)C=CC=2.[Cl:15][C:16]1[CH:25]=[CH:24][C:23]([N+:26]([O-:28])=[O:27])=[C:22]2[C:17]=1[CH:18]=[CH:19][N:20]=[CH:21]2.Cl.CC1N=CC2C(C=1)=C([N+]([O-])=O)C=CC=2.[IH:44]>>[IH:44].[Cl:15][C:16]1[CH:25]=[CH:24][C:23]([N+:26]([O-:28])=[O:27])=[C:22]2[C:17]=1[CH:18]=[CH:19][N:20]=[CH:21]2 |f:2.3,5.6| US03930837  1976";
  EXPECT_TRUE(_rxn.build_from_reaction_smiles(buffer, 0));
  EXPECT_EQ(_rxn.number_reagents(), 2);
  EXPECT_EQ(_rxn.reagent(0).number_fragments(), 3);
  EXPECT_EQ(_rxn.remove_fragments_not_participating(), 3);
  EXPECT_EQ(_rxn.reagent(0).number_fragments(), 2);
  set_include_atom_map_with_smiles(1);
  _rxn.reagent(0).invalidate_smiles();
  EXPECT_EQ(_rxn.reagent(0).unique_smiles(), "I.Clc1c2c(c([N+](=O)[O-])cc1)c[n]cc2");
}
}  // namespace
