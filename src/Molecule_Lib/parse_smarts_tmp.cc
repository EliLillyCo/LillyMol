#include <stdlib.h>
#include <iostream>

#include "substructure.h"
#include "parse_smarts_tmp.h"

using std::cerr;
using std::endl;

Parse_Smarts_Tmp::Parse_Smarts_Tmp ()
{
  _last_query_atom_created = -1;

  return;
}

Parse_Smarts_Tmp::~Parse_Smarts_Tmp ()
{
  DeletePointers();
  return;
}

int
Parse_Smarts_Tmp::set_natoms (int n)
{
  assert (n > 0);

//cerr << "Parse_Smarts_Tmp:set_natoms: natoms " << n << endl;

  return 1;
}

ThreeDots::ThreeDots(int a1, int a2) : AtomsAndQualifier(a1, a2) {
}

int
ThreeDots::BuildFromProto(const SubstructureSearch::NoMatchedAtomsBetween& proto)
{
  if (! proto.has_a1() || ! proto.has_a2()) {
    cerr << "ThreeDots::BuildFromProto:proto missing atom(s) " << proto.ShortDebugString() << endl;
    return 0;
  }

  _matched_atom_1 = proto.a1();
  _matched_atom_2 = proto.a2();

  if (_matched_atom_1 == _matched_atom_2) {
    cerr << "ThreeDots::BuildFromProto:Invalid matched atoms " << proto.ShortDebugString() << endl;
    return 0;
  }

  if (proto.has_qualifier()) {
    _qualifier = proto.qualifier();
  }

  return 1;
}

template <typename T>
void
DeleteAll(resizable_array<T*>& values) {
  for (T * v : values) {
    delete v;
  }
  values.resize(0);
}

void
Parse_Smarts_Tmp::DeletePointers() {
  DeleteAll(_no_matched_atoms_between);
  DeleteAll(_link_atom);
  DeleteAll(_three_dots);
}

DownTheBondInfo::DownTheBondInfo(int a1) {
  _matched_atom_1 = a1;
}

int
Parse_Smarts_Tmp::DownTheBondSpecificationsComplete() const {
  for (const down_the_bond::DownTheBond* dtb : _down_the_bond) {
    if (! dtb->Complete()) {
      return 0;
    }
  }

  return 1;
}
