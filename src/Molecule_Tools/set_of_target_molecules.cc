#include <iostream>
#include <memory>

#include "Utilities/GFP_Tools/gfp_standard.h"

#include "set_of_target_molecules.h"

using std::cerr;

template <typename F>
Set_of_Target_Molecules<F>::Set_of_Target_Molecules() : _lfp(&_lfpd)
{
  _fp = nullptr;

  _atom_typing_specification.build("UST:ARY");

  _turn_off_if_within = -1.0f;

  _match_already_found = nullptr;

  return;
}

template <typename F>
Set_of_Target_Molecules<F>::~Set_of_Target_Molecules()
{
  if (nullptr != _fp)
    delete [] _fp;

  if (nullptr != _match_already_found)
    delete [] _match_already_found;

  return;
}

template <typename F>
int
Set_of_Target_Molecules<F>::build(const char * fname)
{
  FileType itype = discern_file_type_from_name(fname);

  if (FILE_TYPE_INVALID == itype)
  {
    cerr << "Set_of_Target_Molecules::build:cannot discern type '" << fname << "', assuming smiles\n";
    itype = FILE_TYPE_SMI;
  }

  data_source_and_type<Molecule> input(itype, fname);

  if (! input.good())
  {
    cerr << "Set_of_Target_Molecules::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return build(input);
}


template <typename F>
int
Set_of_Target_Molecules<F>::set_turn_off_molecules_matched_to_within(const float d)
{
  assert (nullptr == _match_already_found);

  const int n = _m.number_elements();

  assert (n > 0);

  if (0 == n)
  {
    cerr << "Set_of_Target_Molecules::set_turn_off_exactly_matched_molecules:empty\n";
    return 0;
  }

  _turn_off_if_within = d;

  _match_already_found = new_int(n);

  return 1;
}

template <typename F>
void
Set_of_Target_Molecules<F>::_compute_fingerprint(Molecule & m,
                                                 F & fp)
{
  _mpr(m, fp.molecular_properties());

  _mk(m, _tmp);
  fp.build_mk(_tmp, _mk.nbits());
  _mk.set_level_2_fingerprint(_tmp);
  fp.build_mk2(_tmp, _mk.nbits());

  const int matoms = m.natoms();

  int * atype = new int[matoms]; std::unique_ptr<int[]> free_atype(atype);

  _atom_typing_specification.assign_atom_types(m, atype);

  LFP::Fixed_Width_Fingerprint_Bits b;
  _lfp.process(m, atype, b);

  fp.build_iwfp(b.begin(), b.nset());    // we assume contiguously stored data

  return;
}

template <typename F>
float
Set_of_Target_Molecules<F>::closest_distance(F & f) const
{
  float highest_similarity = f.tanimoto(_fp[0]);

  const int n = _m.number_elements();

  for (int i = 1; i < n; ++i)
  {
    const float s = f.tanimoto(_fp[i]);

    if (s > highest_similarity)
      highest_similarity = s;
  }

//cerr << "Set_of_Target_Molecules::closest_distance:scanned " << n << " fingerprints, highest_similarity " << highest_similarity << '\n';

  return 1.0f - highest_similarity;
}

template <typename F>
int
Set_of_Target_Molecules<F>::build(data_source_and_type<Molecule> & input)
{
  Molecule * m;

  while (nullptr != (m = input.next_molecule()))
  {
    _m.add(m);
  }

  const int n = _m.number_elements();

  if (0 == n)
  {
    cerr << "Set_of_Target_Molecules::build:no molecules\n";
    return 0;
  }

  _fp = new F[n];

  for (int i = 0; i < n; ++i)
  {
    _compute_fingerprint(*_m[i], _fp[i]);
  }

  return 1;
}

template <typename F>
float
Set_of_Target_Molecules<F>::closest_distance(Molecule & m)
{
  F fp;

  _compute_fingerprint(m, fp);

  return closest_distance(fp);
}

template <typename F>
void
Set_of_Target_Molecules<F>::closest_distance(Molecule & m,
                                             Accumulator<float> & acc)
{
  F fp;
  _compute_fingerprint(m, fp);

  const int n = _m.number_elements();

//cerr << "Set_of_Target_Molecules::closest_distance:n = " << n << '\n';

  if (nullptr == _match_already_found)
  {
    for (int i = 0; i < n; ++i) {
      acc.extra(1.0f - fp.tanimoto(_fp[i]));
    }
  }
  else
  {
    for (int i = 0; i < n; ++i) {
      if (_match_already_found[i]) {
        continue;
      }

      const float d = 1.0f - fp.tanimoto(_fp[i]);
      if (d < _turn_off_if_within) {
        _match_already_found[i] = 1;
      }

      acc.extra(d);
    }
  }
}

template class Set_of_Target_Molecules<GFP_Standard>;
