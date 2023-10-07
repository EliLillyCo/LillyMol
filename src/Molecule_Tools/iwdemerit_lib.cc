#include <iostream>
#include "Molecule_Lib/target.h"

#include "iwdemerit_lib.h"

namespace lilly_medchem_rules {

using std::cerr;
using std::endl;

template <typename T>
float 
Fraction(int numerator, int denominator) {
  return static_cast<T>(numerator) / static_cast<T>(denominator);
}

void
DisplayAtomCutoffOptions(std::ostream & output)
{
  output << " -c hmin=<n>           hard min cutoff - molecules with fewer atoms rejected\n";
  output << " -c smin=<h>           soft min cutoff - molecules with fewer atoms demerited\n";
  output << " -c smax=<h>           soft max cutoff - molecules with more  atoms demerited\n";
  output << " -c hmax=<h>           hard max cutoff - molecules with more  atoms rejected\n";
}

MCDemerit::MCDemerit() {
}

int
MCDemerit::Build(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    const_IWSubstring c;
    for (int i = 0; cl.value('c', c, i); ++i) {
      if ("help" == c) {
        DisplayAtomCutoffOptions(cerr);
        return 0;
      }

      const_IWSubstring directive;
      int dvalue;
      if (! c.split_into_directive_and_value(directive, '=', dvalue)) {
        cerr << "Invalid -c directive '" << c << "'\n";
        return 0;
      }

      if (dvalue < 0) {
        cerr << "INvalid numeric directive, cannot be negative '" << c << "'\n";
        return 0;
      }

      if ("hmin" == directive) {
        _hard_lower_atom_count_cutoff = dvalue;
      } else if ("smin" == directive) {
        _soft_lower_atom_count_cutoff = dvalue;
      } else if ("smax" == directive) {
        _soft_upper_atom_count_cutoff = dvalue;
      } else if ("hmax" == directive) {
        _hard_upper_atom_count_cutoff = dvalue;
      } else if ("mindmrt" == directive) {
        ;
      } else if ("maxdmrt" == directive) {
        ;
      } else {
        cerr << "Unrecognised -c qualifier '" << c << "'\n";
        return 0;
      }
    }

    if (0 == _soft_lower_atom_count_cutoff && 0 == _hard_lower_atom_count_cutoff)
      ;
    else if (_hard_lower_atom_count_cutoff <= _soft_lower_atom_count_cutoff)
      ;
    else {
      cerr << "Invalid lower atom count cutoff values soft:" << _soft_lower_atom_count_cutoff << " hard:" << _hard_lower_atom_count_cutoff << endl;
      return 0;
    }

    if (0 == _soft_upper_atom_count_cutoff && 0 == _hard_upper_atom_count_cutoff)
      ;
    else if (_soft_upper_atom_count_cutoff < _hard_upper_atom_count_cutoff)
      ;
    else {
      cerr << "Invalid upper atom count cutoff values soft:" << _soft_upper_atom_count_cutoff << " hard:" << _hard_upper_atom_count_cutoff << endl;
      return 5;
    }

    if (0 == _soft_lower_atom_count_cutoff && 0 == _soft_upper_atom_count_cutoff)
      ;
    else if (_soft_lower_atom_count_cutoff < _soft_upper_atom_count_cutoff)
      ;
    else {
      cerr << "Soft cutoffs invalid, lower:" << _soft_lower_atom_count_cutoff << " upper:" << _soft_upper_atom_count_cutoff << endl;
      return 0;
    }

//  Should do more checks here...
  }

  return 1;
}

Demerit
MCDemerit::Process(Molecule& m,
                   int demerit_numeric_value_index) {
  Demerit result;
  DoAtomCountDemerits(m, result);
  if (result.rejected()) {
    return result;
  }

  if (_whole_molecule.size()) {
    Molecule_to_Match target(&m);
    RunASetOfQueries(target, _whole_molecule, demerit_numeric_value_index, result);
    if (result.rejected()) {
      return result;
    }
  }

  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
    Molecule_to_Match target(&m);

    RunASetOfQueries(target, _largest_fragment, demerit_numeric_value_index, result);

    if (result.rejected()) {
      return result;
    }
  }

  return result;
}

int
OneIfZero(int d) {
  if (d == 0) {
    return 1;
  }
  return d;
}

void
MCDemerit::DoAtomCountDemerits(Molecule & m, Demerit & demerit) const {
  const int nf = m.number_fragments();

  int matoms;
  if (1 == nf)
    matoms = m.natoms();
  else
  {
    matoms = 0;
    for (int i = 0; i < nf; i++)
    {
      int f = m.atoms_in_fragment(i);
      if (f > matoms)
        matoms = f;
    }
  }

  if (_hard_lower_atom_count_cutoff > 0 && matoms <= _hard_lower_atom_count_cutoff)
  {
    demerit.extra(_lower_atom_count_demerit, "too_few_atoms");

    return;
  }

// Feb 2005. Heuristic for adding demerits beyond the hard atom count cutoff

  if (_hard_upper_atom_count_cutoff > 0 && matoms >= _hard_upper_atom_count_cutoff)
  {
    demerit.extra(rejection_threshold() + 6 * (matoms - _hard_upper_atom_count_cutoff), "too_many_atoms");

//  cerr << hard_upper_atom_count_cutoff << " hard_upper_atom_count_cutoff, matoms = " << matoms << endl;
    return;
  }

  if (matoms > _hard_lower_atom_count_cutoff && matoms < _soft_lower_atom_count_cutoff) {
    float r = Fraction<float>(_soft_lower_atom_count_cutoff - matoms,
                       _soft_lower_atom_count_cutoff - _hard_lower_atom_count_cutoff);
    int d = static_cast<int>(_lower_atom_count_demerit * r);

    demerit.extra(OneIfZero(d), "too_few_atoms");
    return;
  }

  if (matoms > _soft_upper_atom_count_cutoff && matoms < _hard_upper_atom_count_cutoff) {
    float r = Fraction<float>(matoms - _soft_upper_atom_count_cutoff,
                       _hard_upper_atom_count_cutoff - _soft_upper_atom_count_cutoff);
    int d = static_cast<int>(_upper_atom_count_demerit * r);
//  cerr << "too_many_atoms: " << soft_upper_atom_count_cutoff << " soft_upper_atom_count_cutoff, matoms = " << matoms << ", d = " << d << endl;
//  cerr << "upper_atom_count_demerit " << upper_atom_count_demerit << endl;
//  cerr << "too_many_atoms: soft_upper_atom_count_cutoff " << soft_upper_atom_count_cutoff << " hard_upper_atom_count_cutoff " << hard_upper_atom_count_cutoff << " matoms " << matoms << " D = " << d << endl;

    demerit.extra(OneIfZero(d), "too_many_atoms");
    return;
  }

  return;
}

void
MCDemerit::RunASetOfQueries(Molecule_to_Match & target,
                 resizable_array_p<Substructure_Hit_Statistics> & queries,
                 int demerit_numeric_value_index,
                 Demerit & demerit) const {
  for (Substructure_Hit_Statistics * q : queries) {
    int nhits = q->substructure_search(target);
    if (0 == nhits) {
      continue;
    }

    if (_verbose > 1) {
      cerr << nhits << " matches to query " << q->comment() << endl;
    }

    double d;
    (void) q->numeric_value(d, demerit_numeric_value_index);
    
    d = d * nhits;

    int intd = static_cast<int>(d * 1.0001);

    if (intd >= rejection_threshold()) {
      demerit.reject(q->comment());
      if (0 == _keep_going_after_rejection) {
        return;
      }
    } else {
      demerit.extra(intd, q->comment());
    }

    if (demerit.rejected() && 0 == _keep_going_after_rejection) {
      return;
    }
  }

  return;
}

// Divide the queries in `queries` into `_largest_fragment` or
// `_whole_molecule` depending on the query's value of
// only_keep_matches_in_largest_fragment.
// Not sure this is really necessary, if the query only matches
// the largest fragment, there is no need to strip fragments.
// Retained as is for historical reasons.
void
MCDemerit::SeparateDependingOnFragmentMatch(resizable_array_p<Substructure_Hit_Statistics> & queries) {
  for (int i = queries.number_elements() - 1; i >= 0; i--)
  {
    Substructure_Hit_Statistics * q = queries.remove_no_delete(i);

    if (q->only_keep_matches_in_largest_fragment())
      _largest_fragment.add(q);
    else
      _whole_molecule.add(q);
  }

  if (_verbose > 1) {
    cerr << "Separated queries into " << _largest_fragment.number_elements() <<
            " and " << _whole_molecule.number_elements() << " queries\n";
  }

  return;
}

}  // namespace lilly_medchem_rules
