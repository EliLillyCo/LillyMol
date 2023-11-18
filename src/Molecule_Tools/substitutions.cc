/*
  Sometimes we want the count the number of substitutions to a substructure
  query match
*/

#include <algorithm>
#include <deque>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <vector>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/histogram/iwhistogram.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/iwminmax.h"
#include "Foundational/iwmisc/minmaxspc.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/mdl_molecule.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

using std::cerr;
using std::endl;

static int add_properties_to_groups = 0;
static const char * prop_header[] = {"natoms", "nrings", "amw", "ro5_ohnh", "ro5_on", "nvrtspsa",  "halogen",  "fcsp3",  "ringsys",  "arring",  "alring", "ringatom",  "mxdst"};
static int prop_number = 13;

static double
novartis_polar_suface_area_nitrogen (Molecule & m,
                                     atom_number_t zatom,
                                     Atom & a,
                                     int is_aromatic)
{
  int ncon = a.ncon();

  int aromatic_bonds = 0;
  int single_bonds = 0;
  int double_bonds = 0;
  int triple_bonds = 0;

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = a[i];

    if (b->is_aromatic())
      aromatic_bonds++;
    else if (b->is_single_bond())
      single_bonds++;
    else if (b->is_double_bond())
      double_bonds++;
    else
      triple_bonds++;
  }

  int hcount = a.implicit_hydrogens();

  formal_charge_t fc = a.formal_charge();

  if (is_aromatic)
  {
    if (0 == fc && 0 == hcount && 2 == ncon && 2 == aromatic_bonds)   // [n](:*):*
      return 12.89;
    if (0 == fc && 0 == hcount && 3 == ncon && 3 == aromatic_bonds)   // [n](:*)(:*):*
      return 4.41;
    if (0 == fc && 0 == hcount && 3 == ncon && 2 == aromatic_bonds && 1 == single_bonds)    // [n](-*)(:*):*
      return 4.93;
    if (0 == fc && 0 == hcount && 3 == ncon && 2 == aromatic_bonds && 1 == double_bonds)    // [n](=*)(:*):*
      return 8.39;
    if (0 == fc && 1 == hcount && 2 == ncon && 2 == aromatic_bonds)   // [nH](:*):*
      return 15.79;
    if (1 == fc && 0 == hcount && 3 == ncon && 3 == aromatic_bonds)   // [n+](:*)(:*):*
      return 4.10;
    if (1 == fc && 0 == hcount && 3 == ncon && 2 == aromatic_bonds && 1 == single_bonds)    // [n+](-*):*):*
      return 3.88;
    if (1 == fc && 1 == hcount && 2 == ncon && 2 == aromatic_bonds)   // [nH+](:*):*
      return 14.14;
  }
  else
  {
    if (0 == fc && 0 == hcount && 3 == ncon && 3 == single_bonds && m.in_ring_of_given_size(zatom, 3))    // case of N1C=C1 maybe handled incorrectly       // [N]1(-*)-*-*-1     do ring queries first
      return 3.01;
    if (0 == fc && 1 == hcount && 2 == ncon && 2 == single_bonds && m.in_ring_of_given_size(zatom, 3))       // [NH]1-*-*-1     do ring queries first
      return 21.94;

    if (0 == fc && 0 == hcount && 3 == ncon && 3 == single_bonds)      // [N](-*)-*
      return 3.24;
    if (0 == fc && 0 == hcount && 2 == ncon && 1 == single_bonds && 1 == double_bonds)   // [N](-*)=*
      return 12.36;
    if (0 == fc && 0 == hcount && 1 == ncon && 1 == triple_bonds)       // [N]#*
      return 23.79;
    if (0 == fc && 0 == hcount && 3 == ncon && 1 == single_bonds && 2 == double_bonds)       // [N](-*)(=*)(=*)
      return 11.68;
    if (0 == fc && 0 == hcount && 2 == ncon && 1 == double_bonds && 1 == triple_bonds)       // [N](=*)#*
      return 13.60;
    if (0 == fc && 1 == hcount && 2 == ncon && 2 == single_bonds)       // [NH](-*)-*
      return 12.03;
    if (0 == fc && 1 == hcount && 1 == ncon && 1 == double_bonds)       // [NH]=*
      return 23.85;
    if (0 == fc && 2 == hcount && 1 == ncon && 1 == single_bonds)       // [NH2]-*
      return 26.02;
    if (1 == fc && 0 == hcount && 4 == ncon && 4 == single_bonds)     // [N+](-*)(-*)(-*)-*
      return 0.0;
    if (1 == fc && 0 == hcount && 3 == ncon && 2 == single_bonds && 1 == double_bonds)     // [N+](-*)(-*)=*
      return 3.01;
    if (1 == fc && 0 == hcount && 2 == ncon && 1 == single_bonds && 1 == triple_bonds)     // [N+](-*)#*
      return 4.36;
    if (1 == fc && 1 == hcount && 3 == ncon && 3 == single_bonds)     // [NH+](-*)(-*)-*
      return 4.44;
    if (1 == fc && 1 == hcount && 2 == ncon && 1 == single_bonds && 1 == double_bonds)     // [NH+](-*)=*
      return 13.97;
    if (1 == fc && 2 == hcount && 2 == ncon && 2 == single_bonds)     // [NH2+](-*)-*
      return 16.61;
    if (1 == fc && 2 == hcount && 1 == ncon && 1 == double_bonds)     // [NH2+]=*
      return 25.59;
    if (1 == fc && 3 == hcount && 1 == single_bonds)     // [NH3+]-*
      return 27.64;
  }

  return 0.0;
}

static double
novartis_polar_suface_area_oxygen (Molecule & m,
                                   atom_number_t zatom,
                                   Atom & a,
                                   int is_aromatic)
{
  int ncon = a.ncon();

  int aromatic_bonds = 0;
  int single_bonds = 0;
  int double_bonds = 0;
  int triple_bonds = 0;

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = a[i];

    if (b->is_aromatic())
      aromatic_bonds++;
    else if (b->is_single_bond())
      single_bonds++;
    else if (b->is_double_bond())
      double_bonds++;
    else
      triple_bonds++;
  }

  int hcount = a.implicit_hydrogens();

  formal_charge_t fc = a.formal_charge();

  if (is_aromatic)
  {
    if (0 == fc && 0 == hcount && 2 == ncon && 2 == aromatic_bonds)   // [o](:*):*
      return 13.14;
  }
  else
  {
    if (0 == fc && 0 == hcount && 2 == ncon && m.in_ring_of_given_size(zatom, 3))     // [O]1-*-*-1     do ring queries first
      return 12.53;
    if (0 == fc && 0 == hcount && 2 == ncon && 2 == single_bonds)     // [O](-*)-*     do ring queries first
      return 9.23;
    if (0 == fc && 0 == hcount && 1 == ncon && 1 == double_bonds)     // [O]=*
      return 17.07;
    if (0 == fc && 1 == hcount && 1 == ncon && 1 == single_bonds)     // [OH]-*
      return 20.23;
    if (-1 == fc && 0 == hcount && 1 == ncon && 1 == single_bonds)    // [O-]-*
      return 23.06;
  }

  return 0;
}

static double
novartis_polar_suface_area_sulphur (Molecule & m,
                                    atom_number_t zatom,
                                    Atom & a,
                                    int is_aromatic)
{
  int ncon = a.ncon();

  int aromatic_bonds = 0;
  int single_bonds = 0;
  int double_bonds = 0;
  int triple_bonds = 0;

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = a[i];

    if (b->is_aromatic())
      aromatic_bonds++;
    else if (b->is_single_bond())
      single_bonds++;
    else if (b->is_double_bond())
      double_bonds++;
    else
      triple_bonds++;
  }

  int hcount = a.implicit_hydrogens();

  formal_charge_t fc = a.formal_charge();

  if (is_aromatic)
  {
    return 0.0;    // the paper doesn't seem to use its own data

    if (0 == fc && 0 == hcount && 2 == ncon && 2 == aromatic_bonds)   // [s](:*):*
      return 28.24;
    if (0 == fc && 0 == hcount && 3 == ncon && 2 == aromatic_bonds && 1 == double_bonds)   // [s](=*)(:*):*
      return 21.70;
  }
  else
  {
    if (0 == fc && 0 == hcount && 2 == ncon && 2 == single_bonds)   // [S](-*)-*
//    return 25.30;    // the paper doesn't seem to use this value, but not using it makes no sense
      return 0.0;
    if (0 == fc && 0 == hcount && 1 == ncon && 1 == double_bonds)   // [S]=*
      return 32.09;
    if (0 == fc && 0 == hcount && 3 == ncon && 2 == single_bonds && 1 == double_bonds)   // [S](-*)(-*)=*
      return 19.21;
    if (0 == fc && 0 == hcount && 4 == ncon && 2 == single_bonds && 2 == double_bonds)   // [S](-*)(-*)(=*)=*
      return 8.38;    // sometimes they don't use this, but I'll use it
    if (0 == fc && 1 == hcount && 1 == ncon && 1 == single_bonds)   // [SH]-*
      return 38.80;
  }

  return 0.0;
}

static double
novartis_polar_suface_area_phosphorus (Molecule & m,
                                       atom_number_t zatom,
                                       Atom & a,
                                       int is_aromatic)
{
  int ncon = a.ncon();

  int aromatic_bonds = 0;
  int single_bonds = 0;
  int double_bonds = 0;
  int triple_bonds = 0;

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = a[i];

    if (b->is_aromatic())
      aromatic_bonds++;
    else if (b->is_single_bond())
      single_bonds++;
    else if (b->is_double_bond())
      double_bonds++;
    else
      triple_bonds++;
  }

  int hcount = a.implicit_hydrogens();

  formal_charge_t fc = a.formal_charge();

  if (0 == fc && 0 == hcount && 3 == ncon && 3 == single_bonds)   // [P](-*)(-*)-*
    return 13.59;
  if (0 == fc && 0 == hcount && 2 == ncon && 1 == single_bonds && 1 == double_bonds)    // [P](-*)=*
    return 34.14;
  if (0 == fc && 0 == hcount && 4 == ncon && 3 == single_bonds && 1 == double_bonds)    // [P](-*)(-*)(-*)=*
    return 9.81;
  if (0 == fc && 1 == hcount && 3 == ncon && 2 == single_bonds && 1 == double_bonds)    // [PH](-*)(-*)=*
    return 23.47;

  return 0.0;
}

static void
quick_prop(Molecule & m,IWString & prop)
{
  prop.resize_keep_storage(0);

  m.compute_aromaticity_if_needed();

  const int matoms = m.natoms();
  const int nr = m.nrings();
  prop <<"\t"<<matoms;
  prop <<"\t"<<nr;
  prop <<"\t"<<m.molecular_weight_ignore_isotopes();
  int ro5_ohnh =0;
  int ro5_on = 0;
  float nvrtspsa = 0.0;
  int halogen = 0; 
  int fcsp3 = 0;
  int ringsys = nr;
  int arring = 0;
  int alring =0;
  int ringatom =0; 
  for (int i = 0; i < matoms; i++)
  {
    Atom * a = const_cast<Atom *>(m.atomi(i));
    atomic_number_t zi = a->atomic_number();
    if (7 == zi || 8 == zi)
    {
      ro5_on ++;
      if (m.hcount(i) > 0)
        ro5_ohnh++;
    }
    if (9 == zi || 17 == zi || 35 == zi || 53 == zi)
      halogen ++;
    if (6 == zi && a->ncon() == a->nbonds())
      fcsp3 ++; 
    if ( m.is_ring_atom(i))
      ringatom ++;
    if ( 7 == zi)
      nvrtspsa += novartis_polar_suface_area_nitrogen(m, i, *a, m.is_aromatic(i));
    else if ( 8 == zi)
      nvrtspsa += novartis_polar_suface_area_oxygen(m, i, *a, m.is_aromatic(i));
    else if (16 == zi)
      nvrtspsa += novartis_polar_suface_area_sulphur(m, i, *a, m.is_aromatic(i));
    else if (15 == zi)
      nvrtspsa += novartis_polar_suface_area_phosphorus(m, i, *a, m.is_aromatic(i));
   
  }
  if(nr > 0)
  {
    int * ring_already_done = new_int(nr); std::unique_ptr<int[]> free_ring_already_done(ring_already_done);
    for (int i = 0; i< nr; ++i)
    {
      const Ring * ri = m.ringi(i);
      if (ri->is_aromatic())
        arring++;
      else
        alring++;
      if( ring_already_done[i] )
        continue;
      if (ri->is_fused())
      {

        for (int j = i + 1; j < nr; j++)
        {
          const Ring * rj = m.ringi(j);
          if (rj->fused_system_identifier() == ri->fused_system_identifier())
          {
            ring_already_done[j] = 1;
            ringsys--;
          }
        }
      }
    }
  }
  prop <<"\t"<<ro5_ohnh; 
  prop <<"\t"<<ro5_on;
  prop <<"\t"<<nvrtspsa;
  prop <<"\t"<<halogen;
  prop <<"\t"<<static_cast<float>(fcsp3)/static_cast<float>(matoms);
  prop <<"\t"<<ringsys;
  prop <<"\t"<<arring;
  prop <<"\t"<<alring;
  prop <<"\t"<<ringatom;
  m.reduce_to_largest_fragment();
  prop <<"\t"<<m.longest_path();
}

namespace {
  typedef std::string String;


/*
  Oct 2015. Add the ability to cap attachment points
*/

  Molecule capping_group;

// We pre-compute the Hydrogen molecular properties. Note that these will change if
// we are doing capping

  String hydrogen_smiles("H");
  String hydrogen_molecular_properties("\t0\t0\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0");

  void adjust_hydrogen_data (Molecule & m)    // we have a capping group
  {
    IWString tmp = m.unique_smiles();
    hydrogen_smiles = tmp.c_str();

    tmp.resize_keep_storage(0);

    quick_prop(m, tmp);
    hydrogen_molecular_properties = tmp.c_str();
//  cerr << "hydrogen_smiles '" << hydrogen_smiles << "'\n";

    return;
  }

  class GroupR {
  public:
    explicit GroupR();
    void push_molecule(const String& smiles, const String& name);
    int get_current_molecule_id() const;

    void push_scaffold(int qid, const String& smiles);
    void push_group(const String& smiles, const String& props);

    int  debug_print(std::ostream &) const;

    void close();

  private:
    std::deque<std::deque<std::pair<String, String> > > _groups;    //the pair consists of smiles and properties
    std::deque<std::pair<String, String> > _current_groups;
    std::deque<std::pair<String, String> > _molecules;        // smiles and name
    std::deque<int> _scaffolds_per_mol;
    std::deque<std::pair<int, String> > _scaffolds;           // query id and smiles
    friend std::ostream& operator<<(std::ostream& stream, const GroupR&);
  };

  GroupR::GroupR() 
  { }

  int GroupR::debug_print(std::ostream & output) const
  {
    output << "GroupR::debug_print:has " << _groups.size() << " groups\n";
    return 1;
  }

  void GroupR::push_molecule(const String& smiles, const String& name) {
    _molecules.push_back(std::make_pair(smiles, name));
    _scaffolds_per_mol.push_back(0);
//  cerr << "push_molecule:got " << smiles << endl;
  }

  int GroupR::get_current_molecule_id() const {
    return _molecules.size()-1;
  }

  void GroupR::push_scaffold(int query_id, const String& smiles) {
    if (!_scaffolds.empty()) {
      _groups.push_back(std::deque<std::pair<String, String> >());
      _groups.back().swap(_current_groups);
    }
    _scaffolds.push_back(std::make_pair(query_id, smiles));
    ++_scaffolds_per_mol.back();
//  cerr << "push_scaffold:got " << smiles << endl;
  }
  
  void GroupR::push_group(const String& smiles, const String& props) {
    _current_groups.push_back(std::make_pair(smiles,props));
//  cerr << "push_group:got " << smiles << " cmp " << _current_groups.back().first << endl;
  }

  void GroupR::close() {
    if (!_scaffolds.empty()) {
      _groups.push_back(std::deque<std::pair<String, String> >());
      _groups.back().swap(_current_groups);
    }
  }

  String replace(const String& source, const String& pattern, const String& data) {
    String::size_type pos1 = source.find(pattern);
    if (pos1 != String::npos) {
      return source.substr(0, pos1) + data + source.substr(pos1+pattern.size());
    }
    else {
      return source;
    }
  }

  String get_scaffold(const String& smiles, const std::set<int>& r2idx) {
    assert(!r2idx.empty());

    Molecule mm;
    mm.build_from_smiles(smiles.c_str());
    cerr << " scaffold smiles " << mm.smiles() << endl;
    int natoms = mm.natoms();

    Set_of_Atoms discard_atoms;
    for (int i = 0; i < natoms; ++i) {
      if (r2idx.find(i) == r2idx.end()) {
        for (int j = 0; j < natoms; ++j) {
          const Atom* atom_j = &mm.atom(j);
          if (atom_j->isotope() == static_cast<isotope_t>(i+1) &&
              atom_j->atomic_number() == 6 && atom_j->size() == 1) {
            discard_atoms.add(j);
          }
        }
      }
    }
    mm.remove_atoms(discard_atoms);
    cerr << "Didscarded " << discard_atoms << " to yield " << mm.smiles() << endl;

    IWString ss = mm.unique_smiles();
    String result = ss.c_str();
    int count_r = 0;
    for (std::set<int>::const_iterator i = r2idx.begin(); i != r2idx.end(); ++i) 
    {
      char buf2[64];
      sprintf(buf2, "[%dCH3]", (*i)+1);
      char buf3[64];
      sprintf(buf3, "[R%d]", ++count_r);
      result = replace(result, buf2, buf3);
    }
    IWString tmp("\t");
    tmp += (mm.natoms() - count_r);
    result+=tmp.c_str();
//  cerr << "Returning scaffold " << result << endl;
    return result;
  }
  
  std::ostream& operator<<(std::ostream& stream, const GroupR& group) {

    int offset = 0;
    std::vector<std::vector<bool> > bits;
//  cerr << "group._molecules.size " << group._molecules.size() << " group._scaffolds_per_mol.size " << group._scaffolds_per_mol.size() << " group._scaffolds.size " << group._scaffolds.size() << endl;
    for (unsigned int j = 0; j < group._molecules.size(); ++j) 
    {
      for (int i = 0; i < group._scaffolds_per_mol[j]; ++i, ++offset) {
        unsigned int qid = group._scaffolds[offset].first;
//      cerr << " j = " << j << " i = " << i << " qid " << qid << " bits.size " << bits.size()  << " offset " << offset << endl;
        if (bits.size() <= qid) {
          bits.resize(qid+1);
        }
//      cerr << " bits[qid].size " << bits[qid].size() << endl;
        cerr.flush();
        
        if (bits[qid].empty()) {
          bits[qid].resize(group._groups[offset].size(), false);
        }
        else {
          assert (bits[qid].size() == group._groups[offset].size());
        }

        for (unsigned int k = 0; k < group._groups[offset].size(); ++k) {
          if (group._groups[offset][k].first != hydrogen_smiles) {
            bits[qid][k] = true;
          }
        }
      }
    }
    
    int max_r = 0;
    for (unsigned int i = 0; i < bits.size(); ++i) {
      int n = 0;
      for (unsigned int j = 0; j < bits[i].size(); ++j) {
        if (bits[i][j]) {
          ++n;
        }      
      }
      if (n > max_r) {
        max_r = n;
      }
    }

    stream << "CompoundName\tSMILES\t";
    stream << "Scaffold\tScaffold_natoms"; // TODO: Dump Scaffold for every molecule
    for (int j = 0; j < max_r; ++j) {
      stream << "\tR" << j+1 << "_SMILES" ;
      if(!add_properties_to_groups)
        continue;
      for(int k =0; k < prop_number; ++k)
        stream <<"\tR" << j+1 << "_"<<prop_header[k];
    }
    stream << "\n";

    offset = 0;
    for (unsigned int j = 0; j < group._molecules.size(); ++j) {
      for (int i = 0; i < group._scaffolds_per_mol[j]; ++i, ++offset) {
        int qid = group._scaffolds[offset].first;
      
        std::set<int> r2idx; 
        for (unsigned int k = 0; k < bits[qid].size(); ++k) {
          if (bits[qid][k]) {
            r2idx.insert(k);
          }
        }
        
        stream << group._molecules[j].second << "\t" << group._molecules[j].first << "\t";
        stream << get_scaffold(group._scaffolds[offset].second, r2idx); // TODO: Dump Scaffold for every molecule
        for (std::set<int>::const_iterator k = r2idx.begin(); k != r2idx.end(); ++k) 
        {
          stream << "\t" << group._groups[offset][*k].first << group._groups[offset][*k].second;
        }
        stream << "\n";
      }
    }
    return stream;
  }

  GroupR* gp_group_r = nullptr;
  String gs_rgroup_file;
}

static Chemical_Standardisation chemical_standardisation;

static Molecule_Output_Object stream_for_substituents;
static Molecule_Output_Object stream_for_molecules_not_matching;

static IWString_and_File_Descriptor stream_for_jw;
static IWString_and_File_Descriptor stream_for_tb;

static int verbose = 0;

static int reduce_to_largest_fragment = 0;

static int make_implicit_hydrogens_explicit=0;

static int remove_non_single_bond_match =0;
static int remove_ring_bond_match =0;

static int molecules_read = 0;

static int molecules_not_hitting_any_queries = 0;

static extending_resizable_array<int> fragments_per_molecule;

static int write_to_stdout = 1;

static resizable_array_p<Substructure_Hit_Statistics> queries;

static int nq = 0;

static extending_resizable_array<int> matches_per_molecule;

static int max_atoms_in_a_fragment = std::numeric_limits<int>::max();

static int fragments_with_too_many_atoms = 0;

static int isotope_for_ring_closures = 0;

static int write_no_match_records = 1;

/*
  We often need to mark atoms as being part of a substituent,
  but we need to exclude those that are part of the embedding.
  To avoid having two separate arrays, we use an unusual
  isotopic value to designate those in the embedding, or
  other atoms that should be excluded from the substituent
*/

#define EMBEDDING_ISO 97531

/*
  Jun 2004. by default, the "scaffold" atom where the
  substituent is attached is included with the substituent.
  That seems wrong, so make that behaviour optional, but
  still the default
*/

static int include_anchor_atom_with_substituent = 1;

/*
  What to do with multiple hits to a query
*/

static int take_first_of_multiple_hits = 0;

static int ignore_queries_hitting_multiple_times = 0;

static int ignore_queries_not_matching = 0;

static int break_after_first_query_match = 0;

static int ok_for_no_queries_to_match = 0;

static int process_all_matching_queries = 0;

static IWString missing_value('.');

static int write_parent_molecule = 1;

/*
  What do we do if our matched atom is part of a ring - and
  there is just one matched atom in that ring
*/

static int treat_rings_as_substitutents = 0;

/*
  Jun 2012. The current implementation allows only one isotopic label designating an attachment point. 
  And it imposes restrictions on where the molecule can be substituted. I suspect nobody is using that
  capability, and it seems wrong. But I am not going to change it. Instead, get better behaviour
  with an option
*/

static int multiple_isotopic_attachment_points = 0;

/*
  At each point of substitution, we need a hash of the fragments found at that point
*/

class Fragments_Found : public IW_STL_Hash_Map_int
{
  private:
  public:

    int extra(const IWString &);

    int report(const IWString &, std::ostream &) const;
};

int
Fragments_Found::extra (const IWString & s)
{
  IW_STL_Hash_Map_int::iterator f = find(s);

  if (f == end())
    insert(std::pair<IWString, int>(s, 1));
  else
    (*f).second++;

  return 1;
}

int
Fragments_Found::report (const IWString & s,
             std::ostream & output) const
{
  output << size() << " attachments found\n";

  for (IW_STL_Hash_Map_int::const_iterator i = begin(); i != end(); ++i)
  {
    output << (*i).first << " " << (*i).second << " occurrences";
    if (s.length())
      output << ' ' << s;
    output << "\n";
  }

  return output.good();
}

/*
  The matched atoms from which to grow is a set of either single atom,
  and/or pairs of atoms - indicating whatever is in between the
  matched atoms.  If we are dealing with a pair, then the entry in the
  _a2 array will be valid
*/

class Matched_Atoms_From_Which_to_Grow : public resizable_array<int>
{
  private:
    int _n;

    atom_number_t * _a1;
    atom_number_t * _a2;

  public:
    Matched_Atoms_From_Which_to_Grow();
    ~Matched_Atoms_From_Which_to_Grow();

    int number_join_points() const { return _n;}

    int build(const const_IWSubstring &);

    int add(const Set_of_Atoms &);

    int is_single_attachment_point(int i) const;

    int set_single_attachment_point (int);
    int set_double_attachment_point (int, int);

    int a1 (int i) const { return _a1[i];}
    int a2 (int i) const { return _a2[i];}
};

/*
  For each query, we can have any number of matched atoms at which we
  look for substituents
*/

static Matched_Atoms_From_Which_to_Grow * matched_atoms_from_which_to_grow = nullptr;

Matched_Atoms_From_Which_to_Grow::Matched_Atoms_From_Which_to_Grow()
{
  _n = 0;

  _a1 = nullptr;
  _a2 = nullptr;

  return;
}

Matched_Atoms_From_Which_to_Grow::~Matched_Atoms_From_Which_to_Grow()
{
  if (nullptr != _a1)
    delete [] _a1;

  if (nullptr != _a2)
    delete [] _a2;

  return;
}

int
Matched_Atoms_From_Which_to_Grow::is_single_attachment_point(int i) const
{
  if (nullptr == _a2)
    return 1;

  return _a2[i] < 0;
}

int 
Matched_Atoms_From_Which_to_Grow::build (const const_IWSubstring & s)
{
  assert(nullptr == _a1);

  if (0 == s.length())    // huh
    return 0;

  _n = s.nwords(',');
  _a1 = new int[_n];
  _a2 = new_int(_n, -1);

  int i = 0;
  const_IWSubstring token;
  int ndx = 0;

  while (s.nextword(token, i, ','))
  {
    const_IWSubstring s1, s2;
    if (token.split(s1, '-', s2))
    {
      if (! s1.numeric_value(_a1[ndx]) || _a1[ndx] < 0 ||
        ! s2.numeric_value(_a2[ndx]) || _a2[ndx] < 0 ||
        _a1[ndx] == _a2[ndx])
      {
        cerr << "Matched_Atoms_From_Which_to_Grow::build:invalid pair specification '" << s << "'\n";
        return 0;
      }

      add(_a1[ndx]);
      add(_a2[ndx]);
    }
    else if (! token.numeric_value(_a1[ndx]) || _a1[ndx] < 0)
    {
      cerr << "Matched_Atoms_From_Which_to_Grow::build:invalid single atom specification '" << token << "'\n";
      return 0;
    }
    else
      resizable_array<int>::add(_a1[ndx]);

    ndx++;
  }

  assert (ndx == _n);

  return ndx;
}

int
Matched_Atoms_From_Which_to_Grow::add (const Set_of_Atoms & s)
{
  assert(0 == _n);

  _n = s.number_elements();
  _a1 = new int[_n];
  _a2 = new_int(_n, -1);

  for (int i = 0; i < _n; i++)
  {
    resizable_array<int>::add(s[i]);
    _a1[i] = s[i];
  }

  return _n;
}

int
Matched_Atoms_From_Which_to_Grow::set_single_attachment_point(int s)
{
  assert(nullptr == _a1);

  _n = 1;
  _a1 = new int[1];
  _a2 = new int[1];

  assert(nullptr != _a1 && nullptr != _a2);

  _a1[0] = s;
  _a2[0] = -1;

  return 1;
}

int
Matched_Atoms_From_Which_to_Grow::set_double_attachment_point(int s1, int s2)
{
  assert(nullptr == _a1);

  _n = 1;
  _a1 = new int[1];
  _a2 = new int[1];

  assert(nullptr != _a1 && nullptr != _a2);

  _a1[0] = s1;
  _a2[0] = s2;

  return 1;
}


/*
  We need a means of identifying each point of substitution. We use
  a string with the query number followed by a period and then the matched
  atom number
*/

IW_STL_Hash_Map<IWString, Fragments_Found *> fragments_found_at_substitution_point;

static int compute_unique_smiles_of_fragments = 0;

static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on

  cerr << "  -q <query>     identify features as queries\n";
  cerr << "  -s <smarts>    identify features as smarts\n";
  cerr << "  -e <smiles>    identify features as smiles (consider explicit hydrogens)\n";
  cerr << "  -m <fname>     queries are isotopically labelled molecules in <fname>\n";
  cerr << "  -h             any number of isotopes allowed\n";
  cerr << "  -z first       take the first match when multiple matches to a query\n";
  cerr << "  -z each        process each of multiple matches to a query\n";
  cerr << "  -z ignmmatch   ignore any query that matches multiple times\n";
  cerr << "  -z ignore      ignore queries that do not match the molecule\n";
  cerr << "  -z oknomatch   ignore molecules that are not matched by any query\n";
  cerr << "  -j <number>    look for substituents at matched atom <number>\n";
  cerr << "                 for matched atoms in same query use '-j 0,4'\n";
  cerr << "  -r             produce substituents that include rings\n";
  cerr << "  -S <fname>     write molecules and substituents to <fname>\n";
  cerr << "  -Z <fname>     write non reacting molecules to <fname>\n";
  cerr << "  -C <number>    max number of atoms in a sidechain\n";
  cerr << "  -c <iso>       isotope to mark ring closing substituents\n";
  cerr << "  -n             suppress normal output (useful for database loads)\n";
  cerr << "  -o <type>      output type(s)\n";
  cerr << "  -k             discard symmetry equivalent matches\n";
  cerr << "  -p             compute the unique smiles of substituents and make counts\n";
  cerr << "  -x             exclude the scaffold atom from the substituent\n";
  cerr << "  -y             do NOT write the parent molecule to the output stream\n";
  cerr << "  -J <fname>     tabular output with unique smiles\n";
  cerr << "  -b             break after first query matches\n";
  cerr << "  -u             do NOT write NO_MATCH data\n";
  cerr << "  -d <fname>     dump substitions r-group into file.\n";
  cerr << "  -B             substitions are single-bond connected to queries.\n";
  cerr << "  -R             substitions are not in same ring as queries.\n";
  cerr << "  -w             calculate properties for each r-group of -d\n";
  cerr << "  -E ...         standard Element options, enter '-E help' for details\n";
  cerr << "  -l             reduce to largest fragment\n";
  cerr << "  -i <type>      input type\n";
  (void) display_standard_aromaticity_options(cerr);
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -v             verbose output\n";

  exit(rc);
}

static int
read_capping_group (Molecule & m,
                    iwstring_data_source & input)
{
  const_IWSubstring buffer;
  if (! input.next_record(buffer))
  {
    cerr << "read_capping_group:cannot read\n";
    return 0;
  }

  if (! m.build_from_smiles(buffer))
  {
    cerr << "read_capping_group:cannot interpret smiles '" << buffer << "'\n";
    return 0;
  }

  return m.natoms();
}

static int
read_capping_group (Molecule & m,
                     const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "read_capping_group::cannot open '" << fname << "'\n";
    return 0;
  }

  return read_capping_group(m, input);
}

static int
do_compute_unique_smiles_of_fragments (const IWString & zkey,
                                       const IWString & usmi)
{
  IW_STL_Hash_Map<IWString, Fragments_Found *>::iterator f = fragments_found_at_substitution_point.find(zkey);
  if (f == fragments_found_at_substitution_point.end())
  {
    Fragments_Found * tmp = new Fragments_Found;
    tmp->extra(usmi);
    fragments_found_at_substitution_point[zkey] = tmp;
  }
  else
  {
    Fragments_Found * ff = (*f).second;
    ff->extra(usmi);
  }

  return 1;
}

#define DEBUG_IDENTIFY_ATOMS_IN_FRAGMENT

static int
identify_atoms_in_fragment (Molecule & m,
              atom_number_t zatom,
              atom_number_t previous_atom,
              int * in_fragment,
              atom_number_t anchor)
{
  in_fragment[zatom] = 1;

  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();

  int rc = 1;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (j == previous_atom)
      continue;

#ifdef DEBUG_IDENTIFY_ATOMS_IN_FRAGMENT
    cerr << "identify_atoms_in_fragment:atom " << zatom << " to " << j << " visited " << in_fragment[j] << endl;
#endif

    if (in_fragment[j])
      continue;

    if (EMBEDDING_ISO == in_fragment[j])    // we've found a ring closing substituent - ring back to the embedding
    {
      if (INVALID_ATOM_NUMBER == previous_atom)    // just passing around the embedding
        continue;

#ifdef DEBUG_IDENTIFY_ATOMS_IN_FRAGMENT
      cerr << "From " << anchor << " found ring closure to atom " << j << " via " << previous_atom << endl;
#endif

      if (isotope_for_ring_closures)
        m.set_isotope(zatom, isotope_for_ring_closures);
      continue;
    }

    rc += identify_atoms_in_fragment(m, j, zatom, in_fragment, anchor);
  }

  return rc;
}

static void
set_isotopes(Molecule & m,
             const Set_of_Atoms & s,
             int iso)
{
  m.set_isotope(s, iso);

  return;
}


static void
identify_unmatched_connections (Molecule & m,
                                const int * in_fragment,
                                atom_number_t zatom,
                                Set_of_Atoms & unmatched_connections,
                                int iso)
{
  const Atom * ai = m.atomi(zatom);

  int acon = ai->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = ai->other(zatom, i);

    if (in_fragment[j])
      continue;

    unmatched_connections.add(j);
    m.set_isotope(j, iso);
  }

  return;
}

static int
identify_ring_between(Molecule & m,
                      atom_number_t previous_atom,
                      atom_number_t zatom,
                      atom_number_t have_ring_if_we_encounter,
                      int * in_fragment,
                      int & found_ring)
{
  in_fragment[zatom] = 1;

  const Atom * ai = m.atomi(zatom);

  int rc = 1;

  int acon = ai->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = ai->other(zatom, i);

    if (j == previous_atom)
      continue;

    if (j == have_ring_if_we_encounter)
      found_ring = 1;

    if (in_fragment[j])
      continue;

    rc += identify_ring_between(m, zatom, j, have_ring_if_we_encounter, in_fragment, found_ring);
  }

  return rc;
}

int
strip_to_matched_atoms(Molecule & parent,
                       const Set_of_Atoms & embedding,
                       Molecule & subset)
{
  subset.add_molecule(&parent);

  const int n = embedding.number_elements();

  for (int i = 0; i < n; ++i)
  {
    subset.set_isotope(embedding[i], i + 1);
  }

  for (int i = subset.natoms() - 1; i >= 0; --i)
  {
    if (0 == subset.isotope(i))
      subset.remove_atom(i);
  }

  return 1;
}

template <typename T>
int
write_thibault_data(Molecule & parent,
                    const Set_of_Atoms & embedding,
                    const int query_number,
                    const int ndx,
                    const int * in_fragment,
                    Molecule & matched_atoms,
                    T & output)
{
//cerr << "write_thibault_data begin " << parent.smiles() << endl;
  Molecule mcopy(parent);

  const int n = embedding.number_elements();

  for (int i = 0; i < n; ++i)
  {
    mcopy.set_isotope(embedding[i], i + 1);
  }

//cerr << "After isotope assignment " << mcopy.smiles() << endl;

  for (int i = 0; i < n; ++i)
  {
    const atom_number_t j = embedding[i];

    const Atom * a = mcopy.atomi(j);

    const int acon = a->ncon();

    int connection_to_fragment = 0;

    for (int k = 0; k < acon; ++k)
    {
      const atom_number_t l = a->other(j, k);

      if (! in_fragment[l])
        continue;

      if (mcopy.isotope(l))    // perhaps a ring???
        continue;

      mcopy.set_isotope(l, i + 1);
      connection_to_fragment++;
    }

//  if (0 == connection_to_fragment)
//    mcopy.set_isotope(j, 0);
  }

  Molecule subset;
  mcopy.create_subset(subset, in_fragment);
  mcopy.remove_atoms(in_fragment);

//output << parent.smiles() << ' ' << mcopy.smiles() << ' ' << parent.name() << ' ' << subset.smiles() << ' ' << query_number << '.' << ndx << '\n';
  output << matched_atoms.smiles() << ' ' << parent.name() << ' ' << subset.smiles() << ' ' << query_number << '.' << ndx << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static int
write_parent_if_needed_and_subset (Molecule & parent,
                                   Molecule & subset,
                                   int & parent_molecule_written,
                                   Molecule_Output_Object & output)
{
  if (! parent_molecule_written)
  {
    output.write(parent);
    parent_molecule_written = 1;
  }

  return output.write(subset);
}

/*
  We need to deal with the case where two single attachment points may be part
  of the same ring. Note that there are severe problems with this because, for
  example, we can't form the unique smiles of a partial aromatic ring.

  This is pretty ugly because the inner loop is so large, but too many
  arguments to pass around if I try to break it up
*/

static int
process_join_points_that_form_rings(Molecule & m,
                                    int * already_done,
                                    int query_number,
                                    const Set_of_Atoms & embedding,
                                    int * in_fragment,
                                    IWString & data_for_jibo,
                                    int & parent_molecule_written)
{
  int matoms = m.natoms();

  const Matched_Atoms_From_Which_to_Grow & ma = matched_atoms_from_which_to_grow[query_number];

  int n = ma.number_join_points();

  for (int i = 0; i < n; i++)
  {
    if (already_done[i])
      continue;

#ifdef DEBUG_PROCESS_JOIN_POINTS_THAT_FORM_RINGS
    cerr << "Can we process " << i << " ? " << ma.is_single_attachment_point(i) << endl;
#endif

    if (! ma.is_single_attachment_point(i))    // too hard to handle these
      continue;

    atom_number_t j1 = embedding[ma.a1(i)];

    for (int k = i + 1; k < n; k++)
    {
      if (already_done[k])
        continue;

#ifdef DEBUG_PROCESS_JOIN_POINTS_THAT_FORM_RINGS
      cerr << "Can we process " << k << " ? " << ma.is_single_attachment_point(k) << endl;
#endif

      if (! ma.is_single_attachment_point(k))
        continue;

      atom_number_t j2 = embedding[ma.a1(k)];

#ifdef DEBUG_PROCESS_JOIN_POINTS_THAT_FORM_RINGS
      cerr << "Atoms " << j1 << " and " << j2 << " in_same_ring " << m.in_same_ring(j1, j2) << endl;
#endif

      if (! m.in_same_ring(j1, j2))
        continue;

#ifdef DEBUG_PROCESS_JOIN_POINTS_THAT_FORM_RINGS
      cerr << "Possible ring substituent, atoms " << j1 << " to " << j2 << endl;
#endif

//    Just because the anchor points are in the same ring, doesn't mean that the substitutents are

      set_vector(in_fragment, matoms, 0);
      embedding.set_vector(in_fragment, EMBEDDING_ISO);

      int found_ring = 0;
      //    cerr << "Looking for ring between atoms " << j1 << " and " << j2 << endl;
      int fragment_size = identify_ring_between(m, j2, j1, j2, in_fragment, found_ring);
      //    cerr << "Found ring? " << found_ring << endl;

      if (! found_ring)    // only interested in these cases
        continue;

      already_done[i] = 1;
      already_done[k] = 1;

      if (0 == fragment_size && ! include_anchor_atom_with_substituent)   // cannot do anything with this
        continue;

      if (fragment_size > max_atoms_in_a_fragment)
      {
        fragments_with_too_many_atoms++;
        continue;
      }

      Set_of_Atoms unmatched_connections;
      if (! include_anchor_atom_with_substituent)
      {
        identify_unmatched_connections(m, in_fragment, j1, unmatched_connections, query_number + 1);
        identify_unmatched_connections(m, in_fragment, j2, unmatched_connections, query_number + 1);
      }

      //    Now translate back all the items that have been set to EMBEDDING_ISO in the array

      if (include_anchor_atom_with_substituent)
      {
        in_fragment[j1] = 1;
        in_fragment[j2] = 1;
        m.set_isotope(j1, query_number + 1);
        m.set_isotope(j2, query_number + 1);
      }
      else
      {
        in_fragment[j1] = 0;
        in_fragment[j2] = 0;
      }

      for (int l = 0; l < matoms; l++)
      {
        if (EMBEDDING_ISO == in_fragment[l])
          in_fragment[l] = 0;
      }

      IWString subset_name;
      subset_name << m.name() << " Q" << query_number << '.' << ma.a1(i) << '@' << ma.a1(k);

      Molecule subset;
      m.create_subset(subset, in_fragment);
      subset.set_name(subset_name);

      if (include_anchor_atom_with_substituent)
      {
        m.set_isotope(j1, 0);
        m.set_isotope(j2, 0);
      }
      else
        m.set_isotope(unmatched_connections, 0);

      if (stream_for_substituents.active())
        write_parent_if_needed_and_subset(m, subset, parent_molecule_written, stream_for_substituents);

      if (compute_unique_smiles_of_fragments)
      {
        IWString zkey;
        zkey << query_number << '.' << ma.a1(i) << '@' << ma.a1(k);
        do_compute_unique_smiles_of_fragments(zkey, subset.unique_smiles());
      }

      if (stream_for_jw.is_open())
        data_for_jibo << " QQ" << query_number << '.' << ma.a1(i) << '@' << ma.a1(k) << '.' << subset.natoms() << ' ' << subset.unique_smiles();
    }
  }

  return 1;
}

static int
apply_isotope_to_first_atom_in_fragment(Molecule & m,
                    const int * in_fragment,
                    atom_number_t zatom,
                    int iso,
                    atom_number_t & c)
{
  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (in_fragment[j])
    {
      m.set_isotope(j, iso);
      c = j;
      return 1;
    }
  }

  return 0;    // should not happen
}

static int 
identify_atoms_between (Molecule & m,
            atom_number_t zatom1,
            atom_number_t zatom2,
            int * in_fragment,
            int flag)
{
  const int matoms = m.natoms();

  int rc = 0;

  int dmin = m.bonds_between(zatom1, zatom2);

  for (int i = 0; i < matoms; i++)
  {
    if (i == zatom1 || i == zatom2)
      continue;

    int d1 = m.bonds_between(zatom1, i);
    int d2 = m.bonds_between(zatom2, i);

    if (d1 < dmin && d2 < dmin)
    {
      in_fragment[i] = flag;
      rc++;
    }
  }

  return rc;
}

static int
identify_substituent_between (Molecule & m,
                              int query_number,
                              const Set_of_Atoms & embedding,
                              int ndx1,
                              int ndx2,
                              int * in_fragment,
                              IWString & data_for_jibo)

{
  assert (embedding.ok_index(ndx1));
  assert (embedding.ok_index(ndx2));

  atom_number_t zatom1 = embedding[ndx1];
  atom_number_t zatom2 = embedding[ndx2];

  set_vector(in_fragment, m.natoms(), 0);

  int fragment_size = identify_atoms_between(m, zatom1, zatom2, in_fragment, 1);

  if (fragment_size <= 1)
  {
    if (stream_for_jw.is_open())
      data_for_jibo << " QQ" << query_number << '.' << ndx1 << '-' << ndx2 << '.' << 0 << ' ' << '*' << '\n';

    return 0;
  }

  if (fragment_size > max_atoms_in_a_fragment)
  {
    fragments_with_too_many_atoms++;
    return 0;
  }

  if (include_anchor_atom_with_substituent)
  {
    in_fragment[zatom1] = 1;
    in_fragment[zatom2] = 1;
  }

  IWString subset_name;
  subset_name << m.name() << " Q" << query_number << '.' << ndx1 << '-' << ndx2;

  m.set_isotope(zatom1, query_number + 1);
  m.set_isotope(zatom2, query_number + 2);

  // We may need to put isotopic labels on the atoms connected to zatom1 and zatom2

  atom_number_t c1 = INVALID_ATOM_NUMBER;
  atom_number_t c2 = INVALID_ATOM_NUMBER;
  if (! include_anchor_atom_with_substituent)
  {
    apply_isotope_to_first_atom_in_fragment(m, in_fragment, zatom1, query_number + 1, c1);
    apply_isotope_to_first_atom_in_fragment(m, in_fragment, zatom2, query_number + 2, c2);
  }

  Smiles_Information smi_info;

  Molecule subset;
  m.create_subset(subset, in_fragment);
  subset.set_name(subset_name);

  m.set_isotope(zatom1, 0);
  m.set_isotope(zatom2, 0);

  if (INVALID_ATOM_NUMBER != c1)
    m.set_isotope(c1, 0);
  if (INVALID_ATOM_NUMBER != c2)
    m.set_isotope(c2, 0);

  if (stream_for_substituents.active())
    stream_for_substituents.write(subset);

  if (compute_unique_smiles_of_fragments)
  {
    IWString zkey;
    zkey << query_number << '.' << ndx1 << '-' << ndx2;
    do_compute_unique_smiles_of_fragments(zkey, subset.unique_smiles());
  }

  if (stream_for_jw.is_open())
    data_for_jibo << " QQ" << query_number << '.' << ndx1 << '-' << ndx2 << '.' << subset.natoms() << ' ' << subset.unique_smiles();
  
  return 1;
}

static void
add_capping_group (Molecule & m, const Molecule & capping_group)
{
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    if (0 == m.isotope(i))
      continue;

    const int current_atoms = m.natoms();

    m.add_molecule(&capping_group);

    m.add_bond(i, current_atoms, SINGLE_BOND);
  }

  return;
}

/*
  this is buggy, we need to differentiate the initial embedding from other atoms that
  may be encountered.
*/

static int
identify_substituent (Molecule & m,
                      const int query_number,
                      const Set_of_Atoms & embedding,
                      const int ndx,
                      int * in_fragment,
                      IWString & data_for_jibo,
                      int & parent_molecule_written,
                      Molecule & matched_atoms)
{
  assert(embedding.ok_index(ndx));

  const atom_number_t zatom = embedding[ndx];

#ifdef DEBUG_IDENTIFY_SUBSTITUENT
  cerr << "Substituent starting with atom " << zatom << ' ' << m.smarts_equivalent_for_atom(zatom) << endl;
#endif

  set_vector(in_fragment, m.natoms(), 0);

  embedding.set_vector(in_fragment, EMBEDDING_ISO);

  Set_of_Atoms unmatched_connections;

  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    const atom_number_t j = a->other(zatom, i);

    if (EMBEDDING_ISO == in_fragment[j])   // part of the embedding
      ;
    else if (! m.in_same_ring(zatom, j))
      unmatched_connections.add(j);
    else if (treat_rings_as_substitutents)
      unmatched_connections.add(j);

    //  cerr << "Atom " << j << " status " << in_fragment[j] << endl;
  }

  cerr << "Around atom " << zatom << ' ' << m.smarts_equivalent_for_atom(zatom) << " unmatched " << unmatched_connections << endl;

  int fragment_size;
  if (unmatched_connections.number_elements() > 0)
    fragment_size = identify_atoms_in_fragment(m, zatom, INVALID_ATOM_NUMBER, in_fragment, zatom);
  else
    fragment_size = 0;   // possibly Hydrogen

  if (verbose > 1)
    cerr << "atom " << zatom << " fragment_size " << fragment_size << endl;

  if (fragment_size <= 1)
  {
    if (stream_for_jw.is_open())
    {
      data_for_jibo << " QQ" << query_number << '.' << ndx << '.' << 0 << ' ';
      if (include_anchor_atom_with_substituent)
        data_for_jibo << '[' << (query_number+1) << m.atomic_symbol(zatom) << ']';
      data_for_jibo << 'H';
    }
    
    if (gp_group_r) {
      IWString prop;
      if (add_properties_to_groups)
        prop << hydrogen_molecular_properties;
      gp_group_r->push_group("H", prop.c_str());

#ifdef BREAKS_THINGS
      IWString s;
      if (include_anchor_atom_with_substituent)
        s << "[1" << m.atomic_symbol(zatom) << "]H";
      else
        s = hydrogen_smiles;
      gp_group_r->push_group(s.c_str(), prop.c_str());
#endif
//    cerr << "Added H substituent " << hydrogen_smiles << endl;
    }

    return 0;
  }

  if (fragment_size > max_atoms_in_a_fragment)
  {
    fragments_with_too_many_atoms++;
    return 0;
  }

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (EMBEDDING_ISO == in_fragment[i])
      in_fragment[i] = 0;
  }

  if (! include_anchor_atom_with_substituent)
    in_fragment[zatom] = 0;

  if (stream_for_tb.is_open())
    write_thibault_data(m, embedding, query_number, ndx, in_fragment, matched_atoms, stream_for_tb);

  IWString subset_name;
  subset_name << m.name() << " Q" << query_number << '.' << ndx;

  if (include_anchor_atom_with_substituent)
    m.set_isotope(zatom, query_number + 1);
  else
    set_isotopes(m, unmatched_connections, query_number + 1);

  Smiles_Information smi_info;

  Molecule subset;
  m.create_subset(subset, in_fragment);
  subset.set_name(subset_name);

  if (include_anchor_atom_with_substituent)
    m.set_isotope(zatom, 0);
  else
    set_isotopes(m, unmatched_connections, 0);

  if (capping_group.natoms())
    add_capping_group(subset, capping_group);

  if (stream_for_substituents.active())
    write_parent_if_needed_and_subset(m, subset, parent_molecule_written, stream_for_substituents);

  if (compute_unique_smiles_of_fragments)
  {
    IWString zkey;
    zkey << query_number << '.' << ndx;
    do_compute_unique_smiles_of_fragments(zkey, subset.unique_smiles());
  }

  if (stream_for_jw.is_open())
    data_for_jibo << " QQ" << query_number << '.' << ndx << '.' << subset.natoms() << ' ' << subset.unique_smiles();

  if (gp_group_r) {
    IWString ss = subset.unique_smiles();
    IWString prop;
     if (add_properties_to_groups)
      quick_prop(subset,prop);
    gp_group_r->push_group(ss.c_str(), prop.c_str());
  }

  return 1;
}

static int
same_atoms(const Matched_Atoms_From_Which_to_Grow & ma,
           const int query_number,
           const Set_of_Atoms & s1,
           const Set_of_Atoms & s2)
{
  int n = ma.number_join_points();

  for (int i = 0; i < n; i++)
  {
    int j = ma.a1(i);

    atom_number_t k1 = s1[j];
    atom_number_t k2 = s2[j];

    if (k1 != k2)    // different atoms hit, embeddings are not the same
      return 0;
  }

  return 1;   // if we get to here, they are identical
}

static int
remove_embeddings_that_result_in_duplicate_atoms(const Matched_Atoms_From_Which_to_Grow & ma,
                          int query_number,
                          Substructure_Results & sresults)
{
  if (! ma.is_single_attachment_point(query_number))   // not bothering to handle these - complex
    return 0;

  int nhits = sresults.number_embeddings();

  if (nhits < 2)
    return 1;

  for (int i = nhits - 1; i >= 0; i--)
  {
    const Set_of_Atoms * si = sresults.embedding(i);

    int remove_embedding_i = 0;

    for (int j = i - 1; j >= 0; j--)
    {
      const Set_of_Atoms * sj = sresults.embedding(j);

      if (same_atoms(ma, query_number, *si, *sj))
        remove_embedding_i = 1;
    }

    if (remove_embedding_i)
      sresults.remove_embedding(i);
  }

  return sresults.number_embeddings();
}

static int
remove_embeddings ( Molecule & m, Substructure_Results & sresults)
{
  int nhits = sresults.number_embeddings();
  for (int i = nhits - 1; i >= 0; i--)
  {
    const Set_of_Atoms * si = sresults.embedding(i);
    int remove_embedding_i = 0;
    for (unsigned int j = 0; j<si->size()&& !remove_embedding_i ;j++)
    {
      int ai = (*si)[j];
      const Atom * a = m.atomi(ai);
      int acon = a->ncon();  
      for (int k = 0; k < acon; k++)
      {
        atom_number_t bi = a->other(ai, k);
        if( si->contains(bi) )
          continue;
        if( remove_non_single_bond_match && !a->bond_to_atom(ai,bi)->is_single_bond())
        {
          remove_embedding_i = 1;
          break;
        }
        //cerr << remove_ring_bond_match << " "<<ai<<" "<<bi<<"\n";
        if( remove_ring_bond_match && m.in_same_ring(ai,bi))
        {
          remove_embedding_i = 1;
          break;
        }
      }
    }  
    if (remove_embedding_i)
      sresults.remove_embedding(i);
  }
  return sresults.number_embeddings();
}

static int
identify_substituents (Molecule & m,
                       const int query_number,
                       const Set_of_Atoms & embedding,
                       int * in_fragment,
                       IWString & data_for_jibo,
                       int & parent_molecule_written)
{
  const Matched_Atoms_From_Which_to_Grow & ma = matched_atoms_from_which_to_grow[query_number];

  int rc = 0;

  const int n = ma.number_join_points();

  int array_size = n;
  if (0 == n)
    array_size = embedding.number_elements();

  int * already_done = new_int(array_size); std::unique_ptr<int[]> free_already_done(already_done);

  //#define DEBUG_IDENTIFY_SUBSTITUENTS
#ifdef DEBUG_IDENTIFY_SUBSTITUENTS
  cerr << "Examining " << n << " join points, embedding " << embedding << endl;
#endif

  if (n > 1)
    process_join_points_that_form_rings(m, already_done, query_number, embedding, in_fragment, data_for_jibo, parent_molecule_written);

  Molecule matched_atoms;

  if (stream_for_tb.is_open())
    strip_to_matched_atoms(m, embedding, matched_atoms);

  if (n > 0)    // just the atoms specified
  {
    for (int i = 0; i < n; i++)
    {
      if (already_done[i])
        continue;

      if (ma.is_single_attachment_point(i))
      {
        const int j = ma.a1(i);

        rc += identify_substituent(m, query_number, embedding, j, in_fragment, data_for_jibo, parent_molecule_written, matched_atoms);
      }
      else
      {
        const int j1 = ma.a1(i);
        const int j2 = ma.a2(i);
        rc += identify_substituent_between(m, query_number, embedding, j1, j2, in_fragment, data_for_jibo);
      }
    }
  }
  else    // All connections to matched atoms
  {
    const int n = embedding.number_elements();
    for (int i = 0; i < n; i++)
    {
      if (already_done[i])
        continue;

      rc += identify_substituent(m, query_number, embedding, i, in_fragment, data_for_jibo, parent_molecule_written, matched_atoms);
    }
  }

  return rc;
}

#ifdef OLD_VERSION_
static int
gp_add_scaffold_iw(const Molecule & m,
                   const int query_number,
                   const Set_of_Atoms & presult)
{
  if (nullptr == gp_group_r)
    return 0;

  const int matoms = m.natoms();

  int * data = new_int(matoms); std::unique_ptr<int[]> free_data(data);
  presult.set_vector(data, 1);

  int * atomInQuery = new_int(matoms, -1); std::unique_ptr<int[]> free_atomInQuery(atomInQuery);

  for (std::size_t j = 0; j < presult.size(); ++j) {
    atomInQuery[presult[j]] = j;
  }

  int * aam = new_int(matoms); std::unique_ptr<int[]> free_aam(aam);    // index of M atom in FRAG

  Molecule frag;
  m.create_subset(frag, data, 1, aam);

  const int natoms_in_frag = frag.natoms(); 

  int * atomInPattern = new_int(natoms_in_frag, -1); std::unique_ptr<int[]> free_atomInPattern(atomInPattern);   // index of atom in FRAG in M

  for (auto j = 0; j < matoms; ++j) {
    if (aam[j] != -1) {
      atomInPattern[aam[j]] = j;
    }
  }

  for (int j = 0; j < natoms_in_frag; ++j)
  {
    Atom* atom_j = new Atom("C");

    assert(atomInPattern[j] >= 0);
    assert(atomInQuery[atomInPattern[j]] >= 0);

    atom_j->set_isotope(atomInQuery[atomInPattern[j]] + 1);
    frag.add(atom_j);
    frag.add_bond(j, natoms_in_frag + j, SINGLE_BOND);
  }

//Molecule mm(m);
//cerr << "From " << mm.smiles() << " generate fragment " << frag.smiles() << endl;

  IWString ss = frag.unique_smiles();
  gp_group_r->push_scaffold(query_number, ss.c_str());

  return 1;
}
#endif

static int
gp_add_scaffold_iw_v2(const Molecule & m,
                      const int query_number,
                      const Set_of_Atoms & presult)
{
  if (nullptr == gp_group_r)
    return 0;

  const int matoms = m.natoms();

#ifdef NEW_SIMPLER_VERSION
  but right now it is putting the isotope on the wrong atom
  Molecule mcopy(m);
  const int n = presult.number_elements();

  for (int i = 0; i < n; ++i)
  {
    const atom_number_t j = presult[i];

    mcopy.set_isotope(j, i + 1);
  }

  for (int i = 0; i < matoms; ++i)
  {
    const Atom * a = mcopy.atomi(i);

    if (a->isotope())
      continue;

    Atom * c = new Atom(6);
    mcopy.add(c);
    mcopy.add_bond(presult[i], mcopy.natoms() - 1, SINGLE_BOND);
  }

  Set_of_Atoms to_remove;

  for (int i = matoms - 1; i >= 0; --i)   // note that we use the old natoms value!
  {
    if (0 == mcopy.isotope(i))
      to_remove.add(i);
  }

  mcopy.remove_atoms(to_remove);

  Molecule xx(m);
  cerr << "from " << xx.smiles() << " form " << mcopy.smiles() << endl;
#endif

  int * data = new_int(matoms); std::unique_ptr<int[]> free_data(data);
  presult.set_vector(data, 1);

  int * atomInQuery = new_int(matoms, -1); std::unique_ptr<int[]> free_atomInQuery(atomInQuery);

  for (std::size_t j = 0; j < presult.size(); ++j) {
    atomInQuery[presult[j]] = j;
  }

  int * aam = new_int(matoms); std::unique_ptr<int[]> free_aam(aam);    // index of M atom in FRAG

  Molecule frag;
  m.create_subset(frag, data, 1, aam);

  const int natoms_in_frag = frag.natoms(); 

  int * atomInPattern = new_int(natoms_in_frag, -1); std::unique_ptr<int[]> free_atomInPattern(atomInPattern);   // index of atom in FRAG in M

  for (auto j = 0; j < matoms; ++j) {
    if (aam[j] != -1) {
      atomInPattern[aam[j]] = j;
    }
  }

  for (int j = 0; j < natoms_in_frag; ++j)
  {
    Atom* atom_j = new Atom("C");

    assert(atomInPattern[j] >= 0);
    assert(atomInQuery[atomInPattern[j]] >= 0);

    atom_j->set_isotope(atomInQuery[atomInPattern[j]] + 1);
    frag.add(atom_j);
    frag.add_bond(j, natoms_in_frag + j, SINGLE_BOND);
  }

  Molecule mm(m);
  cerr << "From " << mm.smiles() << " generate fragment " << frag.smiles() << endl;

  IWString ss = frag.unique_smiles();
  gp_group_r->push_scaffold(query_number, ss.c_str());

  return 1;
}


#ifdef NOT_USED_ERQWER
static int
gp_add_scaffold (const Molecule & m,
                 int query_number,
                 const Set_of_Atoms * presult)
{
  if (0 == gp_group_r)
    return 0;

  Molecule mm(m);
  std::vector<int> data(mm.natoms());
  presult->set_vector(&data[0], 1);    // likely wrong, possible uninitialised values

  std::map<int, int> atomInQuery;
  for (std::size_t j = 0; j < presult->size(); ++j) {
    atomInQuery[presult->item(j)] = j;
  }

  const int matoms = mm.natoms();
  int * aam = new_int(matoms); std::unique_ptr<int[]> free_aam(aam);
//    std::vector<int> aam(mm.natoms(), 0);
  Molecule frag;
//    mm.create_subset(frag, &data[0], 1, &aam[0]);
  mm.create_subset(frag, &data[0], 1, aam);
  int natoms_in_frag = frag.natoms(); 

  std::map<int, int> atomInPattern;
  for (int j = 0; j < matoms; ++j) {
    if (aam[j] != -1) {
      atomInPattern[aam[j]] = j;
    }
  }

  for (int j = 0; j < natoms_in_frag; ++j) {
    Atom* atom_j = new Atom("C");

    assert(atomInPattern.find(j) != atomInPattern.end());
    assert(atomInQuery.find(atomInPattern[j]) != atomInQuery.end());

    atom_j->set_isotope(atomInQuery[atomInPattern[j]] + 1);
    frag.add(atom_j);
    frag.add_bond(j, natoms_in_frag + j, SINGLE_BOND);
  }

  cerr << "From " << mm.smiles() << " generate fragment " << frag.smiles() << endl;

  IWString ss = frag.unique_smiles();
  gp_group_r->push_scaffold(query_number, ss.c_str());

  return 1;
}
#endif

static int
identify_substituents (Molecule & m,
                       const int query_number,
                       int & fragments_this_molecule,
                       const int nhits,
                       const Substructure_Results & sresults,
                       IWString & data_for_jibo,
                       int & parent_molecule_written)
{
  if (verbose > 1)
    cerr << "identify_substituents:processing " << nhits << " embeddings\n";

  int * in_fragment = new int[m.natoms()]; std::unique_ptr<int[]> free_in_fragment(in_fragment);

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    if (nullptr != gp_group_r)
      gp_add_scaffold_iw_v2(m,query_number,*e);

    identify_substituents(m, query_number, *e, in_fragment, data_for_jibo, parent_molecule_written);
  }

  return 1;
}

static int
count_unmatched_connections (const Atom * a,
                             const atom_number_t zatom,
                             const Set_of_Atoms & embedding)
{
  int rc = 0;

  const int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (! embedding.contains(j))
      rc++;
  }

  return rc;
}

static int
count_substitutions (Molecule & m,
           const Set_of_Atoms & embedding,
           const Matched_Atoms_From_Which_to_Grow & attachment_points)
{
  int rc = 0;

  int n = attachment_points.number_elements();

  for (int i = 0; i < n; i++)
  {
    int j = attachment_points[i];

    atom_number_t k = embedding[j];

    rc += count_unmatched_connections(m.atomi(k), k, embedding);
  }

  return rc;
}

static int
count_substitutions (Molecule & m,
                     const Set_of_Atoms & embedding)
{
  int ns = embedding.number_elements();

  int rc = 0;
  for (int i = 0; i < ns; i++)
  {
    atom_number_t ai = embedding[i];

    rc += count_unmatched_connections(m.atomi(ai), ai, embedding);
  }

  return rc;
}

/*
*/

static int
reduce_hits_to_those_with_most_substituents(Molecule & m, 
                                            Substructure_Results & sresults,
                                             const Matched_Atoms_From_Which_to_Grow & attachment_points)
{
  int n = sresults.number_embeddings();

  int * number_substituents = new_int(n); std::unique_ptr<int[]> free_number_substituents(number_substituents);

  for (int i = 0; i < n; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    if (attachment_points.number_elements())
      number_substituents[i] = count_substitutions(m, *e, attachment_points);
    else
      number_substituents[i] = count_substitutions(m, *e);
  }

  const auto largest = std::max_element(number_substituents, number_substituents + n) - number_substituents;

  for (int i = n - 1; i >= 0; i--)
  {
    if (i != largest)
      sresults.remove_embedding(i);
  }

  assert(1 == sresults.number_embeddings());

  return 1;
}

static int
count_substitutions(Molecule & m,
                    const Substructure_Results & sresults)
{
  int nhits = sresults.number_embeddings();

#ifdef ONLY_HANDLE_NHITS_BEING_1
  if (1 == nhits)
    ;
  else if (take_first_of_multiple_hits)
    nhits = 1;
  else
  {
    cerr << "Don't know how to handle nhits = " << nhits << ", '" << m.name() << "', taking first hit\n";
    nhits = 1;
  }
#endif

  int rc = 0;

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    rc += count_substitutions(m, *e);
  }

  return rc;
}

static int
substitutions (Molecule & m,
               int * zresults,
               IWString & data_for_jibo)
{
  int parent_molecule_written = 0;

  Molecule_to_Match target(&m);

  int fragments_this_molecule = 0;

  int queries_matching_this_molecule = 0;

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (verbose > 2)
      cerr << nhits << " hits to query " << i << endl;
          
    if(remove_non_single_bond_match || remove_ring_bond_match )
      remove_embeddings(m,sresults);
  
    if(! process_all_matching_queries)
      remove_embeddings_that_result_in_duplicate_atoms(matched_atoms_from_which_to_grow[i], i, sresults);

    nhits = sresults.number_embeddings();
     //cerr << nhits << " hits to query after remove " << i << endl;

    if (1 == nhits)    // hopefully the most common case
      ;
    else if (0 == nhits)
    {
      if (ignore_queries_not_matching)
        continue;

      cerr << "Query " << i << " '" << queries[i]->comment() << "' no hits to '" << m.name() << "'\n";
      return 0;
    }
    else if (ignore_queries_hitting_multiple_times)
    {
      cerr << "Skipped because of " << ignore_queries_hitting_multiple_times << endl;
      continue;
    }
    else if (take_first_of_multiple_hits)
      nhits = reduce_hits_to_those_with_most_substituents(m, sresults, matched_atoms_from_which_to_grow[i]);
    else if (process_all_matching_queries)   // do them all!
      ;
    else
    {
      cerr << "Don't know how to handle nhits " << nhits << " hits to query '" << m.name() << "'\n";
      return 0;
    }

    if (verbose > 3)
      cerr << " nhits now " << nhits << endl;


    if (stream_for_substituents.active() || stream_for_jw.active() || gp_group_r || stream_for_tb.is_open())
      identify_substituents(m, i, fragments_this_molecule, nhits, sresults, data_for_jibo, parent_molecule_written);

    if (write_to_stdout)
      zresults[i] = count_substitutions(m, sresults);

    queries_matching_this_molecule++;

    if (break_after_first_query_match)
      break;
  }

  fragments_per_molecule[fragments_this_molecule]++;

  if (0 == queries_matching_this_molecule)
  {
    if (! ok_for_no_queries_to_match || verbose)
      cerr << "None of " << nq << " queries matched '" << m.name() << "'\n";

    molecules_not_hitting_any_queries++;
    if (stream_for_molecules_not_matching.active())
      stream_for_molecules_not_matching.write(m);

    return ok_for_no_queries_to_match;
  }

  return 1;
}

static int
substitutions(Molecule & m,
              std::ostream & output)
{
  int * zresults = new_int(nq); std::unique_ptr<int[]> free_zresults(zresults);

  IWString data_for_jibo;

  if (gp_group_r) {
    IWString ss1 = m.smiles(), ss2 = m.name();
    gp_group_r->push_molecule(ss1.c_str(), ss2.c_str());
  }

  int rc = substitutions(m, zresults, data_for_jibo);

  if (! stream_for_jw.is_open())   // nothing to write
    ;
  else if (data_for_jibo.length())  // information added
  {
    stream_for_jw << m.smiles() << ' ' << m.name(); // << '\n';
    stream_for_jw << data_for_jibo;
    stream_for_jw << '\n';
    stream_for_jw.write_if_buffer_holds_more_than(32768);
  }
  else if (write_no_match_records)
  {
    data_for_jibo << m.smiles() << ' ' << m.name() << " NO_MATCH\n";
    stream_for_jw << data_for_jibo;
  }

  if (! write_to_stdout)
    return rc;

  write_space_suppressed_string(m.name(), output);

  for (int i = 0; i < nq; i++)
  {
    output << ' ' << zresults[i];
  }

  output << '\n';

  return rc;
}

static int
preprocess_molecule(Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process (m);
  if (make_implicit_hydrogens_explicit)  
    m.make_implicit_hydrogens_explicit();
  return 1;
}

static int
substitutions(data_source_and_type<Molecule> & input,
         std::ostream & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    if (verbose > 1)
      cerr << molecules_read << " '" << m->name() << "'\n";

    (void) preprocess_molecule(*m);

    if (0 == m->natoms())
    {
      cerr << "Ignoring empty molecule '" << m->name() << "'\n";
      continue;
    }

    if (! substitutions(*m, output))
      return 0;
  }

  return 1;
}

static int
substitutions(const char * fname,
         FileType input_type,
         std::ostream & output)
{
  if (FILE_TYPE_INVALID == input_type)
    input_type = discern_file_type_from_name(fname);

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }
  
  return substitutions(input, output);
}

static int
add_any_no_matched_atoms_between_directives (Substructure_Hit_Statistics & q,
                       const Matched_Atoms_From_Which_to_Grow & matched_atoms_from_which_to_grow)
{
  int n = matched_atoms_from_which_to_grow.number_join_points();

  for (int i = 0; i < n; i++)
  {
    if (matched_atoms_from_which_to_grow.is_single_attachment_point(i))
      continue;

    int a1 = matched_atoms_from_which_to_grow.a1(i);
    int a2 = matched_atoms_from_which_to_grow.a2(i);

    if (! q[0]->add_no_matched_atoms_between_initial_atom_numbers(a1, a2))
    {
      cerr << "Huh, cannot add no matched atoms directive between atoms " << a1 << " and " << a2 << endl;
      return 0;
    }
  }

  return 1;
}

static int
add_any_no_matched_atoms_between_directives (resizable_array_p<Substructure_Hit_Statistics> & queries,
                       const Matched_Atoms_From_Which_to_Grow * matched_atoms_from_which_to_grow)
{
  int n = queries.number_elements();

  for (int i = 0; i < n; i++)
  {
    if (! add_any_no_matched_atoms_between_directives(*(queries[i]), matched_atoms_from_which_to_grow[i]))
      return 0;
  }

  return 1;
}

static int
process_isotopically_labelled_query_molecule (MDL_Molecule & m,
                        const Set_of_Atoms & iso,
                        resizable_array_p<Substructure_Hit_Statistics> & queries)
{
  matched_atoms_from_which_to_grow[queries.number_elements()].add(iso);

  m.transform_to_non_isotopic_form();

  Molecule_to_Query_Specifications mqs;

  Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics();

  if (! q->create_from_molecule(m, mqs))
  {
    cerr << "Cannot create query from molecule\n";
    return 0;
  }

  queries.add(q);

  return 1;
}

static int
create_from_smarts (const_IWSubstring & smarts,
          resizable_array_p<Substructure_Hit_Statistics> & queries)
{
  Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics;

  if (! q->create_from_smarts(smarts))
  {
    delete q;

    cerr << "Cannot parse smarts '" << smarts << "'\n";
    return 0;
  }

  queries.add(q);

  return 1;
}

static int
read_isotopically_labelled_query_molecule (MDL_Molecule & m,
                       resizable_array_p<Substructure_Hit_Statistics> & queries)
{
  int ndx = queries.number_elements();   // index into matched_atoms_from_which_to_grow array

  Set_of_Atoms iso;
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    if (m.isotope(i) > 0)
      iso.add(i);
  }

  if (iso.empty())
  {
    cerr << "No isotopic atoms in '" << m.name() << "'\n";
    return 0;
  }

  if (multiple_isotopic_attachment_points)
    return process_isotopically_labelled_query_molecule(m, iso, queries);

  if (1 == iso.number_elements())
    matched_atoms_from_which_to_grow[ndx].set_single_attachment_point(iso[0]);
  else if (2 == iso.number_elements())
    matched_atoms_from_which_to_grow[ndx].set_double_attachment_point(iso[0], iso[1]);
  else
  {
    cerr << "Sorry, don't know how to handle " << iso.number_elements() << " substitution points\n";
    return 0;
  }

  Molecule_to_Query_Specifications mqs;

  mqs.set_substituents_only_at_isotopic_atoms(1);

  //m.build(m);

  //m.only_allow_substitutions_at_isotopic_atoms(mqs);

  Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics();

  if (! q->create_from_molecule(m, mqs))
  {
    cerr << "Cannot create query from molecule\n";
    return 0;
  }

  queries.add(q);

  //q->write_msi(cerr);

  return queries.number_elements();
}

static int
read_isotopically_labelled_query_molecules (data_source_and_type<MDL_Molecule> & input,
                      resizable_array_p<Substructure_Hit_Statistics> & queries)
{
  MDL_Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<MDL_Molecule> free_m(m);

    if (! read_isotopically_labelled_query_molecule(*m, queries))
    {
      cerr << "Fatal error processing query molecule '" << m->name() << "'\n";
      return 0;
    }
  }

  return queries.number_elements();
}

static int
read_isotopically_labelled_query_molecules(const const_IWSubstring & fname,
                      resizable_array_p<Substructure_Hit_Statistics> & queries)
{
  assert(queries.empty());

  const FileType input_type = discern_file_type_from_name(fname);
  assert (FILE_TYPE_INVALID != input_type);

  data_source_and_type<MDL_Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << "cannot open '" << fname << "'\n";
    return 0;
  }

  nq = input.molecules_remaining();      // this logic breaks if multiple files specified..
  if (0 == nq)
  {
    cerr << "Empty file\n";
    return 0;
  }

  matched_atoms_from_which_to_grow = new Matched_Atoms_From_Which_to_Grow[nq];

  return read_isotopically_labelled_query_molecules(input, queries);
}

static void
build_header (IWString & header)
{
  for (int i = 0; i < nq; i++)
  {
    if (i > 0)
      header += ' ';

    header += "SUBS";

    header << i;
  }

  return;
}

static int
substitutions (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "e:wBRvi:E:A:g:lq:s:z:nS:o:j:kC:prxbyJ:d:c:um:Z:hP:T:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (! process_elements(cl))
  {
    usage (2);
  }

  if (! process_standard_aromaticity_options(cl, verbose))
  {
    cerr << "Cannot process aromaticity options (-A)\n";
    usage (5);
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }


  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;
    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }
  if (cl.option_present('w'))
  {
    add_properties_to_groups = 1;
    set_default_iwstring_float_concatenation_precision(4);
    if (verbose)
      cerr << "Will add properties to R groups of -d option\n";
  }
  if (cl.option_present('B'))  
  {
    remove_non_single_bond_match =1;
  }
  if (cl.option_present('R'))  
  {
    remove_ring_bond_match =1;
  }
  if (cl.option_present('b'))
  {
    break_after_first_query_match = 1;
    if (verbose)
      cerr << "Will break after first query match\n";
  }

  if (cl.option_present('n'))
  {
    write_to_stdout = 0;

    if (verbose)
      cerr << "Normal output suppressed\n";
  }

  if (cl.option_present('y'))
  {
    write_parent_molecule = 0;

    if (verbose)
      cerr << "Will NOT write the parent molecule to the output stream\n";
  }

  if (cl.option_present('z'))
  {
    int i = 0;
    const_IWSubstring m;
    while (cl.value('z', m, i++))
    {
      if ("first" == m || 'f' == m)
      {
        take_first_of_multiple_hits = 1;
        if (verbose)
          cerr << "Will take the first of multiple hits for user query\n";
      }
      else if ("ignore" == m || "i" == m)
      {
        ignore_queries_not_matching = 1;
        if (verbose)
          cerr << "Will ignore queries not matching\n";
      }
      else if ("nom" == m || "ignmmatch" == m)
      {
        ignore_queries_hitting_multiple_times = 1;
        if (verbose)
          cerr << "Will ignore multiple query matches in the user query\n";
      }
      else if ("oknomatch" == m)
      {
        ok_for_no_queries_to_match = 1;
        if (verbose)
          cerr << "It is OK if none of the queries match\n";
      }
      else if ("each" == m)
      {
        process_all_matching_queries = 1;
        if (verbose)
          cerr << "Will process all matching queries\n";
      }
      else
      {
        cerr << "Unrecognised -z qualifier '" << m << "'\n";
        usage(17);
      }
    }

    if (process_all_matching_queries && take_first_of_multiple_hits)
    {
      cerr << "Specifying '-z first' and '-z each' is incompatible\n";
      usage(4);
    }
  }
  
  if (cl.option_present('e'))  
  {
    int i = 0;
    const_IWSubstring smiles;
    while (cl.value('e', smiles, i++))
    {
      Molecule m;
      if (! m.build_from_smiles(smiles))
      {
        cerr << "build_query_from_smiles:invalid smiles '" << smiles << "'\n";
        usage(1);
      }

      Molecule_to_Query_Specifications mqs;
      mqs.set_make_embedding(1);
      mqs.set_condense_explicit_hydrogens_to_anchor_atoms(1);

      Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics();

      if (! q->create_from_molecule(m, mqs))
      {
        cerr << "build_query_from_smiles:invalid molecule?? '" << smiles << "'\n";
        delete q;
        usage(1);
      }

      queries.add(q);
    }
  }

  if (cl.option_present('s'))
  {
    int i = 0;
    const_IWSubstring smarts;
    while (cl.value('s', smarts, i++))
    {
      if (! create_from_smarts(smarts, queries))
      {
        cerr << "Cannot parse -s option '" << smarts << "'\n";
        return 1;
      }
    }
  }

  if (cl.option_present('q'))
  {
    int i = 0;
    const_IWSubstring q;

    while (cl.value('q', q, i++))
    {
      if (! process_cmdline_token('q', q, queries, verbose))
      {
        cerr << "Cannot parse -q option '" << q << "'\n";
        return 4;
      }
    }
  }

  if (cl.option_present('m'))
  { 
    if (cl.option_present('j'))
    {
      cerr << "Cannot use the -j option with the -m option - connection points determined automatically\n";
      return 5;
    }

    if (cl.option_present('h'))
    {
      multiple_isotopic_attachment_points = 1;

      if (verbose)
        cerr << "Isotopic labels designate attachment points - any number\n";
    }

    int i = 0;
    const_IWSubstring m;
    while (cl.value('m', m, i++))
    {
      if (! read_isotopically_labelled_query_molecules(m, queries))
      {
        cerr << "Cannot read query molecules '" << m << "'\n";
        return 4;
      }
    }
  }

  nq = queries.number_elements();

  if (verbose)
    cerr << "Defined " << nq << " user defined queries\n";

  if (0 == nq)
  {
    cerr << "No queries defined, cannot continue\n";
    usage(4);
  }

  // If there is just one query, and they've entered '-z i' that also means that it
  // is OK to ignore the case of all queries not matching

  if (1 == nq && ignore_queries_not_matching)
    ok_for_no_queries_to_match = 1;

  if (cl.option_present('k'))
  {
    for (int i = 0; i < nq; i++)
    {
      queries[i]->set_do_not_perceive_symmetry_equivalent_matches(1);
    }
  }

  if (cl.option_present('r'))
  {
    treat_rings_as_substitutents = 1;

    if (verbose)
      cerr << "Will treat rings as substituents\n";
  }

  if (cl.option_present('P'))
  {
    const char * fname = cl.option_value('P');

    if (! read_capping_group(capping_group, fname))
    {
      cerr << "Cannot read capping group '" << fname << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Fragments capped with " << capping_group.unique_smiles() << "\n";

    adjust_hydrogen_data(capping_group);
  }

  if (cl.option_present('c'))
  {
    if (! cl.value('c', isotope_for_ring_closures) || isotope_for_ring_closures < 0)
    {
      cerr << "The isotope for ring closing substituents value (-c) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will label atoms in substituents that close rings with " << isotope_for_ring_closures << endl;
  }

  if (cl.option_present('d')) {
    const char * j = cl.option_value('d');
    gs_rgroup_file = j;
    gp_group_r = new GroupR();
    set_display_abnormal_valence_messages(0);
  }

  if (cl.option_present('u'))
  {
    write_no_match_records = 0;
    if (verbose)
      cerr << "Will NOT write NO_MATCH data\n";
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  // If we are reading isotopically labelled molecules, the
  // matched_atoms_from_which_to_grow array will already have been set up

  if (nullptr == matched_atoms_from_which_to_grow)
    matched_atoms_from_which_to_grow = new Matched_Atoms_From_Which_to_Grow[nq];

  // We can optionally keep track of the substituents found

  if (cl.option_present('S'))
  {
    if (cl.option_present('o'))
    {
      if (! stream_for_substituents.determine_output_types (cl, 'o'))
      {
        cerr << "Cannot determine output type(s)\n";
        return 8;
      }
    }
    else
      stream_for_substituents.add_output_type(FILE_TYPE_SMI);

    const_IWSubstring s = cl.string_value('S');

    if (stream_for_substituents.would_overwrite_input_files(cl, s))
    {
      cerr << "Cannot overwrite input files '" << s << "' is invalid\n";
      return 4;
    }

    if (! stream_for_substituents.new_stem(s))
    {
      cerr << "Cannot create output file(s) with stem '" << s << "'\n";
      return 8;
    }

    if (verbose)
      cerr << "Substituents written to '" << s << "'\n";
  }

  if (cl.option_present('Z'))
  {
    if (cl.option_present('o'))
    {
      if (! stream_for_molecules_not_matching.determine_output_types(cl, 'o'))
      {
        cerr << "Cannot determine output type(s)\n";
        return 8;
      }
    }
    else
      stream_for_molecules_not_matching.add_output_type(FILE_TYPE_SMI);

    const_IWSubstring s = cl.string_value('Z');

    if (stream_for_molecules_not_matching.would_overwrite_input_files(cl, s))
    {
      cerr << "Cannot overwrite input files '" << s << "' is invalid\n";
      return 4;
    }

    if (! stream_for_molecules_not_matching.new_stem(s))
    {
      cerr << "Cannot create output file(s) with stem '" << s << "'\n";
      return 8;
    }

    if (verbose)
      cerr << "Non matching molecules written to written to '" << s << "'\n";
  }

  if (cl.option_present('S') || cl.option_present('J') || cl.option_present('d') || cl.option_present('T'))
  {
    if (cl.option_present('j'))
    {
      if (cl.option_count('j') > nq)
      {
        cerr << "Can be just one -j option for each query, nq = " << nq << endl;
        return 4;
      }

      int i = 0;    // index into cl
      const_IWSubstring s;
      int j = 0;    // index into matched_atoms_from_which_to_grow array
      while (cl.value('j', s, i++))
      {
        if (! matched_atoms_from_which_to_grow[j].build(s))
        {
          cerr << "Invalid embedding member(s) specification '" << s << "'\n";
          return 5;
        }
        j++;
      }
    }

    if (cl.option_present('j'))
    {
      if (! add_any_no_matched_atoms_between_directives(queries, matched_atoms_from_which_to_grow))
      {
        cerr << "Fatal error setting up no matched atoms between directives\n";
        return 5;
      }
    }

    if (cl.option_present('C'))
    {
      if (! cl.value ('C', max_atoms_in_a_fragment) || max_atoms_in_a_fragment < 1)
      {
        cerr << "Invalid value for the max atoms in a fragment option (-C)\n";
        usage(5);
      }

      if (verbose)
        cerr << "Will not write fragments with more than " << max_atoms_in_a_fragment << " atoms\n";
    }

    if (cl.option_present('p'))
    {
      compute_unique_smiles_of_fragments = 1;

      if (verbose)
        cerr << "Will compute the unique smiles of the substituents\n";
    }

    if (cl.option_present('x'))
    {
      include_anchor_atom_with_substituent = 0;

      if (verbose)
        cerr << "Will exclude the scaffold atom from the substituent fragment\n";
    }
  }
  else if (cl.option_present('p') || cl.option_present('C'))
  {
    cerr << "The -p and -C options only make sense with the -S option\n";
    usage(6);
  }

  for (int i = 0; i < nq; i++)
  {
    queries[i]->set_respect_initial_atom_numbering(1);
  }

  if (cl.option_present('J'))
  {
    const char * j = cl.option_value('J');

    if (! stream_for_jw.open(j))
    {
      cerr << "Cannot open tabular file '" << j << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Tabular output written to '" << j << "'\n";
  }

  if (cl.option_present('T'))
  {
    IWString t = cl.string_value('T');

    if (! t.ends_with(".smi"))
      t << ".smi";

    if (! stream_for_tb.open(t.null_terminated_chars()))
    {
      cerr << "Cannot open Thibault output file '" << t << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Tibault output written to '" << t << "'\n";
  }

  if (cl.empty())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  // Write the header for descriptor files

  IWString header;

  build_header(header);

  if (write_to_stdout)
    std::cout << "Name " << header << '\n';

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! substitutions(cl[i], input_type, std::cout))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";

    for (int i = 0; i < matches_per_molecule.number_elements(); i++)
    {
      if (matches_per_molecule[i])
        cerr << matches_per_molecule[i] << " molecules matched " << i << " queries\n";
    }

    for (int i = 0; i < nq; i++)
    {
      queries[i]->report(cerr, verbose);
    }

    if (molecules_not_hitting_any_queries)
      cerr << molecules_not_hitting_any_queries << " hit none of the queries\n";

    if (fragments_per_molecule.number_elements())
    {
      for (int i = 0; i < fragments_per_molecule.number_elements(); i++)
      {
        if (fragments_per_molecule[i])
          cerr << fragments_per_molecule[i] << " molecules had " << i << " fragments\n";
      }
    }

    if (fragments_with_too_many_atoms)
      cerr << fragments_with_too_many_atoms << " with more than " << max_atoms_in_a_fragment << " atoms suppressed\n";
  }

  if (fragments_found_at_substitution_point.size() > 0)
  {
    for (int i = 0; i < nq; i++)
    {
      IWString prefix;
      prefix << i << '.';

      cerr << "Attachment summary query " << i;
      if (queries[i]->comment().length())
        cerr << " '" << queries[i]->comment().length() << "':";
      else
        cerr << ':';

      for (IW_STL_Hash_Map<IWString, Fragments_Found *>::const_iterator i = fragments_found_at_substitution_point.begin(); i != fragments_found_at_substitution_point.end(); ++i)
      {
        const IWString & s = (*i).first;
        if (s.starts_with(prefix))
          cerr << ' ' << (*i).first << " (" << (*i).second->size () << ")";
      }
      cerr << endl;
    }

    for (IW_STL_Hash_Map<IWString, Fragments_Found *>::const_iterator i = fragments_found_at_substitution_point.begin(); i != fragments_found_at_substitution_point.end(); ++i)
    {
      cerr << "Attachments at '" << (*i).first << "', N = " << (*i).second->size() << endl;
      if(verbose)
        (*i).second->report((*i).first, cerr);
    }
  }

  if (gp_group_r) {
    gp_group_r->close();
    std::ofstream ss(gs_rgroup_file.c_str());
    ss << *gp_group_r;
    ss.close();
    delete gp_group_r;
  }

  if(stream_for_jw.is_open())
    stream_for_jw.close();

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = substitutions(argc, argv);

  return rc;
}
