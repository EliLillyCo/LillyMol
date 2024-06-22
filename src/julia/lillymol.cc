#include <ranges>
#include <string>
#include <tuple>
#include <vector>

#include "Eigen/Core"

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/tuple.hpp"
//#include "jlcxx/stl.hpp"

#include "Molecule_Lib/chiral_centre.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/ring_data.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/xlogp.h"
#include "Molecule_Tools/alogp.h"

#include "julia/resizable_array_holder.h"

namespace lillymol_julia {

std::string
greet() {
  return "hello world";
}

struct World
{
  World(const std::string& message = "default hello") : msg(message){}
  void set(const std::string& msg) { this->msg = msg; }
  std::string greet() { return msg; }
  std::string msg;
  ~World() { std::cout << "Destroying World with message " << msg << std::endl; }
};

JLCXX_MODULE define_types_module(jlcxx::Module& types)
{
  types.add_type<World>("World")
    .constructor<const std::string&>()
    .method("set", &World::set)
    .method("greet", &World::greet)
  ;
}

enum BondType {
  kInvalidBond = 0,
  kSingleBond = 1,
  kDoubleBond = 2,
  kTripleBond = 3,
  kAromaticBond = 4
};

BondType
ToBondType(const Bond& b) {
  if (b.is_aromatic()) {
    return kAromaticBond;
  }
  if (b.is_single_bond()) {
    return kSingleBond;
  }
  if (b.is_double_bond()) {
    return kDoubleBond;
  }
  if (b.is_triple_bond()) {
    return kTripleBond;
  }

  return kInvalidBond;
}

bond_type_t
BtypeEnumToBtype(const BondType btenum) {
  switch (btenum) {
    case kSingleBond:
      return SINGLE_BOND;
    case kDoubleBond:
      return DOUBLE_BOND;
    case kTripleBond:
      return TRIPLE_BOND;
    case kAromaticBond:
      return AROMATIC_BOND;
    default:
      return NOT_A_BOND;
  }
}

class SetOfRings : public ResizableArrayHolder<Ring> {
  private:
  public:
    SetOfRings(const resizable_array_p<Ring>& r) : ResizableArrayHolder<Ring>(r) {
    }
    const Ring* operator[](int ndx) const {
      return _ref[ndx];
    }
};

// Make it convenient to extract all the ring atoms from a molecule.
class RingAtoms {
  private:
    int _nrings;
    Set_of_Atoms* _rings;
  public:
    RingAtoms();
    ~RingAtoms();

    int GatherRings(Molecule& m);

    int nrings() const {
      return _nrings;
    }

    const Set_of_Atoms& operator[](int ndx) const {
      return _rings[ndx];
    }
};

RingAtoms::RingAtoms() {
  _nrings = 0;
  _rings = nullptr;
}

RingAtoms::~RingAtoms() {
  _nrings = -1;
  delete [] _rings;
}

int
RingAtoms::GatherRings(Molecule& m) {
  if (_rings) {
    delete [] _rings;
    _rings = nullptr;
  }

  _nrings = m.nrings();
  if (_nrings == 0) {
    return 1;
  }
  _rings = new Set_of_Atoms[_nrings];
  for (int i = 0; i < _nrings; ++i) {
    _rings[i] = *m.ringi(i);
  }

  return _nrings;
}


template <typename T>
IWString
AtomsAsString(const T & s, const char* name) {
  IWString result;
  result << name << " : N=" << s.size() << " [";
  bool need_space = false;
  for (atom_number_t a : s) {
    if (need_space) {
      result << ' ';
    } else {
      need_space = true;
    }
    result << a;
  }
  result << ']';
  return result;
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{

  mod.method("greet", &greet);

  mod.add_type<World>("World")
    .constructor<const std::string&>()
    .method("set", &World::set)
    .method("greet", &World::greet)
  ;

  mod.add_bits<BondType>("BondType", jlcxx::julia_type("CppEnum"));
  mod.set_const("SINGLE_BOND", kSingleBond);
  mod.set_const("DOUBLE_BOND", kDoubleBond);
  mod.set_const("TRIPLE_BOND", kTripleBond);
  mod.set_const("AROMATIC_BOND", kAromaticBond);

  mod.add_bits<FileType>("FileType", jlcxx::julia_type("CppEnum"));
  mod.set_const("SMI", FILE_TYPE_SMI);
  mod.set_const("SDF", FILE_TYPE_SDF);

  mod.add_bits<lillymol::RIProperty>("RIProperty", jlcxx::julia_type("CppEnum"));
  mod.set_const("RIP_NONE", lillymol::RIProperty::kNone);
  mod.set_const("RIP_AROMATIC", lillymol::RIProperty::kAromatic);
  mod.set_const("RIP_FUSED", lillymol::RIProperty::kFused);

  mod.add_type<lillymol::RingInformation>("RingInformation")
    .constructor<>()
    .method("nrings", &lillymol::RingInformation::nrings)
    .method("aromatic", &lillymol::RingInformation::aromatic)
    .method("fused", &lillymol::RingInformation::fused)
  ;
    
//  mod.add_bits<FileType>("FileType", jlcxx::julia_type("CppEnum"));
//  mod.set_const("SMI", FILE_TYPE_SMI);
//  mod.set_const("SDF", FILE_TYPE_SDF);

  mod.add_type<data_source_and_type<Molecule>>("MoleculeReader")
    .constructor<FileType, std::string&>()
    .method("next_molecule",
      [](data_source_and_type<Molecule>& input, jlcxx::BoxedValue<Molecule>& boxed_molecule)->bool{
        Molecule& m = jlcxx::unbox<Molecule&>(boxed_molecule);
        if (input.next_molecule(m)) {
          return true;
        }

        return false;
      }
    )
    .method("next_molecule",
      [](jlcxx::BoxedValue<data_source_and_type<Molecule>>& boxed_input, jlcxx::BoxedValue<Molecule>& boxed_molecule)->bool{
        data_source_and_type<Molecule>& input = jlcxx::unbox<data_source_and_type<Molecule>&>(boxed_input);
        Molecule& m = jlcxx::unbox<Molecule&>(boxed_molecule);
        return input.next_molecule(m);
      }
    )
    .method("set_connection_table_errors_allowed",
      [](data_source_and_type<Molecule>& input) {
        input.set_connection_table_errors_allowed(std::numeric_limits<int>::max());
      },
      "Allow unlimited connection table errors"
    )
    .method("set_connection_table_errors_allowed",
      [](data_source_and_type<Molecule>& input, int s) {
        input.set_connection_table_errors_allowed(s);
      },
      "Specify the max number of connection table errors to ignore"
    )
    .method("connection_table_errors_encountered",
      [](const data_source_and_type<Molecule>& input) {
        return input.connection_table_errors_encountered();
      },
      "The number of connection table errors encountered"
    )
    .method("set_skip_first",
      [](data_source_and_type<Molecule>& input, int s) {
        input.set_skip_first(s);
      },
      "Skip the first <s> records in the file"
    )
    .method("read_only",
      [](data_source_and_type<Molecule>& input, int s) {
        input.set_do_only(s);
      },
      "Only read the first <s> records"
    )
    .method("molecules_remaining",
      [](data_source_and_type<Molecule>& input) {
        return input.molecules_remaining();
      },
      "The number of molecules not yet read - does not work with pipes"
    )

    .method("molecules_read",
      [](const data_source_and_type<Molecule>& input) {
        return input.molecules_read();
      },
      "The number of molecules read"
    )
  ;

  mod.add_type<Chiral_Centre>("ChiralCentre")
    .constructor<atom_number_t>()
    .method("centre",
      [](const Chiral_Centre& c)->atom_number_t{
        return c.centre();
      }
    )
    .method("centre",
      [](const jlcxx::BoxedValue<const Chiral_Centre>& boxed_chiral_centre)->atom_number_t{
      const Chiral_Centre& c = jlcxx::unbox<const Chiral_Centre&>(boxed_chiral_centre);
      return c.centre();
      }
    )
    .method("top_front", &Chiral_Centre::top_front)
    .method("top_front",
      [](const jlcxx::BoxedValue<const Chiral_Centre>& boxed_chiral_centre)->atom_number_t{
      const Chiral_Centre& c = jlcxx::unbox<const Chiral_Centre&>(boxed_chiral_centre);
      return c.top_front();
      }
    )
    .method("top_back", &Chiral_Centre::top_back)
    .method("top_back",
      [](const jlcxx::BoxedValue<const Chiral_Centre>& boxed_chiral_centre)->atom_number_t{
      const Chiral_Centre& c = jlcxx::unbox<const Chiral_Centre&>(boxed_chiral_centre);
      return c.top_back();
      }
    )
    .method("left_down", &Chiral_Centre::left_down)
    .method("left_down",
      [](const jlcxx::BoxedValue<const Chiral_Centre>& boxed_chiral_centre)->atom_number_t{
      const Chiral_Centre& c = jlcxx::unbox<const Chiral_Centre&>(boxed_chiral_centre);
      return c.left_down();
      }
    )
    .method("right_down", &Chiral_Centre::right_down)
    .method("right_down",
      [](const jlcxx::BoxedValue<const Chiral_Centre>& boxed_chiral_centre)->atom_number_t{
      const Chiral_Centre& c = jlcxx::unbox<const Chiral_Centre&>(boxed_chiral_centre);
      return c.right_down();
      }
    )
  ;

  mod.method("is_chiral_implicit_hydrogen",
    [](int c)->bool {
      return IsChiralImplicitHydrogen(c);
    },
    "True if chiral connection is an implicit hydrogen"
  );

  mod.add_type<Set_of_Atoms>("SetOfAtoms")
    .constructor<>()
    // Does not work
    // LoadError: No appropriate factory for type St16initializer_listIiE
    //.constructor<std::initializer_list<atom_number_t>>()
    .constructor([] (jlcxx::ArrayRef<int64_t,1> args) {
        std::vector<int64_t> tmp;
        tmp.reserve(args.size());
        for (auto a : args) {
          tmp.push_back(a);
        }
        return new Set_of_Atoms(tmp);
      }
    )
    .method("push!",
      [](Set_of_Atoms& s, atom_number_t zextra) {
        s << zextra;
      }
    )
    .method("contains",
      [](const Set_of_Atoms& s, atom_number_t a)->bool{
        return s.contains(a);
      }
    )
  ;

  mod.set_override_module(jl_base_module);
  mod.method("getindex",
    [](const Set_of_Atoms& s, int ndx)->atom_number_t{
      return s[ndx - 1];
    }
  );
  mod.method("setindex!",
    [](Set_of_Atoms& s, int32_t value, uint64_t ndx) {
      s[ndx - 1] = value;
    }
  );
  mod.unset_override_module();

  mod.method("soa_getindex",
    [](const Set_of_Atoms& s, int ndx) {
      //std::cerr << "Getting item " << ndx << " from thing with " << s.size() << " atoms\n";
      assert(ndx >= 0 && ndx < static_cast<int>(s.size()));
      return s[ndx];
    },
    "Needed to disambiguate things"
  );

  mod.method("add!",
    [](Set_of_Atoms& s, atom_number_t a) {
      s.add(a);
    }
  );

  mod.method("set_of_atoms_show_text",
    [](const Set_of_Atoms& s)->std::string{
      return AtomsAsString(s, "SetOfAtoms").AsString();
    }
  );
  mod.method("set_of_atoms_show_text",
    [](const Set_of_Atoms* s)->std::string {
      return AtomsAsString(*s, "SetOfAtoms").AsString();
    }
  );
  mod.method("equals",
    [](const Set_of_Atoms& s, const jlcxx::ArrayRef<int64_t> v)->bool {
      // std::cerr << "Set of atoms equals sizes " << s.size() << " and " << v.size() << '\n';
      if (s.size() != v.size()) {
        return false;
      }

      for (uint32_t i = 0; i < s.size(); ++i) {
        // std::cerr << " cmp " << s[i] << " and " << v[i] << '\n';
        if (s[i] != v[i]) {
          return false;
        }
      }

      return true;
    }
  );

  mod.add_type<Ring>("Ring")
    .constructor<>()
    .method("atoms_in_ring",
      [](const Ring& r) {
        return r.number_elements();
      }
    )
    .method("ring_number", &Ring::ring_number)
    .method("atoms_in_ring",
      [](const Ring* r) {
        return r->number_elements();
      }
    )
    .method("fragment_membership", &Ring::fragment_membership)
    .method("fused_system_identifier", &Ring::fused_system_identifier)
    .method("fused_ring_neighbours", &Ring::fused_ring_neighbours)
    .method("fused_neighbour", &Ring::fused_neighbour)
    .method("strongly_fused_ring_neighbours", &Ring::strongly_fused_ring_neighbours)
    .method("contains_bond", &Ring::contains_bond)
    .method("contains_both", &Ring::contains_both)
    .method("is_fused_to",
      [](const Ring* r, const jlcxx::BoxedValue<Ring*>& boxed_ring)->bool{
        const Ring* unboxed_ring = jlcxx::unbox<Ring*>(boxed_ring);
        std::cerr << "unboxed_ring " << unboxed_ring->ring_number() << " size " << unboxed_ring->size() << '\n';
        std::cerr << &unboxed_ring << '\n';
        return r->is_fused_to(unboxed_ring);
      }
    )
    .method("is_fused_to",
      [](const Ring* r, const Ring* other)->bool {
        return r->is_fused_to(other);
      }
    )

    .method("is_fused",
      [](const Ring* r)->bool{
        return r->is_fused();
      }
    )
    .method("is_aromatic",
      [](const Ring* r)->bool{
        return r->is_aromatic();
      }
    )
    .method("largest_number_of_bonds_shared_with_another_ring",
      [](const Ring* r) {
        return r->largest_number_of_bonds_shared_with_another_ring();
      }
    )
    .method("ring_show_text",
      [](const Ring* r)->std::string{
        IWString result = AtomsAsString(*r, "Ring");
        if (r->is_aromatic()) {
          result << " arom";
        } else {
          result << " aliph";
        }
        if (r->is_fused()) {
          result << " fused";
        }
        result << " frag " << r->fragment_membership();
        result << " fsid " << r->fused_system_identifier();
        result << " fused " << r->is_fused();
        return result.AsString();
      }
    )
    .method("ring_show_text",
      [](const jlcxx::BoxedValue<Ring>& boxed_ring)->std::string{
        const Ring& r = jlcxx::unbox<Ring&>(boxed_ring);
        IWString result = AtomsAsString(r, "Ring");
        if (r.is_aromatic()) {
          result << " arom";
        } else {
          result << " aliph";
        }
        if (r.is_fused()) {
          result << " fused";
        }
        result << " frag " << r.fragment_membership();
        result << " fsid " << r.fused_system_identifier();
        result << " fused " << r.is_fused();
        return result.AsString();
      }
    )
  ;
  mod.method("ring_equals_vector",
    [](const Set_of_Atoms* s, const jlcxx::ArrayRef<int64_t> v)->bool {
      if (s->size() != v.size()) {
        return false;
      }

      for (uint32_t i = 0; i < s->size(); ++i) {
        if (s->item(i) != v[i]) {
          return false;
        }
      }

      return true;
    }
  );

  mod.set_override_module(jl_base_module);
  mod.method("getindex",
    [](const Ring& a, int i)->atom_number_t{
      assert(i > 0 && i <= a.number_elements());
      return a[i - 1];
    }
  );
  mod.method("length",
    [](const Set_of_Atoms& s){
      return s.number_elements();
    }
  );
  mod.method("length",
    [](const jlcxx::BoxedValue<Set_of_Atoms>& boxed_set_of_atoms){
      const Set_of_Atoms& s = jlcxx::unbox<Set_of_Atoms&>(boxed_set_of_atoms);
      return s.number_elements();
    }
  );
#ifdef IN_RING__
  mod.method("in",
    [](const jlcxx::BoxedValue<Ring>& boxed_ring, atom_number_t atom){
      const Ring& r = jlcxx::unbox<Ring&>(boxed_ring);
      return r.contains(atom);
    }
  );
  mod.method("in",
    [](const jlcxx::BoxedValue<Set_of_Atoms>& boxed_set_of_atoms, atom_number_t atom){
      const Set_of_Atoms& s = jlcxx::unbox<Set_of_Atoms&>(boxed_set_of_atoms);
      return s.contains(atom);
    }
  );
  mod.method("in",
    [](const Ring& r, atom_number_t a) ->bool {
      return r.contains(a);
    }
  );
  mod.method("in",
    [](const Ring* r, atom_number_t a) ->bool {
      return r->contains(a);
    }
  );
#endif


// collect not working, not sure why...
//mod.method("collect",
//  [](const jlcxx::BoxedValue<Ring>& boxed_ring)->std::vector<atom_number_t>{
//    const Ring& r = jlcxx::unbox<Ring&>(boxed_ring);
//    std::vector<atom_number_t> result;
//    result.reserve(r.size());
//    for (atom_number_t a : r) {
//      result.push_back(a);
//    }
//    return result;
//  }
//);
  mod.unset_override_module();

  mod.add_type<Mol2Graph>("Mol2Graph")
    .method("set_exclude_triple_bonds_from_graph_reduction", &Mol2Graph::set_exclude_triple_bonds_from_graph_reduction)
    .method("set_revert_all_directional_bonds_to_non_directional", &Mol2Graph::set_revert_all_directional_bonds_to_non_directional)
    .method("set_preserve_cc_double_bonds_no_heteroatoms ", &Mol2Graph::set_preserve_cc_double_bonds_no_heteroatoms )
    .method("set_preserve_cc_double_bonds_saturated ", &Mol2Graph::set_preserve_cc_double_bonds_saturated )
    .method("set_append_molecular_formula ", &Mol2Graph::set_append_molecular_formula )
    .method("set_aromatic_distinguishing_formula", &Mol2Graph::set_aromatic_distinguishing_formula)
    .method("set_remove_chiral_centres ", &Mol2Graph::set_remove_chiral_centres )
  ;

  mod.add_type<Chemical_Standardisation>("ChemicalStandardisation")
    .method("activate_all", &Chemical_Standardisation::activate_all)
    .method("activate",
      [](Chemical_Standardisation& standardise, const std::string& directive, bool verbose)->bool {
        IWString tmp(directive);
        return standardise.Activate(tmp, verbose);
      },
      "Activate specific transformation(s)"
    )
    .method("process", 
      [](Chemical_Standardisation& s, jlcxx::BoxedValue<Molecule>& boxed_mol)->int {
        Molecule& m = jlcxx::unbox<Molecule&>(boxed_mol);
        return s.process(m);
      }
    )
  ;
    
  mod.add_type<Bond>("Bond")
    .method("a1",
      [](const Bond& b)->int64_t {
        return b.a1();
      }
    )
    .method("a1",
      [](const jlcxx::BoxedValue<const Bond>& boxed_bond)->int64_t {
        const Bond& b = jlcxx::unbox<const Bond&>(boxed_bond);
        return b.a1();
      }
    )
    .method("a2",
      [](const Bond& b)->int64_t {
        return b.a2();
      }
    )
    .method("a2",
      [](const jlcxx::BoxedValue<const Bond>& boxed_bond)->int64_t {
        const Bond& b = jlcxx::unbox<const Bond&>(boxed_bond);
        return b.a2();
      }
    )
    .method("btype",
      [](const Bond& b)->BondType{
        return ToBondType(b);
      }
    )
    .method("other", &Bond::other)
    .method("is_single_bond",
      [](const Bond& b)->bool{
        return b.is_single_bond();
      }
    )
    .method("is_single_bond",
      [](const jlcxx::BoxedValue<Bond>& boxed_bond)->bool{
        const Bond& b = jlcxx::unbox<const Bond&>(boxed_bond);
        return b.is_single_bond();
      }
    )
    .method("is_double_bond",
      [](const Bond& b)->bool{
        return b.is_double_bond();
      }
    )
    .method("is_double_bond",
      [](const jlcxx::BoxedValue<Bond>& boxed_bond)->bool{
        const Bond& b = jlcxx::unbox<const Bond&>(boxed_bond);
        return b.is_double_bond();
      }
    )
    .method("is_triple_bond",
      [](const Bond& b)->bool{
        return b.is_triple_bond();
      }
    )
    .method("is_triple_bond",
      [](const jlcxx::BoxedValue<Bond>& boxed_bond)->bool{
        const Bond& b = jlcxx::unbox<const Bond&>(boxed_bond);
        return b.is_triple_bond();
      }
    )
    .method("is_aromatic",
      [](const Bond& b)->bool{
        return b.is_aromatic();
      }
    )
    .method("is_aromatic",
      [](const jlcxx::BoxedValue<Bond>& boxed_bond)->bool{
        const Bond& b = jlcxx::unbox<const Bond&>(boxed_bond);
        return b.is_aromatic();
      }
    )
    .method("is_aromatic_bond",
      [](const Bond& b)->bool{
        return b.is_aromatic();
      }
    )
    .method("is_aromatic_bond",
      [](const jlcxx::BoxedValue<Bond>& boxed_bond)->bool{
        const Bond& b = jlcxx::unbox<const Bond&>(boxed_bond);
        return b.is_aromatic();
      }
    )
    .method("nrings",
      [](const Bond& b) {
        return b.nrings();
      }
    )
    .method("nrings",
      [](const Bond* b) {
        return b->nrings();
      }
    )
    .method("bond_number", &Bond::bond_number)
    .method("either_atom_set_in_array", &Bond::either_atom_set_in_array)
    .method("atoms",
      [](const Bond& b)->std::tuple<atom_number_t, atom_number_t>{
        return std::make_tuple(b.a1(), b.a2());
      }
    )
    .method("involves",
      [](const Bond& a, atom_number_t o)->bool{
        return a.involves(o);
      }
    )
    .method("involves",
      [](const Bond& a, atom_number_t o1, atom_number_t o2)->bool{
        return a.involves(o1, o2);
      }
    )
    .method("joins",
      [](const Bond& b1, const Bond& b2)->bool{
        return b1.joins(&b2);
      }
    )
    .method("bond_show_text",
      [](const Bond& b) ->std::string {
        IWString result;
        result << b.a1();
        if (b.is_aromatic()) {
          result << ':';
        } else if (b.is_single_bond()) {
          result << '-';
        } else if (b.is_double_bond()) {
          result << '=';
        } else if (b.is_triple_bond()) {
          result << '#';
        }
        result << b.a2();

        return result.AsString();
      }
    )
  ;

  mod.add_type<Atom>("Atom")
    // .constructor<jlcxx::cxxint_t>(true) // with finalizer
    .method("valence_ok", 
      [](Atom& a)->bool {
        return a.valence_ok();
      },
      "True if the valence is OK"
    )
    // It is unfortunate that the valence_ok method is non const.
    .method("valence_ok", 
      [](const Atom& a)->bool {
        Atom* me = const_cast<Atom*>(&a);
        return me->valence_ok();
      },
      "True if the valence is OK"
    )
    .method("ncon",
      [](const Atom& a)->int64_t {
        return a.ncon();
      }
    )
    .method("isotope", &Atom::isotope)
    .method("lone_pair_count",
      [](Atom&a)->int64_t{
        int lp;
        if (a.lone_pair_count(lp)) {
          return lp;
        }
        return 0;
      }
    )
    .method("nbonds",
      [](const Atom& a)->int64_t {
        return a.nbonds();
      }
    )
    .method("nbonds",
      [](Atom& a)->int64_t {
        return a.nbonds();
      }
    )
    .method("bond_to_atom",
      [](const Atom& a, atom_number_t o)->Bond{
        return *a.bond_to_atom(o);
      }
    )
    .method("is_bonded_to",
      [](const Atom& a, atom_number_t o)->bool{
        return a.is_bonded_to(o);
      }
    )
    .method("is_halogen",
      [](const Atom& a)->bool {
        return a.is_halogen();
      },
      "True if the underlying element is a halogen"
    )
    .method("atomic_number", &Atom::atomic_number)
    .method("atomic_symbol",
      [](const Atom& a)->std::string{
        return a.atomic_symbol().AsString();
      }
    )
    .method("implicit_hydrogens", &Atom::implicit_hydrogens)
    .method("formal_charge", &Atom::formal_charge)
    .method("other", &Atom::other)
    .method("connections",
      [](const Atom& a, atom_number_t me)->Set_of_Atoms{
        return a.connections(me);
      }
    )
    .method("atomic_weight", &Atom::atomic_weight)
    .method("saturated", 
      [](const Atom& a)->bool{
        return a.fully_saturated();
      }
    )
    .method("unsaturation", 
      [](const Atom& a)->int64_t {
        return a.unsaturation();
      },
      "nbonds() - ncon()"
    )
    .method("atom_show_text",
      [](const Atom& a) -> std::string {
        IWString result;
        result << "Atom: " << a.atomic_number() << " ncon " << a.ncon();
        return result.AsString();
      }
    )
  ;

    mod.set_override_module(jl_base_module);
    mod.method("getindex",
      [](const Atom& a, int i)->const Bond{
        assert(i >  0 && i <= a.ncon());
        return *a[i - 1];
      }
    );
    mod.method("getindex",
      [](Atom& a, int i)->const Bond{
        assert(i > 0 && i <= a.ncon());
        return *a[i - 1];
      }
    );
    mod.unset_override_module();
  ;
  mod.method("internal_get_item",
    [](const Atom& a, int64_t ndx)->const Bond* {
      return a[ndx - 1];
    }
  );

  mod.add_type<SetOfRings>("SetOfRings")
    .method("size",
      [](const SetOfRings& r) {
        return r.size();
      }
    )
    .method("rings_in_set",
      [](const SetOfRings&r)->int64_t {
        return r.size();
      }
    )
  ;

  mod.add_type<SetOfChiralCentres>("SetOfChiralCentres")
    .method("size",
      [](const SetOfChiralCentres& s)->int64_t {
        return s.size();
      }
    )
    .method("internal_items_in_set",
      [](const SetOfChiralCentres& s)->int64_t {
        return s.size();
      }
    )
    .method("items_in_set",
      [] (const SetOfChiralCentres& s)->int64_t {
        return s.size();
      }
    )
  ;

  mod.set_override_module(jl_base_module);
  mod.method("getindex",
    [](const SetOfRings& r, int i)->const Ring*{
      // std::cerr << "getindex SetOfRings item " << i << " cmp " << r.size() << '\n';
      assert(i >= 1 && i <= static_cast<int>(r.size()));
      return r[i - 1];
    }
  );
  mod.method("getindex",
    [](const SetOfChiralCentres& r, int i)->const Chiral_Centre*{
      assert(i >= 1 && i <= static_cast<int>(r.size()));
      return r[i - 1];
    }
  );
  mod.unset_override_module();

  mod.add_type<Bond_list>("BondList")
    .constructor<>()
  ;
  mod.method("bonds_in_set",
    [](const Bond_list& b)->int64_t {
      return b.size();
    }
  );

  mod.set_override_module(jl_base_module);
  mod.method("getindex",
    [](const jlcxx::BoxedValue<Bond_list>& boxed_bond_list, int ndx)->const Bond*{
      const Bond_list& b = jlcxx::unbox<const Bond_list&>(boxed_bond_list);
      return b[ndx];
    }
  );
  mod.method("getindex",
    [](const Bond_list& blist, int i)->const Bond*{
      //std::cerr << "Bond_list::getindex length " << blist.size() << " ndx " << i << '\n';
      return blist[i - 1];
    }
  );
#ifdef DOES_NOT_SEEM_NECESSARY
  mod.method("getindex",
    [](Bond_list& blist, int i)->const Bond*{
      //std::cerr << "Bond_list::getindex length " << blist.size() << " ndx " << i << '\n';
      return blist[i];
    }
  );
#endif
  mod.unset_override_module();

  // Not sure why this was necessary, but for some reason the returned BondList object
  // generated an ambiguous getitem invocation, so this is a supporting function to
  // support that functionality. Do not use.
  mod.method("internal_get_item",
    [](const Bond_list& blist, int ndx)->const Bond* {
      return blist[ndx - 1];
    }
  );

  mod.add_type<RingAtoms>("RingAtoms")
    .constructor<>()
    .method("nrings", &RingAtoms::nrings)
  ;

  mod.set_override_module(jl_base_module);
  mod.method("getindex",
    [](const RingAtoms& rings, int ndx)->const Set_of_Atoms&{
      return rings[ndx - 1];
    }
  );
  mod.unset_override_module();

  mod.add_type<Molecule>("Molecule")
    .constructor<>()
    .constructor<Molecule>()
    .method("ok",
      [](const Molecule& m)->bool{
        return m.ok();
      }
    )
    .method("debug_string", &Molecule::debug_string)
    .method("empty",
      [](const Molecule& m)->bool{
        return m.empty();
      }
    )
    .method("x",
      [](const Molecule& m, atom_number_t a){
        return m.x(a);
      }
    )
    .method("y",
      [](const Molecule& m, atom_number_t a){
        return m.y(a);
      }
    )
    .method("z",
      [](const Molecule& m, atom_number_t a){
        return m.z(a);
      }
    )
    .method("setx!",
      [](Molecule& m, atom_number_t a, float x){
        return m.setx(a, x);
      }
    )
    .method("setx!",
      [](Molecule& m, atom_number_t a, double x) {
        return m.setx(a, x);
      }
    )
    .method("sety!",
      [](Molecule& m, atom_number_t a, float y){
        return m.sety(a, y);
      }
    )
    .method("sety!",
      [](Molecule& m, atom_number_t a, double y){
        return m.sety(a, y);
      }
    )
    .method("setz!",
      [](Molecule& m, atom_number_t a, float z){
        return m.setz(a, z);
      }
    )
    .method("setz!",
      [](Molecule& m, atom_number_t a, double z){
        return m.setz(a, z);
      }
    )
    .method("setxyz!", &Molecule::setz)
    // Maybe it would be more efficient to just extract the values atom at a time
    // and not do a transpose.
    .method("get_coordinates",
      [](const Molecule& m)-> jlcxx::ArrayRef<float, 2> {
        const int matoms = m.natoms();
        std::unique_ptr<float[]> coords = m.GetCoords();

        typedef Eigen::Matrix<float, Eigen::Dynamic, 3> myvector;
        Eigen::Map<myvector, Eigen::Aligned> as_eigen(coords.get(), matoms, 3);
        as_eigen.transposeInPlace();

        static constexpr bool kJuliaOwned = true;
        jlcxx::ArrayRef<float, 2> result(kJuliaOwned, coords.get(), matoms, 3);
        coords.release();
        return result;
      }
    )
    // Remember matrix ordering
    .method("set_xyz!",
      [](Molecule& m, const jlcxx::ArrayRef<float, 2> coords) {
        const int matoms = m.natoms();
        for (int i = 0; i < matoms; ++i) {
          m.setx(i, coords[i]);
          m.sety(i, coords[matoms + i]);
          m.setz(i, coords[matoms + matoms + i]);
        }
      }
    )
    .method("dihedral_scan",
      [](Molecule& m, atom_number_t a2, atom_number_t a3, double delta, double bump_check) -> jlcxx::ArrayRef<float, 3> {
        const int matoms = m.natoms();
        std::vector<std::unique_ptr<float[]>> b = m.DihedralScan(a2, a3, delta, bump_check);

        int nconf = b.size();
        std::unique_ptr<float[]> tmp = std::make_unique<float[]>(nconf * matoms * 3);

        for (int i = 0; i < nconf; ++i) {
          float* start = tmp.get() + i * matoms;
          for (int j = 0; j < matoms; ++j) {
            start[j] = b[i].get()[j * matoms];
            start[matoms + j] = b[i].get()[j * matoms + 1];
            start[matoms + matoms + j] = b[i].get()[j * matoms + 2];
          }
        }

        constexpr int kJuliaOwned = true;

        jlcxx::ArrayRef<float, 3> result(kJuliaOwned, tmp.get(), nconf * matoms * 3);
        tmp.release();

        return result;
      }
    )

    .method("jl_add_bond!",
      [](Molecule& m, atom_number_t a1, atom_number_t a2, BondType bt)->bool{
        return m.add_bond(a1, a2, BtypeEnumToBtype(bt));
      }
    )
    .method("are_bonded",
      [](const Molecule& m, atom_number_t a1, atom_number_t a2)->bool{
        return m.are_bonded(a1, a2);
      }
    )
    .method("formal_charge",
      [](const Molecule& m, atom_number_t a)->int64_t{
        return m.formal_charge(a);
      }
    )
    .method("isotope",
      [](const Molecule& m, atom_number_t a){
        return m.isotope(a);
      }
    )
    .method("set_isotope!",
      [](Molecule& m, atom_number_t zatom, isotope_t iso)->bool {
        return m.set_isotope(zatom, iso);
      }
    )
    .method("set_isotopes!",
      [](Molecule& m, const Set_of_Atoms& s, isotope_t iso) {
        return m.set_isotope(s, iso);
      }
    )
    .method("set_isotopes!",
      [](Molecule& m, const Set_of_Atoms* s, isotope_t iso) {
        return m.set_isotope(*s, iso);
      }
    )
    .method("number_isotopic_atoms",
      [](const Molecule& m)->int64_t {
        return m.number_isotopic_atoms();
      }
    )
    .method("number_isotopic_atoms",
      [](const Molecule& m, const isotope_t iso)->int64_t {
        return m.number_isotopic_atoms(iso);
      }
    )
    .method("remove_isotopes!",
      [](Molecule& m)->int64_t {
        return m.transform_to_non_isotopic_form();
      }
    )
    .method("number_formally_charged_atoms", &Molecule::number_formally_charged_atoms)
    .method("net_formal_charge", &Molecule::net_formal_charge)
    .method("number_formal_charges", &Molecule::number_formally_charged_atoms)
    .method("has_formal_charges", 
      [](const Molecule& m)->bool {
        return m.has_formal_charges();
      },
      "True if the molecule has any formal charges"
    )

    .method("set_name!",
      [](Molecule& m, const std::string& s) {
        return m.set_name(s);
      }
    )
    .method("name",
      [](const Molecule& m)->std::string{
        return m.name().AsString();
      }
    )
    .method("molecule_show_text",
      [](Molecule& m)->std::string{
        IWString result;
        result << m.name() << ' ' << m.natoms() << " atoms";
        return result.AsString();
      }
    )
    .method("molecule_show_text",
      [](Molecule* m)->std::string {
        IWString result;
        result << m->name() << ' ' << m->natoms() << " atoms";
        return result.AsString();
      }
    )
    .method("molecular_formula",
      [](Molecule& m)->std::string{
        IWString tmp;
        m.isis_like_molecular_formula_dot_between_fragments(tmp);
        return std::string(tmp.data(), tmp.size());
      }
    )
    .method("natoms",
      [](const Molecule& m)->int64_t{
        return m.natoms();
      }
    )
    .method("natoms",
      [](const Molecule& m, atomic_number_t z)->int64_t {
        return m.natoms(z);
      }
    )
    .method("natoms",
      [](const Molecule& m, std::string& asymbol)->int64_t{
        return m.natoms(asymbol.c_str());
      }
    )
    .method("nedges", &Molecule::nedges)
    .method("atomic_number",
      [](const Molecule& m, atom_number_t a) {
        return m.atomic_number(a);
      }
    )
    .method("atomic_symbol",
      [](const Molecule& m, atom_number_t a)->std::string{
        return m.atomic_symbol(a).AsString();
      }
    )
    .method("nrings",
      [](Molecule& m)->int64_t {
        return m.nrings();
      }
    )
    .method("nrings",
      [](Molecule& m, atom_number_t a)->int64_t {
        return m.nrings(a);
      }
    )
    .method("nrings",
      [](Molecule& m, atom_number_t a, int rsize)->int64_t {
        return m.nrings(a, rsize);
      }
    )
    .method("is_ring_atom",
      [](Molecule& m, atom_number_t a)->bool{
        return m.is_ring_atom(a);
      }
    )
    .method("ring_bond_count", 
      [](Molecule& m, atom_number_t a)->int64_t {
        return m.ring_bond_count(a);
      }
    )
    .method("in_ring_of_given_size", 
      [](Molecule& m, atom_number_t zatom, int rsize)->bool {
        return m.in_ring_of_given_size(zatom, rsize);
      }
    )
    .method("largest_ring_size", 
      [](Molecule & m)->int64_t {
        return m.LargestRingSize();
      }
    )
    .method("fused_system_size",
      [](Molecule& m, atom_number_t a)->int64_t {
        return m.fused_system_size(a);
      }
    )
    .method("rings_with_fused_system_identifier", &Molecule::rings_with_fused_system_identifier)
    .method("fused_system_identifier",
      [](Molecule& m, atom_number_t a)->int64_t {
        return m.fused_system_identifier(a);
      }
    )
    .method("in_same_ring",
      [](Molecule& m, atom_number_t a1, atom_number_t a2)->bool{
        return m.in_same_ring(a1, a2);
      }
    )
    .method("in_same_aromatic_ring",
      [](Molecule& m, atom_number_t a1, atom_number_t a2)->bool{
        return m.in_same_aromatic_ring(a1, a2);
      }
    )
    .method("in_same_ring_system",
      [](Molecule& m, atom_number_t a1, atom_number_t a2)->bool{
        return m.in_same_ring_system(a1, a2);
      }
    )
    .method("number_ring_systems", &Molecule::number_ring_systems)
    .method("ring_membership",
      [](Molecule& m)->std::vector<int>{
        const int* r = m.ring_membership();
        return std::vector<int>(r, r + m.natoms());
      }
    )
    .method("rings_containing_both",
      [](Molecule& m, atom_number_t a1, atom_number_t a2) {
        return m.in_same_rings(a1, a2);
      }
    )
    .method("is_part_of_fused_ring_system",
      [](Molecule& m, atom_number_t a)->bool{
        return m.is_part_of_fused_ring_system(a);
      }
    )
    .method("ring",
      [](Molecule& m, int rnum)->const Ring*{
        return m.ringi(rnum);
      }
    )
    .method("ring_containing_atom",
      [](Molecule& m, atom_number_t a)->Ring{
        const Ring* r = m.ring_containing_atom(a);
        if (r == nullptr) {
          return Ring();
        }
        return *r;
      }
    )

    .method("sssr_rings",
      [](Molecule& m)->SetOfRings{
        SetOfRings result(m.sssr_rings());
        return result;
      }
    )

    // Same as sssr_rings
    .method("rings",
      [](Molecule& m)->SetOfRings{
        SetOfRings result(m.sssr_rings());
        return result;
      }
    )

    .method("label_atoms_by_ring_system",
      [](Molecule& m)->std::vector<int>{
        const int matoms = m.natoms();
        std::unique_ptr<int[]> tmp = std::make_unique<int[]>(matoms);
        m.label_atoms_by_ring_system(tmp.get());
        return std::vector<int>(tmp.get(), tmp.get() + matoms);
      }
    )
    .method("label_atoms_by_ring_system_including_spiro_fused",
      [](Molecule& m)->std::vector<int>{
        const int matoms = m.natoms();
        std::unique_ptr<int[]> tmp = std::make_unique<int[]>(matoms);
        m.label_atoms_by_ring_system_including_spiro_fused(tmp.get());
        return std::vector<int>(tmp.get(), tmp.get() + matoms);
      }
    )
    .method("nrings_including_non_sssr_rings",
      [](Molecule& m, atom_number_t a){
        return m.nrings_including_non_sssr_rings(a);
      }
    )
    .method("non_sssr_rings", &Molecule::non_sssr_rings)
    .method("non_sssr_ring",
      [](Molecule& m, int rnum)->Ring{
        return *m.non_sssr_ring(rnum);
      }
    )
    .method("is_spiro_fused",
      [](Molecule& m, atom_number_t a)->bool{
        return m.is_spiro_fused(a);
      }
    )
    .method("is_halogen", 
      [](const Molecule& m, atom_number_t a)->bool{
        return m.is_halogen(a);
      }
    )
    .method("ncon",
      [](const Molecule& m, atom_number_t a)->int64_t {
        return m.ncon(a);
      }
    )
    .method("nbonds",
      [](const Molecule& m, atom_number_t a)->int64_t {
        return m.nbonds(a);
      }
    )
    .method("maximum_connectivity", &Molecule::maximum_connectivity)
    .method("other",
      [](Molecule& m, atom_number_t atom, int ndx){
        return m.other(atom, ndx);
      }
    )
    .method("connections",
      [](const Molecule& m, atom_number_t a)->std::vector<atom_number_t>{
        std::vector<atom_number_t> result;
        for (atom_number_t o : m.connections(a)) {
          result.push_back(o);
        }
        return result;
      }
    )

    .method("build_from_smiles", 
      [](Molecule& m, const std::string& s)->bool{
        return m.build_from_smiles(s);
      }
    )
    .method("smiles", 
      [](Molecule& m)->std::string{
        return m.smiles().AsString();
      }
    )
    .method("internal_smiles",
      [](Molecule* m)->std::string {
        return m->smiles().AsString();
      }
    )
    .method("unique_smiles", 
      [](Molecule& m)->std::string{
        return m.unique_smiles().AsString();
      }
    )
    .method("random_smiles", 
      [](Molecule& m)->std::string{
        return m.random_smiles().AsString();
      }
    )
    .method("smiles_starting_with_atom", 
      [](Molecule& m, atom_number_t zatom)->std::string{
        return m.smiles_starting_with_atom(zatom).AsString();
      }
    )
    .method("unique_kekule_smiles", 
      [](Molecule& m)->std::string{
        return m.UniqueKekuleSmiles().AsString();
      }
    )
    .method("isotopically_labelled_smiles",
      [](Molecule& m)->std::string{
        return m.isotopically_labelled_smiles().AsString();
      }
    )
    .method("aromatic_smiles",
      [](Molecule& m)->std::string {
        return m.aromatic_smiles().AsString();
      },
      "non-unique smiles including aromaticity"
    )
    .method("is_aromatic", 
      [](Molecule& m, atom_number_t a)->bool{
        return m.is_aromatic(a);
      }
    )
    .method("atom", &Molecule::atom)
    .method("set_atomic_number!",
      [](Molecule& m, atom_number_t zatom, atomic_number_t q) {
        return m.set_atomic_number(zatom, q);
      }
    )
    .method("occursin",
      [](const std::string& s, const Molecule& m)->bool {
        const Element* e = get_element_from_symbol_no_case_conversion(s.data(), s.size());
        if (e == nullptr) {
          return false;
        }
        for (const Atom* a : m) {
          if (a->element() == e) {
            return true;
          }
        }
        return false;
      }
    )
    .method("contains",
      [](const Molecule& m, const std::string& s)->bool {
        const Element* e = get_element_from_symbol_no_case_conversion(s.data(), s.size());
        if (e == nullptr) {
          return false;
        }
        for (const Atom* a : m) {
          if (a->element() == e) {
            return true;
          }
        }
        return false;
      }
    )

    .method("set_bond_type_between_atoms!",
      [](Molecule& m, atom_number_t a1, atom_number_t a2, BondType bt) {
        //std::cerr << "Ttyoe will be " << BtypeEnumToBtype(bt) << '\n';
        return m.set_bond_type_between_atoms(a1, a2, BtypeEnumToBtype(bt));
      }
    )
    .method("change_to_graph_form",
      [](Molecule& m) {
        return m.change_to_graph_form();
      }
    )
    .method("to_scaffold!",
      [](Molecule& m) {
        return m.ToScaffold();
      }
    )
    .method("compute_aromaticity_if_needed", &Molecule::compute_aromaticity_if_needed)
    .method("change_to_graph_form",
      [](Molecule& m, const Mol2Graph& mol2graph) {
        return m.change_to_graph_form(mol2graph);
      }
    )
    .method("smiles_atom_order",
      [](Molecule& m)->std::vector<int>{
        std::vector<int> result(m.natoms());
        m.smiles_atom_order(result.data());
        return result;
      }
    )
    .method("atom_order_in_smiles",
      [](Molecule& m)->std::vector<int>{
        const resizable_array<int>& order = m.atom_order_in_smiles();
        const int matoms = m.natoms();
        std::vector<int> result;
        result.reserve(matoms);
        for (int i = 0; i < matoms; ++i) {
          result.push_back(order[i]);
        }
        return result;
      }
    )
    .method("bond", 
      [](const Molecule& m, int ndx) {
        return *m.bondi(ndx);
      }
    )
    .method("bond_between_atoms",
      [](Molecule& m, atom_number_t a1, atom_number_t a2){
        return *m.bond_between_atoms(a1, a2);
      }
    )

    //.method("compute_canonical_ranking", &Molecule::compute_canonical_ranking)
    .method("canonical_rank",
      [](Molecule& m, atom_number_t a) {
        return m.canonical_rank(a);
      }
    )
    .method("canonical_ranks",
      [](Molecule& m)->std::vector<int> {
        const int * c = m.canonical_ranks();
        return std::vector<int>(c, c + m.natoms());
      }
    )
    .method("symmetry_class", 
      [](Molecule& m, atom_number_t a) {
        return m.symmetry_class(a);
      }
    )
    .method("number_symmetry_classes", &Molecule::number_symmetry_classes)
    .method("symmetry_equivalents",
      [](Molecule& m, atom_number_t a)->Set_of_Atoms {
        Set_of_Atoms result;
        m.symmetry_equivalents(a, result);
        return result;
      }
    )
    .method("symmetry_classes",
      [](Molecule& m)->std::vector<int>{
        std::vector<int> result;
        const int matoms = m.natoms();
        result.reserve(matoms);
        const int* sym = m.symmetry_classes();
        std::copy(sym, sym + matoms, std::back_inserter(result));
        return result;
      }
    )
    .method("attached_heteroatom_count", 
      [](const Molecule& m, atom_number_t a) {
        return m.attached_heteroatom_count(a);
      }
    )
    //.method("multiple_bond_to_heteroatom", &Molecule::multiple_bond_to_heteroatom)

    .method("bond_length",
      [](const Molecule& m, atom_number_t a1, atom_number_t a2){
        return m.bond_length(a1, a2);
      }
    )
    .method("bond_angle",
      [](const Molecule& m, atom_number_t a1, atom_number_t a2, atom_number_t a3){
        return m.bond_angle(a1, a2, a3);
      }
    )
    .method("dihedral_angle",
      [](const Molecule& m, atom_number_t a1, atom_number_t a2, atom_number_t a3, atom_number_t a4){
        return m.dihedral_angle(a1, a2, a3, a4);
      }
    )
    .method("signed_dihedral_angle",
      [](const Molecule& m, atom_number_t a1, atom_number_t a2, atom_number_t a3, atom_number_t a4){
        return m.signed_dihedral_angle(a1, a2, a3, a4);
      }
    )
    .method("set_bond_length",
      [](Molecule& m, atom_number_t a1, atom_number_t a2, float dist){
        return m.set_bond_length(a1, a2, dist);
      }
    )
    .method("set_bond_angle",
      [](Molecule& m, atom_number_t a1, atom_number_t a2, atom_number_t a3, angle_t angle){
        return m.set_bond_angle(a1, a2, a3, angle);
      }
    )
    .method("set_dihedral_angle",
      [](Molecule& m, atom_number_t a1, atom_number_t a2, atom_number_t a3, atom_number_t a4, angle_t angle){
        return m.set_dihedral(a1, a2, a3, a4, angle);
      }
    )
    .method("bond_list",
      [](const Molecule& m)->const Bond_list& {
        return m.bond_list();
      }
    )
    .method("jl_bonds",
      [](const Molecule& m)->const Bond_list& {
        return m.bond_list();
      }
    )
//  .method("bonds",
//    [](const Molecule* m)->const Bond_list& {
//      return m->bond_list();
//    }
//  )
    .method("add!",
      [](Molecule& m, const Molecule& rhs) {
        m.add_molecule(&rhs);
      }
    )
    .method("add_atom!",
      [](Molecule& m, atomic_number_t z)->bool {
        const Element* e = get_element_from_atomic_number(z);
        if (e == nullptr) {
          return false;
        }

        m.add(e);

        return true;
      }
    )
    .method("add!",
      [](Molecule& m, atomic_number_t z)->bool {
        const Element* e = get_element_from_atomic_number(z);
        if (e == nullptr) {
          return false;
        }

        m.add(e);

        return true;
      }
    )
    .method("has_partial_charges",
      [](const Molecule& m)->bool{
        return m.has_partial_charges();
      }
    )
    .method("set_formal_charge!",
      [](Molecule& m, atom_number_t a, formal_charge_t q) {
        return m.set_formal_charge(a, q);
      }
    )
    .method("formal_charge", 
      [](Molecule& m, atom_number_t a) {
        return m.formal_charge(a);
      }
    )
    .method("remove_atom!",
      [](Molecule& m, atom_number_t a) {
        return m.remove_atom(a);
      }
    )
    .method("remove_atoms!",
      [](Molecule& m, Set_of_Atoms& s) {
        return m.remove_atoms(s);
      }
    )
    .method("remove_atoms!",
      [](Molecule& m, const jlcxx::ArrayRef<int64_t> to_remove) {
        return m.remove_atoms(to_remove.data());
      }
    )
    .method("delete_fragment!",
      [](Molecule& m, int frag) {
        return m.delete_fragment(frag);
      }
    )
    .method("remove_fragment_containing_atom!",
      [](Molecule& m, atom_number_t a) {
        return m.remove_fragment_containing_atom(a);
      }
    )
    .method("remove_all!",
      [](Molecule& m, atomic_number_t z) {
        return m.remove_all(z);
      }
    )
    .method("remove_all_non_natural_elements!", &Molecule::remove_all_non_natural_elements)
    .method("valence_ok",
      [](Molecule& m)->bool{
        return m.valence_ok();
      }
    )
    .method("remove_explicit_hydrogens!", &Molecule::remove_explicit_hydrogens)
    .method("chop!", &Molecule::chop)
    .method("remove_bonds_to_atom!",
      [](Molecule& m, atom_number_t a) {
        return m.remove_bonds_to_atom(a);
      }
    )
    .method("remove_bond!", &Molecule::remove_bond)
    .method("remove_bond_between_atoms!", &Molecule::remove_bond_between_atoms)
    .method("remove_all_bonds!", &Molecule::remove_all_bonds)
    .method("molecular_weight",
      [](const Molecule& m) {
        return m.molecular_weight();
      }
    )
    .method("amw",
      [](const Molecule& m)->molecular_weight_t {
        return m.molecular_weight();
      }
    )
    .method("molecular_weight_count_isotopes", &Molecule::molecular_weight_count_isotopes)
    .method("molecular_weight_ignore_isotopes", &Molecule::molecular_weight_ignore_isotopes)
    .method("highest_coordinate_dimensionality", &Molecule::highest_coordinate_dimensionality)
    .method("exact_mass",
      [](const Molecule& m) {
        return m.exact_mass();
      }
    )

    .method("translate",
      [](Molecule& m, float dx, float dy, float dz) {
        return m.translate_atoms(dx, dy, dz);
      }
    )
    .method("number_fragments",
      [](Molecule& m) {
        return m.number_fragments();
      }
    )
    .method("fragment_membership",
      [](Molecule& m, atom_number_t a) {
        return m.fragment_membership(a);
      }
    )
    .method("fragment_membership",
      [](Molecule& m, jlcxx::ArrayRef<int64_t> fragm) {
        std::unique_ptr<int[]> f = m.fragment_membership();
        std::copy(f.get(), f.get() + m.natoms(), fragm.data());
        return 1;
      }
    )
    .method("fragment_membership",
      [](Molecule& m) -> jlcxx::ArrayRef<int32_t> {
        std::unique_ptr<int32_t[]> f = m.fragment_membership();
        jlcxx::ArrayRef<int32_t> result(true, f.get(), m.natoms());
        f.release();
        return result;
      }
    )
    .method("atoms_in_fragment", 
      [](Molecule& m, int frag) {
        return m.atoms_in_fragment(frag);
      }
    )
    .method("get_atoms_in_fragment",
      [](Molecule& m, int frag)->Set_of_Atoms{
        Set_of_Atoms result;
        m.atoms_in_fragment(result, frag);
        return result;
      }
    )
    .method("largest_fragment", &Molecule::largest_fragment)
    .method("identify_spinach",
      [](Molecule& m)->std::vector<int64_t>{
        const int matoms = m.natoms();
        std::unique_ptr<int[]> spch = std::make_unique<int[]>(matoms);
        m.identify_spinach(spch.get());
        std::vector<int64_t> result(matoms);
        for (int i = 0; i < matoms; ++i) {
          result[i] = spch[i];
        }
        return result;
      }
    )
    .method("rings_in_fragment", &Molecule::rings_in_fragment)

    .method("create_subset",
      [](Molecule& m, jlcxx::ArrayRef<int64_t> subset)->Molecule{
        const int matoms = m.natoms();
        std::unique_ptr<int[]> tmp = std::make_unique<int[]>(matoms);
        for (int i = 0; i < matoms; ++i) {
          if (subset[i] > 0) {
            tmp[i] = 1;
          } else {
            tmp[i] = 0;
          }
        }
        Molecule result;
        m.create_subset(result, tmp.get(), 1);
        return result;
      }
    )
    .method("create_subset",
      [](Molecule& m, const Set_of_Atoms& subset){
        return m.create_subset(subset);
      }
    )
    .method("reduce_to_largest_fragment!", &Molecule::reduce_to_largest_fragment)
    .method("reduce_to_largest_organic_fragment!", &Molecule::reduce_to_largest_organic_fragment)
    .method("reduce_to_largest_fragment_carefully!", &Molecule::reduce_to_largest_fragment_carefully)

    .method("contains_non_periodic_table_elements",
      [](const Molecule& m)->bool {
        return m.contains_non_periodic_table_elements();
      }
    )
    .method("organic_only",
      [](const Molecule& m)->bool {
        return m.organic_only();
      }
    )
    .method("longest_path", &Molecule::longest_path)
    .method("bonds_between",
      [](Molecule& m, atom_number_t a1, atom_number_t a2) {
        return m.bonds_between(a1, a2);
      }
    )
    .method("atoms_between",
      [](Molecule& m, atom_number_t a1, atom_number_t a2)->Set_of_Atoms{
        Set_of_Atoms result;
        m.atoms_between(a1, a2, result);
        return result;
      }
    )
    .method("most_distant_pair",
      [](Molecule& m)->std::tuple<int, int> {
        atom_number_t a1 = INVALID_ATOM_NUMBER;
        atom_number_t a2 = INVALID_ATOM_NUMBER;
        int longest_distance = 0;
        const int matoms = m.natoms();
        for (int i = 0; i < matoms; ++i) {
          for (int j = i + 1; j < matoms; ++j) {
            if (m.fragment_membership(i) != m.fragment_membership(j)) {
              continue;
            }
            const int d = m.bonds_between(i, j);
            if (d > longest_distance) {
              a1 = i;
              a2 = j;
              longest_distance = d;
            }
          }
        }
        return std::make_tuple(a1, a2);
      }
    )
    .method("implicit_hydrogens",
      [](Molecule& m, atom_number_t a) {
        return m.implicit_hydrogens(a);
      }
    )
    .method("unset_all_implicit_hydrogen_information", &Molecule::unset_all_implicit_hydrogen_information)
    .method("remove_hydrogens_known_flag_to_fix_valence_errors", &Molecule::remove_hydrogens_known_flag_to_fix_valence_errors)
    .method("explicit_hydrogens", &Molecule::explicit_hydrogens)
    .method("hcount", &Molecule::hcount)
    .method("make_implicit_hydrogens_explicit!",
      [](Molecule& m) {
        return m.make_implicit_hydrogens_explicit();
      }
    )
    .method("move_hydrogens_to_end_of_connection_table!", &Molecule::move_hydrogens_to_end_of_connection_table)
    .method("pi_electrons",
      [](Molecule& m, atom_number_t a) {
        int result;
        if (! m.pi_electrons(a, result)) {
          return 0;
        }
        return result;
      }
    )
    .method("lone_pair_count",
      [](Molecule& m, atom_number_t a) {
        int result;
        if (! m.lone_pair_count(a, result)) {
          return 0;
        }

        return result;
      }
    )
    .method("aromatic_atom_count", &Molecule::aromatic_atom_count)
    .method("aromatic_ring_count", &Molecule::aromatic_ring_count)
    .method("distance_between_atoms", &Molecule::distance_between_atoms)
    .method("longest_intra_molecular_distance", &Molecule::longest_intra_molecular_distance)
    .method("saturated", 
      [](const Molecule& m, atom_number_t a)->bool{
        return m.saturated(a);
      }
    )
    .method("unsaturation", &Molecule::unsaturation)
    .method("chiral_centres",
      [](Molecule& m) -> const SetOfChiralCentres& {
        return m.ChiralCentres();
      }
    )
    .method("number_chiral_centres", &Molecule::chiral_centres)
    .method("chiral_centre_at_atom", &Molecule::chiral_centre_at_atom)
    .method("chiral_centre_in_molecule_not_indexed_by_atom_number", &Molecule::chiral_centre_in_molecule_not_indexed_by_atom_number)
    .method("remove_chiral_centre_at_atom!", &Molecule::remove_chiral_centre_at_atom)
    .method("remove_all_chiral_centres!", &Molecule::remove_all_chiral_centres)
    .method("invert_chirality_on_atom!", &Molecule::invert_chirality_on_atom)
    .method("smarts_equivalent_for_atom",
      [](Molecule& m, atom_number_t a)->std::string{
        return m.smarts_equivalent_for_atom(a).AsString();
      }
    )
    .method("smarts",
      [](Molecule& m)->std::string{
        return m.smarts().AsString();
      }
    )
    .method("atom_map_number", &Molecule::atom_map_number)
    .method("set_atom_map_number!", &Molecule::set_atom_map_number)
    .method("atom_with_atom_map_number", &Molecule::atom_with_atom_map_number)
    .method("reset_all_atom_map_numbers!", &Molecule::reset_all_atom_map_numbers)
    .method("reset_atom_map_numbers!", &Molecule::reset_all_atom_map_numbers)

    .method("unset_unnecessary_implicit_hydrogens_known_values!", &Molecule::unset_unnecessary_implicit_hydrogens_known_values)
    .method("discern_chirality_from_3d_structure", &Molecule::discern_chirality_from_3d_structure)

    .method("sort_atoms!",
      [](Molecule& m, jlcxx::ArrayRef<int64_t> order, int direction){
        return m.sort(order.data(), direction);
      }
    )


  ;

  mod.method("gather_rings",
    [](Molecule& m, RingAtoms& ring_atoms)->int {
      return ring_atoms.GatherRings(m);
    },
    "Copy the rings from `m` to `ring_atoms`"
  );

  mod.method("gather_rings",
    [](Molecule& m, lillymol::RingInformation& rinf, lillymol::RIProperty& properties)  {
      return rinf.GatherRings(m, properties);
    },
    "Copy ring atoms and optionally ring properties from `m` to ring_information"
  );

  mod.set_override_module(jl_base_module);
  mod.method("getindex",
    [](const lillymol::RingInformation& rings, int ndx)->const Set_of_Atoms&{
      return rings[ndx - 1];
    }
  );
  mod.unset_override_module();

  mod.add_type<Components>("Components")
    .constructor<>()
    .method("length",
      [](const Components& c) {
        return c.size();
      },
      ""
      )
  ;

  mod.set_override_module(jl_base_module);
  mod.method("getindex",
    [](const Components& components, int ndx)->Molecule&{
      return *components[ndx - 1];
    }
  );
  mod.unset_override_module();

  mod.method("create_components",
    [](Molecule& m, Components& components) {
      return m.create_components(components);
    },
    ""
  );
  mod.method("create_components_internal",
    [](Molecule& m, Components& components) {
      return m.create_components(components);
    },
    ""
  );
  mod.method("jl_create_components",
    [](Molecule& m)->Components {
      Components result;
      m.create_components(result);
      return result;
    },
    ""
  );

  mod.method("set_auto_create_new_elements", &set_auto_create_new_elements);
  mod.method("set_atomic_symbols_can_have_arbitrary_length", &set_atomic_symbols_can_have_arbitrary_length);
  mod.method("set_include_atom_map_with_smiles", &set_include_atom_map_with_smiles);
  mod.method("set_copy_name_in_molecule_copy_constructor", &set_copy_name_in_molecule_copy_constructor);
  mod.method("set_display_smiles_interpretation_error_messages", 
    [](int s) {
      set_display_smiles_interpretation_error_messages(s);
    },
    "Suppress display of smiles parsing errors"
  );

  mod.method("returns_vector",
    [](int i, int j)->std::vector<int64_t>{
      std::vector<int64_t> result;
      result.push_back(i);
      result.push_back(j);
      return result;
    }
  );

  mod.set_override_module(jl_base_module);
  mod.method("getindex",
    [](const Molecule& m, int i)->const Atom&{
      return m[i];
    }
  );
  mod.unset_override_module();

  mod.method("MolFromSmiles",
    [](const std::string& smiles)->Molecule{
      Molecule result;
      result.build_from_smiles(smiles);
      return result;
    }
  );
  mod.method("LillyMolFromSmiles",
    [](const std::string& smiles)->Molecule{
      Molecule result;
      result.build_from_smiles(smiles);
      return result;
    }
  );

  mod.method("xlogp",
    [](Molecule& m) ->std::tuple<float, bool>{
      auto rc = xlogp::XLogP(m);
      if (rc) {
        return std::make_tuple<float, bool>(*rc, true);
      }
      return std::make_tuple<float, bool>(0.0f, false);
    }
  );

  mod.add_type<alogp::ALogP>("alogp")
    .method("logp",
      [](alogp::ALogP& a, Molecule& m)->float {
        std::optional<float> res = a.LogP(m);
        if (! res) {
          // Silently ignore failures for now...
          return 0.0f;
        }
        return *res;
      },
      "Computes alogp"
    )
  ;

  // Substructre Related
  mod.add_type<SetOfEmbeddings>("SetOfEmbeddings")
    .constructor<>()
    .method("size",
      [](const SetOfEmbeddings& r) {
        return r.size();
      }
    )
    .method("items_in_set",
      [](const SetOfEmbeddings&r) {
        return r.size();
      }
    )
    .method("set_of_embeddings_show_text",
      [](const SetOfEmbeddings& soe)->std::string {
        IWString result;
        result << "SetOfEmbeddings with " << soe.size() << " embeddings";
        return result.AsString();
      }
    )
  ;

  mod.set_override_module(jl_base_module);
  mod.method("getindex",
    [](const SetOfEmbeddings& s, int i)->const Set_of_Atoms*{
      return s[i - 1];
    }
  );
  mod.unset_override_module();

  mod.method("set_of_embeddings_getindex",
    [](const SetOfEmbeddings& s, int ndx)->const Set_of_Atoms* {
      return s[ndx - 1];
    }
  );

  //jlcxx::stl::apply_stl<Molecule*>(mod);
  mod.add_type<Molecule_to_Query_Specifications>("MoleculeToQuerySpecifications")
    .constructor<>()
  ;

  mod.add_type<Substructure_Results>("SubstructureResults")
    .constructor<>()
    .method("number_embeddings", &Substructure_Results::number_embeddings)
    .method("embeddings",
      [](const Substructure_Results& sresults)->const SetOfEmbeddings& {
        return sresults.embeddings();
      },
      "Iterable collection of the embeddings found during substructure search"
    )
    .method("each_embedding_set_vector",
      [](const Substructure_Results& sresults, int n, int value)->jlcxx::ArrayRef<int> {
        std::unique_ptr<int[]>tmp = std::make_unique<int[]>(n);
        std::fill_n(tmp.get(), n, 0);
        sresults.each_embedding_set_vector(tmp.get(), value);

        jlcxx::ArrayRef<int> result(true, tmp.get(), n);
        tmp.release();
        return result;
      },
      "Generate array with all matched atoms set to a fixed value"
    )
  ;

  mod.add_type<Molecule_to_Match>("MoleculeToMatch")
    .constructor<Molecule*>()
  ;

  mod.add_type<Substructure_Query>("SubstructureQuery")
    .constructor<>() 
    .method("name",
      [](const Substructure_Query& q)->std::string {
        return q.comment().AsString();
      },
      "Name of query"
    )
    .method("set_only_keep_matches_in_largest_fragment!",
      [](Substructure_Query& q, bool s) {
        q.set_only_keep_matches_in_largest_fragment(s);
      },
      "Only retain matched found in the largest frgament"
    )
    .method("set_embeddings_can_overlap!",
      [](Substructure_Query& q, bool s) {
        q.set_embeddings_do_not_overlap(!s);
      },
      "By default, embeddings can overlap"
    )
    .method("read_msi",
      [](Substructure_Query& q, const std::string& fname)->bool {
        IWString tmp(fname);
        return q.read(tmp);
      },
      "Read query in MSI format"
    )
    .method("read_proto",
      [](Substructure_Query& q, const std::string& fname)->bool {
        IWString tmp(fname);
        return q.ReadProto(tmp);
      },
      "Read query from textproto file"
    )
    .method("write_msi",
      [](Substructure_Query& q, const std::string& fname)->bool {
        IWString tmp(fname);
        return q.write_msi(tmp);
      },
      "Write as MSI form"
    )
    .method("write_proto",
      [](Substructure_Query& q, const std::string& fname)->bool {
        IWString tmp(fname);
        return q.WriteProto(tmp);
      },
      "Write as textproto form"
    )

    .method("set_find_one_embedding_per_atom!", &Substructure_Query::set_find_one_embedding_per_atom)
    .method("set_find_unique_embeddings_only!", &Substructure_Query::set_find_unique_embeddings_only)
    .method("set_min_matches_to_find!", &Substructure_Query::set_min_matches_to_find)
    .method("set_max_matches_to_find!", &Substructure_Query::set_max_matches_to_find)
    .method("set_perceive_symmetry_equivalent_matches!",
      [](Substructure_Query& q, bool s)->void {
        q.set_perceive_symmetry_equivalent_matches(s);
      },
      "Symmetry related matches are found by default"
    )
    .method("set_save_matched_atoms!",
      [](Substructure_Query& q, bool s)->void {
        q.set_save_matched_atoms(s);
      },
      "Searches can be faster if matched atoms are not saved, but this inhibits other functionality"
    )
    .method("max_query_atoms_matched_in_search", &Substructure_Query::max_query_atoms_matched_in_search)

    .method("matches",
      [](Substructure_Query& q, Molecule& m)->bool {
        return q.substructure_search(&m) > 0;
      },
      "True if the query matches `m`"
    )
    .method("substructure_search",
      [](Substructure_Query& q, Molecule& m)->uint32_t {
        return q.substructure_search(&m);
      },
      "The number of matches to `q` in `m`"
    )
    .method("substructure_search",
      [](Substructure_Query& q, Molecule& m, Substructure_Results& sresults) {
        return q.substructure_search(m, sresults);
      },
      "Perform substructure search, embeddings in `sresults`"
    )
    .method("substructure_search",
      [](Substructure_Query& q, Molecule_to_Match& target, Substructure_Results& sresults) {
        return q.substructure_search(target, sresults);
      },
      "Perform substructure search, embeddings in `sresults`"
    )

    .method("build_from_smiles",
      [](Substructure_Query& q, const std::string& smiles)->bool{
        IWString tmp(smiles);
        return q.create_from_smiles(tmp);
      },
      "interpret `smiles` as smarts - mostly works but there are gotcha's"
    )
    .method("build_from_smarts",
      [](Substructure_Query& q, const std::string& smarts)->bool{
        IWString tmp(smarts);
        return q.create_from_smarts(tmp);
      },
      "Build query from `smarts`"
    )
    .method("build_from_molecule",
      [](Substructure_Query& q, Molecule& m, Molecule_to_Query_Specifications& mqs) {
        return q.create_from_molecule(m, mqs);
      },
      "Create query from molecule under control of `mqs`"
    )
    .method("substructure_search_as_vector",
      [](Substructure_Query& q, Molecule& m)->jlcxx::ArrayRef<int32_t, 2>{
        Substructure_Results sresults;
        const int nhits = q.substructure_search(m, sresults);
        if (nhits == 0) {
          return jlcxx::ArrayRef<int32_t, 2>(nullptr, 0, 0);
        }

        const int matched_atoms = sresults.embedding(0)->size();
        std::unique_ptr<int[]> tmp = std::make_unique<int[]>(nhits * matched_atoms);

#ifdef WHEN_WE_HAVE_CPP20
        for (auto [i, embedding] : std::views::enumerate(sresults.embeddings())) {
          for (auto [j, a] : std::views::enumerate(*embedding)) {
            tmp[j * nhits + i] = a;
          }
        }
#else
        for (int i = 0; i < nhits; ++i) {
          const Set_of_Atoms* embedding = sresults.embedding(i);
          for (int j = 0; j < embedding->number_elements(); ++j) {
            tmp[j * nhits + i] = embedding->item(j);
            // std::cerr << i << ' ' << j << " fill " << (j * matched_atoms + i) << " value " << embedding->item(j) << '\n';
          }
        }
#endif

        static constexpr int kJuliaOwned = true;

        jlcxx::ArrayRef<int32_t, 2> result(kJuliaOwned, tmp.get(), nhits, matched_atoms);
        tmp.release();
        return result;
      },
      "Return atoms hit as a matrix"
    )
  ;


  mod.add_type<Reaction_Iterator>("ReactionIterator")
    .constructor<>()

    .method("active",
      [](const Reaction_Iterator& rxnit)->bool{
        return rxnit.active();
      },
      "True if still active"
    )
    .method("increment",
      [](Reaction_Iterator& rxnit) {
        rxnit++;
      },
      "Move to next reagent"
    )
    .method("reagent", &Reaction_Iterator::reagent)
  ;


  mod.add_type<Sidechain_Match_Conditions>("SidechainMatchConditions")
    .constructor<>()
    .method("set_make_new_reagent_for_each_hit", &Sidechain_Match_Conditions::set_make_new_reagent_for_each_hit)
    .method("set_max_matches_to_find", &Sidechain_Match_Conditions::set_max_matches_to_find)
    .method("set_strip_reagents_to_largest_fragment", &Sidechain_Match_Conditions::set_strip_reagents_to_largest_fragment)
    .method("set_ignore_not_reacting", &Sidechain_Match_Conditions::set_ignore_not_reacting)
    .method("set_find_unique_embeddings_only", &Sidechain_Match_Conditions::set_find_unique_embeddings_only)
    .method("set_one_embedding_per_start_atom", &Sidechain_Match_Conditions::set_one_embedding_per_start_atom)
    .method("set_ignore_symmetry_related_matches", &Sidechain_Match_Conditions::set_ignore_symmetry_related_matches)
  ;

  mod.add_type<IWReaction>("Reaction")
    .constructor<>() 
    .method("name",
      [](const IWReaction& q)->std::string {
        return q.comment().AsString();
      },
      "Name of reaction"
    )
    .method("number_reagents", &IWReaction::number_reagents)
    .method("number_sidechains", &IWReaction::number_sidechains)
    .method("set_one_embedding_per_start_atom", &IWReaction::set_one_embedding_per_start_atom)
    .method("add_sicechain_reagents",
      [](IWReaction& rxn, int sidechain, const char* fname, FileType file_type, Sidechain_Match_Conditions& smc)->bool{
        return rxn.add_sidechain_reagents(sidechain, fname, file_type, smc);
      },
      "Add reagents to a sidechain"
    )
    .method("add_sidechain",
      [](IWReaction& rxn, int sidechain, const Molecule& m, Sidechain_Match_Conditions& smc)->bool {
        Sidechain_Reaction_Site* s = rxn.sidechain(sidechain);
        if (s == nullptr) {
          return false;
        }
        return s->add_reagent(m, smc);
      },
      "Add a sidechain reagent"
    )
    .method("read_textproto",
      [](IWReaction& rxn, const std::string& fname)->bool{
        IWString tmp(fname);
        return rxn.Read(tmp);
      },
      "read textproto reaction"
    )
    .method("substructure_search",
      [](IWReaction& rxn, Molecule& m, Substructure_Results& sresults) {
        return rxn.substructure_search(m, sresults);
      }
    )
    .method("in_place_transformation",
      [](IWReaction& rxn, Molecule& m)->bool{
        return rxn.in_place_transformations(m);
      },
      "apply reaction to 'm'"
    )
    .method("perform_reaction",
      [](IWReaction& rxn, Molecule& scaffold, const Set_of_Atoms* embedding, Molecule& product)->bool{
        return rxn.perform_reaction(&scaffold, embedding, product);
      },
      "perform reaction"
    )
#ifdef CANNOT_FIGURE_OUT_OPTIONAL_IN_CXXWRAP
    .method("perform_reaction",
      [](IWReaction& rxn, Molecule& scaffold, const Set_of_Atoms* embedding)->std::optional<Molecule>{
        Molecule product;
        if (! rxn.perform_reaction(&scaffold, embedding, product)) {
          return std::nullopt;
        }
        return product;
      },
      "Perform reaction with particular set of matched atoms"
    )
    .method("perform_reaction",
      [](IWReaction& rxn, const Molecule& scaffold, const Set_of_Atoms* embedding,
         const Reaction_Iterator& iter)->std::optional<Molecule>{
        Molecule product;
        if (! rxn.perform_reaction(&scaffold, embedding, iter, product)) {
          return std::nullopt;
        }
        return product;
      },
      "generate product based on embedding and iter"
    )
#endif
#ifdef ITER_LATER
#endif
    ;
      

  mod.method("initialise",
    [](Reaction_Iterator& rit, IWReaction& rxn) {
      return rit.initialise(rxn);
    }
  );
  mod.method("initialise",
    [](Reaction_Iterator& rit, const jlcxx::BoxedValue<IWReaction>& boxed_reaction)->void {
      IWReaction& rxn = jlcxx::unbox<IWReaction&>(boxed_reaction);
      rit.initialise(rxn);
    }
  );
  mod.method("increment!",
    [](Reaction_Iterator& iter) {
      iter++;
    },
    "maybe we could use ++"
  );

  mod.method("perform_reaction",
    [](IWReaction& rxn, const Molecule& scaffold, const Set_of_Atoms* embedding,
       const Reaction_Iterator& iter, Molecule& product)->bool{
      return rxn.perform_reaction(&scaffold, embedding, iter, product);
    },
    "generate product based on embedding and iter"
  );
}


}  // namespace lillymol_julia
