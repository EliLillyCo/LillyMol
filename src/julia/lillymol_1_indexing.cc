#include <string>
#include <vector>

#include "jlcxx/jlcxx.hpp"

#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"

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
  types.add_bits<FileType>("FileType", jlcxx::julia_type("CppEnum"));
  types.set_const("SMI", FILE_TYPE_SMI);
  types.set_const("SDF", FILE_TYPE_SDF);

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

void
IncrementAtomNumbers(Bond& b) {
  b.set_a1(b.a1() + 1);
  b.set_a2(b.a2() + 1);
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
  mod.set_const("SINGLE", kSingleBond);
  mod.set_const("DOUBLE", kDoubleBond);
  mod.set_const("TRIPLE", kTripleBond);
  mod.set_const("AROMATIC", kAromaticBond);
    
  mod.add_bits<FileType>("FileType", jlcxx::julia_type("CppEnum"));
  mod.set_const("SMI", FILE_TYPE_SMI);
  mod.set_const("SDF", FILE_TYPE_SDF);

  mod.add_type<Set_of_Atoms>("SetOfAtoms")
    .constructor<>()
    .method("contains", &Set_of_Atoms::contains)
  ;

  mod.add_type<Ring>("Ring")
    .constructor<>()
    .method("atoms_in_ring",
      [](const Ring& r) {
        return r.number_elements();
      }
    )
    .method("ring_number", &Ring::ring_number)
    .method("fragment_membership", &Ring::fragment_membership)
    .method("fused_system_identifier", &Ring::fused_system_identifier)
    .method("fused_ring_neighbours", &Ring::fused_ring_neighbours)
    .method("fused_neighbour", &Ring::fused_neighbour)
    .method("largest_number_of_bonds_shared_with_another_ring", &Ring::largest_number_of_bonds_shared_with_another_ring)
    .method("strongly_fused_ring_neighbours", &Ring::strongly_fused_ring_neighbours)
    .method("contains_bond", &Ring::contains_bond)
    .method("contains_both", &Ring::contains_both)

    .method("is_fused",
      [](const Ring& r)->bool{
        return r.is_fused();
      }
    )
    .method("is_aromatic",
      [](const Ring& r)->bool{
        return r.is_aromatic();
      }
    )
    .method("contains",
      [](const Set_of_Atoms& s, atom_number_t a)->bool{
        return s.contains(a - 1);
      }
    )
  ;

  mod.set_override_module(jl_base_module);
  mod.method("getindex",
    [](const Ring& a, int i)->atom_number_t{
      return a[i-1];
    }
  );
  mod.method("getindex",
    [](const jlcxx::BoxedValue<Ring>& boxed_ring, int ndx)->atom_number_t{
      const Ring& r = jlcxx::unbox<Ring&>(boxed_ring);
      return r[ndx - 1];
    }
  );
#ifdef RING_LENGTH
  mod.method("length",
    [](const Ring& r){
      return r.number_elements();
    }
  );
  mod.method("length",
    [](const jlcxx::BoxedValue<Ring>& boxed_ring){
      const Ring& r = jlcxx::unbox<Ring&>(boxed_ring);
      return r.number_elements();
    }
  );
#endif
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
#endif
  mod.method("in",
    [](const jlcxx::BoxedValue<Set_of_Atoms>& boxed_set_of_atoms, atom_number_t atom){
      const Set_of_Atoms& s = jlcxx::unbox<Set_of_Atoms&>(boxed_set_of_atoms);
      return s.contains(atom);
    }
  );
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
    .method("process", 
      [](Chemical_Standardisation& s, jlcxx::BoxedValue<Molecule>& boxed_mol)->int {
        Molecule& m = jlcxx::unbox<Molecule&>(boxed_mol);
        return s.process(m);
      }
    )
  ;
    
  // Atom numbers from Bond objects are not offset, since Bonds have only come from
  // Molecules and will have been adjusted when given to Julia
  mod.add_type<Bond>("Bond")
    .method("a1",
      [](const Bond& b) {
        return b.a1();
      }
    )
    .method("a1",
      [](const jlcxx::BoxedValue<const Bond>& boxed_bond) {
        const Bond& b = jlcxx::unbox<const Bond&>(boxed_bond);
        return b.a1();
      }
    )
    .method("a2",
      [](const Bond& b) {
        return b.a2();
      }
    )
    .method("a2",
      [](const jlcxx::BoxedValue<const Bond>& boxed_bond) {
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
  ;

  mod.add_type<Atom>("Atom")
    // .constructor<jlcxx::cxxint_t>(true) // with finalizer
    .method("valence_ok", &Atom::valence_ok)
    .method("ncon",
      [](const Atom& a) {
        return a.ncon();
      }
    )
    .method("isotope", &Atom::isotope)
    .method("lone_pair_count",
      [](Atom&a){
        int lp;
        if (a.lone_pair_count(lp)) {
          return lp;
        }
        return 0;
      }
    )
    .method("nbonds",
      [](const Atom& a) {
        return a.nbonds();
      }
    )
    .method("nbonds",
      [](Atom& a) {
        return a.nbonds();
      }
    )
    .method("bond_to_atom",
      [](const Atom& a, atom_number_t o)->Bond{
        Bond result = *a.bond_to_atom(o);
        IncrementAtomNumbers(result);
        return result;
      }
    )
    .method("is_bonded_to",
      [](const Atom& a, atom_number_t o)->bool{
        return a.is_bonded_to(o);
      }
    )
    .method("is_halogen", &Atom::is_halogen)
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
        Set_of_Atoms result = a.connections(me);
        result.EachAtomIncrement(1);
        return result;
      }
    )
    .method("atomic_weight", &Atom::atomic_weight)
    .method("saturated", 
      [](const Atom& a)->bool{
        return a.fully_saturated();
      }
    );

    mod.set_override_module(jl_base_module);
    mod.method("getindex",
      [](const Atom& a, int i)->Bond{
        Bond result = *a[i-1];
        IncrementAtomNumbers(result);
        return result;
      }
    );
    mod.unset_override_module();
  ;

  mod.add_type<Bond_list>("BondList")
  ;

  mod.add_type<Molecule>("Molecule")
    .constructor<>()
    //.constructor<jlcxx::cxxint_t>(false) // no finalizer
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
        return m.x(a-1);
      }
    )
    .method("y",
      [](const Molecule& m, atom_number_t a){
        return m.y(a-1);
      }
    )
    .method("z",
      [](const Molecule& m, atom_number_t a){
        return m.z(a-1);
      }
    )
    .method("set_x",
      [](Molecule& m, atom_number_t a, float x){
        return m.setx(a-1, x);
      }
    )
    .method("set_y",
      [](Molecule& m, atom_number_t a, float y){
        return m.sety(a-1, y);
      }
    )
    .method("set_z",
      [](Molecule& m, atom_number_t a, float z){
        return m.setz(a-1, z);
      }
    )

    .method("add_bond",
      [](Molecule& m, atom_number_t a1, atom_number_t a2, BondType bt)->bool{
        return m.add_bond(a1 - 1, a2 - 1, BtypeEnumToBtype(bt));
      }
    )
    .method("are_bonded",
      [](const Molecule& m, atom_number_t a1, atom_number_t a2)->bool{
        return m.are_bonded(a1 - 1, a2 - 1);
      }
    )
    .method("has_formal_charges",
      [](const Molecule& m)->bool {
        return m.has_formal_charges();
      }
    )
    .method("formal_charge",
      [](const Molecule& m, atom_number_t a){
        return m.formal_charge(a - 1);
      }
    )
    .method("isotope",
      [](const Molecule& m, atom_number_t a){
        return m.isotope(a - 1);
      }
    )
    .method("number_formally_charged_atoms", &Molecule::number_formally_charged_atoms)
    .method("net_formal_charge", &Molecule::net_formal_charge)

    .method("set_name",
      [](Molecule& m, const std::string& s) {
        return m.set_name(s);
      }
    )
    .method("name",
      [](const Molecule& m)->std::string{
        return m.name().AsString();
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
      [](const Molecule& m)->int{
        return m.natoms();
      }
    )
    .method("natoms",
      [](const Molecule& m, atomic_number_t z)->int {
        return m.natoms(z);
      }
    )
    .method("natoms",
      [](const Molecule& m, std::string& asymbol)->int{
        return m.natoms(asymbol.c_str());
      }
    )
    .method("nedges", &Molecule::nedges)
    .method("atomic_number",
      [](const Molecule& m, atom_number_t a) {
        return m.atomic_number(a - 1);
      }
    )
    .method("nrings",
      [](Molecule& m) {
        return m.nrings();
      }
    )
    .method("nrings",
      [](Molecule& m, atom_number_t a) {
        return m.nrings(a - 1);
      }
    )
    .method("nrings",
      [](Molecule& m, atom_number_t a, int rsize) {
        return m.nrings(a - 1, rsize);
      }
    )
    .method("is_ring_atom",
      [](Molecule& m, atom_number_t a)->bool{
        return m.is_ring_atom(a - 1);
      }
    )
    .method("ring_bond_count", 
      [](Molecule& m, atom_number_t a) {
        return m.ring_bond_count(a - 1);
      }
    )
    .method("fused_system_size",
      [](Molecule& m, atom_number_t a){
        return m.fused_system_size(a - 1);
      }
    )
    .method("rings_with_fused_system_identifier", &Molecule::rings_with_fused_system_identifier)
    .method("fused_system_identifier",
      [](Molecule& m, atom_number_t a) {
        return m.fused_system_identifier(a - 1);
      }
    )
    .method("in_same_ring",
      [](Molecule& m, atom_number_t a1, atom_number_t a2)->bool{
        return m.in_same_ring(a1 - 1, a2 - 1);
      }
    )
    .method("in_same_aromatic_ring",
      [](Molecule& m, atom_number_t a1, atom_number_t a2)->bool{
        return m.in_same_aromatic_ring(a1 - 1, a2 - 1);
      }
    )
    .method("in_same_ring_system",
      [](Molecule& m, atom_number_t a1, atom_number_t a2)->bool{
        return m.in_same_ring_system(a1 - 1, a2 - 1);
      }
    )
    .method("ring_membership",
      [](Molecule& m)->std::vector<int>{
        const int* r = m.ring_membership();
        return std::vector<int>(r, r + m.natoms());
      }
    )
    .method("rings_containing_both",
      [](Molecule& m, atom_number_t a1, atom_number_t a2) {
        return m.in_same_rings(a1 - 1, a2 - 1);
      }
    )
    .method("is_part_of_fused_ring_system",
      [](Molecule& m, atom_number_t a)->bool{
        return m.is_part_of_fused_ring_system(a - 1);
      }
    )
    .method("ring",
      [](Molecule& m, int rnum)->Ring{
        Ring r(*m.ringi(rnum - 1));
        r.EachAtomIncrement(1);
        return r;
      }
    )
    .method("ring_containing_atom",
      [](Molecule& m, atom_number_t a)->Ring{
        const Ring* r = m.ring_containing_atom(a - 1);
        if (r == nullptr) {
          return Ring();
        }
        Ring result(*r);
        result.EachAtomIncrement(1);
        return result;
      }
    )

#ifdef NEEDS_TO_BE_FIXED
    // Not working LoadError: No appropriate factory for type St6vectorI4RingSaIS0_EE
    .method("sssr_rings",
      [](Molecule& m)->std::vector<Ring>{
        std::vector<Ring> result;
        result.reserve(m.nrings());
        for (const Ring* r : m.sssr_rings()) {
          Ring tmp(*r);
          tmp.EachAtomIncrement(1);
          result.push_back(tmp);
        }
        return result;
      }
    )
#endif // NEEDS_TO_BE_FIXED

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
        return m.nrings_including_non_sssr_rings(a - 1);
      }
    )
    .method("non_sssr_rings", &Molecule::non_sssr_rings)
    .method("non_sssr_ring",
      [](Molecule& m, int rnum)->Ring{
        Ring result(*m.non_sssr_ring(rnum - 1));
        result.EachAtomIncrement(1);
        return result;
      }
    )
    .method("is_spiro_fused",
      [](Molecule& m, atom_number_t a)->bool{
        return m.is_spiro_fused(a - 1);
      }
    )
    .method("is_halogen", 
      [](const Molecule& m, atom_number_t a)->bool{
        return m.is_halogen(a-1);
      }
    )
    .method("ncon",
      [](const Molecule& m, atom_number_t a) {
        return m.ncon(a - 1);
      }
    )
    .method("nbonds",
      [](const Molecule& m, atom_number_t a) {
        return m.nbonds(a - 1);
      }
    )
    .method("maximum_connectivity", &Molecule::maximum_connectivity)
    .method("other",
      [](Molecule& m, atom_number_t atom, int ndx){
        return m.other(atom - 1, ndx) + 1;
      }
    )
    .method("connections",
      [](const Molecule& m, atom_number_t a)->std::vector<atom_number_t>{
        std::vector<atom_number_t> result;
        for (atom_number_t o : m.connections(a - 1)) {
          result.push_back(o + 1);
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
    .method("is_aromatic", 
      [](Molecule& m, atom_number_t a)->bool{
        return m.is_aromatic(a - 1);
      }
    )
    .method("atom", &Molecule::atom)

    .method("change_to_graph_form",
      [](Molecule& m) {
        return m.change_to_graph_form();
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
          result.push_back(order[i] + 1);
        }
        return result;
      }
    )
    .method("bond", 
      [](const Molecule& m, int ndx) {
        Bond result = *m.bondi(ndx - 1);
        IncrementAtomNumbers(result);
        return result;
      }
    )
    .method("bond_between_atoms",
      [](Molecule& m, atom_number_t a1, atom_number_t a2){
        Bond result = *m.bond_between_atoms(a1 - 1, a2 - 1);
        IncrementAtomNumbers(result);
        return result;
      }
    )

    //.method("compute_canonical_ranking", &Molecule::compute_canonical_ranking)
    .method("canonical_rank",
      [](Molecule& m, atom_number_t a) {
        return m.canonical_rank(a - 1);
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
        return m.symmetry_class(a - 1);
      }
    )
    .method("number_symmetry_classes", &Molecule::number_symmetry_classes)
    .method("symmetry_equivalents",
      [](Molecule& m, atom_number_t a)->Set_of_Atoms {
        Set_of_Atoms result;
        m.symmetry_equivalents(a - 1, result);
        result.EachAtomIncrement(1);
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
        return m.attached_heteroatom_count(a-1);
      }
    )
    //.method("multiple_bond_to_heteroatom", &Molecule::multiple_bond_to_heteroatom)

    .method("bond_length",
      [](const Molecule& m, atom_number_t a1, atom_number_t a2){
        return m.bond_length(a1 - 1, a2 - 1);
      }
    )
    .method("bond_angle",
      [](const Molecule& m, atom_number_t a1, atom_number_t a2, atom_number_t a3){
        return m.bond_angle(a1 - 1, a2 - 1, a3 - 1);
      }
    )
    .method("dihedral_angle",
      [](const Molecule& m, atom_number_t a1, atom_number_t a2, atom_number_t a3, atom_number_t a4){
        return m.dihedral_angle(a1 - 1, a2 - 1, a3 - 1, a4 - 1);
      }
    )
    .method("signed_dihedral_angle",
      [](const Molecule& m, atom_number_t a1, atom_number_t a2, atom_number_t a3, atom_number_t a4){
        return m.signed_dihedral_angle(a1 - 1, a2 - 1, a3 - 1, a4 - 1);
      }
    )
    .method("set_bond_length",
      [](Molecule& m, atom_number_t a1, atom_number_t a2, float dist){
        return m.set_bond_length(a1 - 1, a2 - 1, dist);
      }
    )
    .method("set_bond_angle",
      [](Molecule& m, atom_number_t a1, atom_number_t a2, atom_number_t a3, angle_t angle){
        return m.set_bond_angle(a1 - 1, a2 - 1, a3 - 1, angle);
      }
    )
    .method("set_dihedral_angle",
      [](Molecule& m, atom_number_t a1, atom_number_t a2, atom_number_t a3, atom_number_t a4, angle_t angle){
        return m.set_dihedral(a1 - 1, a2 - 1, a3 - 1, a4 - 1, angle);
      }
    )
    .method("bond_list",
      [](const Molecule& m)->Bond_list{
        const Bond_list& bonds = m.bond_list();
        Bond_list result;
        for (const Bond* b : m.bond_list()) {
          Bond* julia_bond = new Bond(*b);
          IncrementAtomNumbers(*julia_bond);
          result << julia_bond;
        }
        return result;
      }
    )
      
  ;

  mod.set_override_module(jl_base_module);
  mod.method("getindex",
    [](const Molecule& m, int i)->const Atom&{
      return m[i-1];
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
}

}  // namespace lillymol_julia
