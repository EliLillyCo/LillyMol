module LillyMol
  using CxxWrap

  # abstract type AbstractSetOfAtoms <: AbstractVector{Int32} end

  import Base: getindex, iterate, in, length, size, push!, show, occursin, contains, ==, eachindex

  so_file = joinpath("bazel-bin/julia/", "lillymol_julia.so")
  #@wrapmodule(() -> joinpath("bazel-bin/julia/", "lillymol_julia.so"))
  @wrapmodule(() -> so_file)

  function __init__()
    @initcxx
  end
  export World
  export greet

  export BondType, SINGLE_BOND, DOUBLE_BOND, TRIPLE_BOND, AROMATIC_BOND
  export FileType, SMI, SDF
  export RIP_NONE, RIP_AROMATIC, RIP_FUSED
  export next_molecule, molecules_read, set_connection_table_errors_allowed, connection_table_errors_encountered
  export set_skip_first, read_only, molecules_remaining, molecules_read

  export Molecule, SetOfAtoms, Atom, Bond, ChemicalStandardisation, BondList, Mol2Graph, ChiralCentre
  export SetOfRings
  # export RingAtoms
  export SetOfChiralCentres
  export RingInformation
  export Components

  iterate(m::Molecule, state=0) = (state >= natoms(m) ? nothing : (m[state], state + 1))
  iterate(a::Atom, state=1) = (state > ncon(a) ? nothing : (a[state], state + 1))
  iterate(b::Bond, state=1) = (state == 1 ? (a1(b), 2) : state == 2 ? (a2(b), 2) : nothing)
  iterate(b::BondList, state=1) = (state > bonds_in_set(b) ? nothing : (b[state], state + 1))
  iterate(r::Ring, state=1) = (state > length(r) ? nothing : (r[state], state + 1))
  iterate(s::SetOfAtoms, state=1) = (state > length(s) ? nothing : (s[state], state + 1))
  iterate(r::SetOfRings, state=1) = (state > length(r) ? nothing : (r[state], state + 1))
  iterate(r::SetOfChiralCentres, state=1) = (state > length(r) ? nothing : (r[state], state + 1))
  iterate(r::RingInformation, state=1) = (state > length(r) ? nothing : (r[state], state + 1))
  iterate(r::Components, state=1) = (state > length(r) ? nothing : (r[state], state + 1))
  in(z::Int, m::Molecule) = (natoms(m, z) > 0)
  in(atom::Int, a::Atom) = involves(a, atom)
  length(m::Molecule) = natoms(m)
  length(r::Ring) = atoms_in_ring(r)
  length(s::SetOfRings) = rings_in_set(s)
  length(s::SetOfChiralCentres) = items_in_set(s)
  length(b::BondList) = bonds_in_set(b)
  #length(r::RingAtoms) = nrings(r)
  length(r::RingInformation) = nrings(r)
  # length(r::Ring) = size(r)
  size(m::Molecule) = natoms(m)
  size(r::Ring) = (atoms_in_ring(r),)
  size(s::SetOfAtoms) = (length(s),)
  size(s::SetOfAtomsAllocated) = (length(s),)
  size(s::SetOfRings) = (rings_in_set(s),)
  size(s::SetOfChiralCentres) = (items_in_set(s),)
  size(b::BondList) = (length(b),)
  #size(r::RingAtoms) = (length(r),)
  eachindex(m::Molecule) = 0:(length(m) - 1)
  eachindex(s::SetOfEmbeddings) = 1:(length(s) - 1)
  eachindex(s::CxxWrap.CxxWrapCore.ConstCxxRef{SetOfEmbeddings}) = 1:(length(s) - 1)
  export getindex
  export iterate
  export length
  export in
  export eachindex
  export distance_between_atoms
  export atoms_in_ring, contains
  export is_fused, is_fused_to, largest_number_of_bonds_shared_with_another_ring
  export natoms, smiles, unique_smiles, nrings, ring_bond_count, build_from_smiles, is_aromatic, set_name!, name
  export aromatic_smiles
  export get_coordinates, set_xyz!, dihedral_scan
  export in_ring_of_given_size
  export random_smiles, smiles_starting_with_atom
  export atomic_number, molecular_formula, nedges, is_ring_atom, fused_system_size, fused_system_identifier
  export rings_with_fused_system_identifier, in_same_ring, in_same_aromatic_ring, in_same_ring_system
  export Ring
  export gather_rings
  export ring_membership, rings_containing_both, is_part_of_fused_ring_system, ring, ring_containing_atom
  export label_atoms_by_ring_system, label_atoms_by_ring_system_including_spiro_fused, number_ring_systems
  export nrings_including_non_sssr_rings, non_sssr_rings, non_sssr_ring, is_spiro_fused, is_halogen
  export maximum_connectivity, connections, isotopically_labelled_smiles, is_aromatic, atom, formal_charge
  export number_formal_charges, has_formal_charges
  export set_formal_charge!
  export set_bond_type_between_atoms!, set_atomic_number!
  export isotope, set_isotope!, set_isotopes!, number_isotopic_atoms, remove_isotopes!
  export number_formally_charged_atoms, net_formal_charge, bond, bond_between_atoms
  export compute_aromaticity_if_needed, bond_between_atoms, number_symmetry_classes, symmetry_class, symmetry_equivalents
  export symmetry_classes, attached_heteroatom_count, add_bond!, are_bonded, bond_between_atoms
  export rings, sssr_rings, rings_in_set, bond_list, bonds, bonds_in_set
  export add!, add_atom!, remove_atom!, remove_atoms!, delete_fragment!, chop!
  export set_copy_name_in_molecule_copy_constructor
  export remove_fragment_containing_atom!, remove_all!, atomic_symbol, remove_all_non_natural_elements!
  export remove_explicit_hydrogens!, valence_ok, remove_bonds_to_atom!, remove_bond!, remove_bond_between_atoms!
  export remove_all_bonds!, molecular_weight, amw, molecular_weight_count_isotopes, molecular_weight_ignore_isotopes
  export bond_length, bond_angle, dihedral_angle, signed_dihedral_angle, highest_coordinate_dimensionality, exact_mass
  export translate, discern_chirality_from_3d_structure
  export centre, top_front, top_back, left_down, right_down
  export is_chiral_implicit_hydrogen
  export number_chiral_centres, chiral_centres, chiral_centre_at_atom, chiral_centre_in_molecule_not_indexed_by_atom_number
  export remove_chiral_centre_at_atom!, remove_all_chiral_centres!, invert_chirality_on_atom!
  export number_fragments, fragment_membership, atoms_in_fragment, get_atoms_in_fragment, largest_fragment
  export identify_spinach, rings_in_fragment, create_components, create_subset
  export reduce_to_largest_fragment!, reduce_to_largest_organic_fragment!, reduce_to_largest_fragment_carefully!
  export organic_only, contains_non_periodic_table_elements
  export longest_path, atoms_between, bonds_between, most_distant_pair
  export implicit_hydrogens, explicit_hydrogens, hcount, move_hydrogens_to_end_of_connection_table!
  export unset_all_implicit_hydrogen_information, remove_hydrogens_known_flag_to_fix_valence_errors
  export make_implicit_hydrogens_explicit!, pi_electrons, lone_pair_count, saturated, unsaturation
  export aromatic_atom_count, aromatic_ring_count, unset_unnecessary_implicit_hydrogens_known_values!
  export smarts_equivalent_for_atom, smarts
  export atom_map_number, set_atom_map_number!, atom_with_atom_map_number, reset_all_atom_map_numbers!, reset_atom_map_numbers!
  export set_include_atom_map_with_smiles
  export to_scaffold!
  export dihedral_scan
  export x, y, z, set_x, setx!, sety!, setz!, setxyz!
  export ncon, nbonds, involves, other
  export is_single_bond, is_double_bond, is_triple_bond, is_aromatic
  export largest_ring_size
  export atomic_number
  export a1, a2
  export sort_atoms!
  export set_display_smiles_interpretation_error_messages
  show(io::IO, s::SetOfAtoms) = print(io, set_of_atoms_show_text(s))
  show(io::IO, s::Atom) = print(io, atom_show_text(s))
  show(io::IO, s::Bond) = print(io, bond_show_text(s))
  # For some reason this does not work, but works for set_of_atoms.
  show(io::IO, r::Ring) = print(io, ring_show_text(r))
  show(io::IO, r::CxxWrap.CxxWrapCore.ConstCxxPtr{Ring}) = print(io, ring_show_text(r))
  show(io::IO, m::Molecule) = print(io, molecule_show_text(m))
  show(io::IO, m::CxxWrap.CxxWrapCore.CxxPtr{Molecule}) = print(io, molecule_show_text(m))
  ==(s::SetOfAtoms, v::Vector{Int}) = equals(s, v)
  ==(s::CxxWrap.CxxWrapCore.ConstCxxRef{SetOfAtoms}, v::Vector{Int}) = equals(s, v)

  # None of these work. Not sure why.
  ==(s::Ring, v::Vector{Int}) = ring_equals_vector(s, v)
  ==(s::CxxWrap.CxxWrapCore.ConstCxxPtr{<:Ring}, v::Vector{Int}) = ring_equals_vector(s, v)
  # Get rid of this once we figure out the ==(Ring problem...
  export ring_equals_vector

  export xlogp

  export activate_all, process
  # Not sure why this was necessary. The returned BondList from bond_list(m) did not
  # generate a proper getindex method, and defining this method was suggested
  # as the solution.
  getindex(blist::Union{CxxWrap.CxxWrapCore.ConstCxxRef{<:BondList}, CxxWrap.CxxWrapCore.CxxRef{<:BondList}}, ndx::Int64) =
    internal_get_item(blist, ndx)
  getindex(atom::Union{CxxWrap.CxxWrapCore.ConstCxxRef{<:Atom}, CxxWrap.CxxWrapCore.CxxRef{<:Atom}}, ndx::Int64) =
    internal_get_item(atom, ndx);
  # There must be another way of doing this, TODO:ianwatson investigate automatic casting.
  smiles(m::CxxWrap.CxxWrapCore.CxxPtr{Molecule}) = internal_smiles(m)
  rings(m::CxxWrap.CxxWrapCore.CxxPtr{Molecule}) = rings(Base.unsafe_convert(CxxWrap.CxxWrapCore.CxxRef{Molecule}, m))
  nrings(m::CxxWrap.CxxWrapCore.CxxPtr{Molecule}) = nrings(Base.unsafe_convert(CxxWrap.CxxWrapCore.CxxRef{Molecule}, m))
  number_fragments(m::CxxWrap.CxxWrapCore.CxxPtr{Molecule}) = number_fragments(Base.unsafe_convert(CxxWrap.CxxWrapCore.CxxRef{Molecule}, m))
# @cxxdereference bonds(m::Molecule)=jl_bonds(m)
  function bonds(m::CxxWrap.reference_type_union(Molecule))
    m = CxxWrap.dereference_argument(m)
    jl_bonds(m)
  end
# @cxxdereference create_components(m:Molecule) = jl_create_components(m)
  function create_components(m::CxxWrap.reference_type_union(Molecule))
    m = CxxWrap.dereference_argument(m)
    jl_create_components(m)
  end
  function add_bond!(m::CxxWrap.reference_type_union(Molecule), a1::Integer, a2::Integer, btype::BondType)
    m = CxxWrap.dereference_argument(m)
    jl_add_bond!(m, a1, a2, btype)
  end

  # Wanted to place substructure in a separate module, but could never make it work. Revisit...
  export SubstructureQuery
  export MoleculeToMatch
  export SubstructureResults
  export MoleculeToQuerySpecifications
  export SetOfEmbeddings

  export build_from_smarts, build_from_smiles, build_from_molecule, read_proto
  export matches, substructure_search
  export set_find_one_embedding_per_atom!, set_find_unique_embeddings_only!, set_min_matches_to_find!, set_max_matches_to_find!
  export set_perceive_symmetry_equivalent_matches!, set_save_matched_atoms!, max_query_atoms_matched_in_search
  export set_embeddings_can_overlap!
  export number_embeddings, embeddings
  iterate(r::SetOfEmbeddings, state=1) = (state > length(r) ? nothing : (r[state], state + 1))
  show(io::IO, s::SetOfEmbeddings) = print(io, set_of_embeddings_show_text(s))
  export each_embedding_set_vector
  export substructure_search_as_vector
  getindex(s::Union{CxxWrap.CxxWrapCore.ConstCxxRef{<:SetOfAtoms}, CxxWrap.CxxWrapCore.CxxRef{<:SetOfAtoms}}, ndx::Int64) = soa_getindex(s, ndx)
  getindex(s::Union{CxxWrap.CxxWrapCore.ConstCxxRef{<:SetOfEmbeddings}, CxxWrap.CxxWrapCore.CxxRef{<:SetOfEmbeddings}}, ndx::Int64) = set_of_embeddings_getindex(s, ndx)

  in(q::SubstructureQuery, m::Molecule) = Bool(matches(q, m))

  export Reaction
  export SidechainMatchConditions
  export ReactionIterator
  export perform_reaction, in_place_transformation
  export in_place_transformation
  export add_sidechain
  export read_textproto
  export initialise
  export increment!
  export active
end

