push!(LOAD_PATH, @__DIR__)
println(LOAD_PATH)

using LillyMol

function boobar()
  println("boobar done")
  true
end
#getindex(m::LillyMol.Molecule, a::Int64)=LillyMol.atom(m, a)

using .LillyMol
println(pathof(LillyMol))
# Never a great idea to use random with tests...
using Random
using Test

# Call greet and show the result
@show LillyMol.greet()

w = LillyMol.World()
println(LillyMol.greet(w));
LillyMol.set(w, "hello")
println(LillyMol.greet(w))

# Mimic the same functionality from GoogleTest
# Implemented after I had finished this, use for new cases.
function unordered_elements_are_array(v1::Vector, v2::Vector)::Bool
  size(v1) != size(v2) && return false
  return sort(v1) == sort(v2)
end

function test_set_of_atoms()::Bool
  s = SetOfAtoms()
  for i in 1:6
    push!(s, i)
  end
  for i in 1:6
    i in s || return false
  end
  for i in 1:6
    s[i] = 6 - i
  end
  collect(s) == [5, 4, 3, 2, 1, 0] || return false

  true
end

function test_set_of_atoms_equals()::Bool
  s = SetOfAtoms()
  push!(s, 0)
  push!(s, 1)
  push!(s, 2)
  s == [0, 1, 2] || return is_failure("soa not eq")
end

function test_empty_molecule()::Bool
  m = Molecule()
  natoms(m) == 0 || return false
  length(name(m)) == 0 || return false
  return true
end
  
function test_methane()::Bool
  m = Molecule()
  build_from_smiles(m, "C methane") || return false
  name(m) == "methane" || return false
  natoms(m) == 1 || return false
  natoms(m, 6) == 1 || return false
  natoms(m, 7) == 0 || return false
  nedges(m) == 0 || return false
  nrings(m) == 0 || return false
  is_ring_atom(m, 0) && return false
  number_fragments(m) == 1 || return false
  atomic_symbol(m, 0) == "C" || return false
  atomic_number(m, 0) == 6 || return false
  molecular_formula(m) == "CH4" || return false
  smiles(m) == "C" || return false
  unique_smiles(m) == "C" || return false
  nrings(m) == 0 || return false
  hcount(m, 0) == 4 || return false
  implicit_hydrogens(m, 0) == 4 || return false
  explicit_hydrogens(m, 0) == 0 || return false
  isotope(m, 0) == 0 || return false
  valence_ok(m) || return false
  highest_coordinate_dimensionality(m) == 0 || return false
  true
end

function test_build_from_smiles()::Bool
  m = Molecule()
  build_from_smiles(m, "c1ccccc1OC foo")
end

function test_set_name()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  set_name!(m, "foo")
  name(m) == "foo"
end

function test_copy_constructor()::Bool
  m1 = Molecule()
  build_from_smiles(m1, "C") || return false
  m2 = Molecule(m1)
  smiles(m1) == smiles(m2) || return false
  add_atom!(m1, 6)
  smiles(m1) == "C.C" || return false
  smiles(m2) == "C" || return false
  true
end

function test_copy_constructor_with_name()::Bool
  m1 = Molecule()
  set_copy_name_in_molecule_copy_constructor(true)
  build_from_smiles(m1, "C methane") || return false
  m2 = Molecule(m1)
  set_copy_name_in_molecule_copy_constructor(false)
  name(m2) == "methane" || return false
  return true
end

function test_copy_operation()::Bool
  m1 = Molecule()
  build_from_smiles(m1, "C") || return false
  m2 = copy(m1)
  smiles(m2) == "C" || return false
  add_atom!(m2, 6)
  smiles(m1) == "C" || return false
  smiles(m2) == "C.C" || return false
  true
end

function test_hydrogen_related()::Bool
  m = Molecule()
  build_from_smiles(m, "C([H])([H])([H])N([H])C([H])([H])C([H])([H])C([H])(C1=C([H])C(=C([H])C(=C1[H])[H])[H])OC1=C([H])C(=C(C(=C1[H])[H])C(F)(F)F)[H] prozac") || return false
  natoms(m) == (22 + 18) || return false
  natoms(m, 1) == 18 || return false

  atomic_number(m, 0) == 6 || return false
  ncon(m, 0) == 4 || return false
  ncon(m, 1) == 1 || return false
  isapprox(amw(m), 309.33, atol=0.01) || return false
  isapprox(exact_mass(m), 309.13404868, atol=0.000001) || return false

  remove_all!(m, 1)
  isapprox(amw(m), 309.33, atol=0.01) || return false
  isapprox(exact_mass(m), 309.13404868, atol=0.000001) || return false

  make_implicit_hydrogens_explicit!(m)
  natoms(m) == (22 + 18) || return false
  natoms(m, 1) == 18 || return false
  isapprox(amw(m), 309.33, atol=0.01) || return false
  isapprox(exact_mass(m), 309.13404868, atol=0.000001) || return false
end

function is_failure(msg, m)::Bool
  println(msg, " ", smiles(m))
  false
end

function is_failure(msg)::Bool
  println(msg)
  false
end

macro is_failure(msg)
  println(msg)
  :false 
end
  
function test_set_bond_type_between_atoms()::Bool
  m = Molecule()
  build_from_smiles(m, "CC") || return is_failure("Bad smiles", m)
  set_bond_type_between_atoms!(m, 0, 1, DOUBLE_BOND)
  smiles(m) == "C=C" || return is_failure("Double bond added", m)
  set_bond_type_between_atoms!(m, 0, 1, TRIPLE_BOND)
  smiles(m) ==  "C#C" || return is_failure("Triple bond added", m)
  set_bond_type_between_atoms!(m, 0, 1, SINGLE_BOND)
  smiles(m) ==  "CC" || return is_failure("single bond added", m)
  true
end
  
function test_set_atomic_number()::Bool
  m = Molecule()
  build_from_smiles(m, "CC") || return is_failure("Bad smiles", m)
  set_atomic_number!(m, 1, 103)
  smiles(m) == "C[Lr]" || return is_failure("Cannot set Lr", m)
  true
end


function test_isotopes()::Bool
  m = Molecule()
  build_from_smiles(m, "[1CH4]") || return is_failure("cannot build smiles, m")
  valence_ok(m) || return is_failure("invalid valence", m)
  isotope(m, 0) == 1 || return is_failure("Not isotope 1", m)
  number_isotopic_atoms(m) == 1 || return is_failure("not 1 isotopic atom", m)
  set_isotope!(m, 0, 2)
  smiles(m) == "[2CH4]" || return is_failure("Bad isotope 2", m)
  set_isotope!(m, 0, 108)
  smiles(m) == "[108CH4]" || return is_failure("Bad isotope 108", m)
  remove_isotopes!(m)
  smiles(m) == "C" || return is_failure("Not C", m)

  build_from_smiles(m, "C[1C][2C][3C]C") || return is_failure("Bad smiles", m)
  valence_ok(m) && return is_failure("valence should be bad", m)
  natoms(m, 1) == 0 || return is_failure("exph present", m)
  hcount(m, 0) == 3 || return is_failure("should be ch3", m)
  hcount(m, 1) == 0 || return is_failure("h on atom 1", m)
  hcount(m, 2) == 0 || return is_failure("h on atom 2", m)
  hcount(m, 3) == 0 || return is_failure("h on atom 3", m)
  hcount(m, 4) == 3 || return is_failure("should be ch3", m)

  for i in eachindex(m)
    unset_all_implicit_hydrogen_information(m, i)
  end

  valence_ok(m) || return is_failure("valence wrong")
  hcount(m, 1) == 2 || return is_failure("should be ch2", m)
  hcount(m, 2) == 2 || return is_failure("should be ch2", m)
  hcount(m, 3) == 2 || return is_failure("should be ch2", m)

  s = SetOfAtoms([1, 2, 3])
  set_isotopes!(m, s, 5)
  isotope(m, 1) == 5 || return is_failure("Should be iso 5", m)
  isotope(m, 2) == 5 || return is_failure("Should be iso 5", m)
  isotope(m, 3) == 5 || return is_failure("Should be iso 5", m)

  remove_isotopes!(m)
  smiles(m) == "CCCCC" || return is_failure("Should be CCCCC")

  build_from_smiles(m, "C[1C][2C][3C]C") || return is_failure("Cannot build smiles", m)
  valence_ok(m) && return is_failure("Bad valence", m)
  remove_hydrogens_known_flag_to_fix_valence_errors(m)
  valence_ok(m) || return is_failure("remove hknown did not fix valence", m)
end


function test_natoms()::Bool
  m = Molecule()
  for i in 1:20
    smiles = "C" ^ i
    build_from_smiles(m, smiles) || return false
    natoms(m) == i || return false
  end
  true
end

function test_natoms_atomic_number()::Bool
  m = Molecule()
  for i in 1:10
    smiles = "CN" ^ i
    build_from_smiles(m, smiles) || return false
    natoms(m, 6) == i || return false
    natoms(m, 7) == i || return false
  end
  true
end

function test_natoms_element()::Bool
  m = Molecule()
  for i in 1:10
    smi = "CN" ^ i
    build_from_smiles(m, smi) || return false
    natoms(m, "C") == i || return false
    natoms(m, "N") == i || return false
  end
  true
end

function test_length_molecule()::Bool
  m = Molecule()
  for i = 1:10
    smi = "CN"^i
    build_from_smiles(m, smi) || return false
    length(m) == (2 * i) || return false
  end
  true
end

function test_atomic_number()::Bool
  atomic_numbers = [9, 6, 7, 8, 92]
  m = Molecule()
  build_from_smiles(m, "FCNO.[U]")
  for (ndx, atom) in enumerate(m)
    atomic_number(m, ndx - 1) == atomic_numbers[ndx] || return false
    atomic_number(atom) == atomic_numbers[ndx] || return false
  end
  true
end

function test_enumerate_atom()::Bool
  m = Molecule()
  build_from_smiles(m, "CC(=N)C") || return is_failure("Bad smiles")
  a = m[1]
  ncon(a) == 3 || return is_failure("Not 3 connected", m)
  for (ndx, bond) in enumerate(a)
    if ndx == 1
      is_single_bond(bond) || return is_failure("1 not single", m)
    elseif ndx == 2
      is_double_bond(bond) || return is_failure("2 not double", m)
    elseif ndx == 3
      is_single_bond(bond) || return is_failure("3 not single", m)
    else
      return is_failure("Invalid index")
    end
  end
  true
end

function test_nedges()::Bool
  mol_and_edges = Dict{String, Int}(""=>0, "[Fe]"=>0, "II"=>1, "III"=>2, "I1II1"=>3)
  for (smi,edges) in mol_and_edges
    m = Molecule()
    build_from_smiles(m, smi) || return false
    nedges(m) == edges || "$(edges) failed" == ""
  end
  true
end

function test_atomic_number_in():Bool
  m = Molecule()
  for i in 1:10
    smi = "CNOC(F)(F)" ^ i
    build_from_smiles(m, smi) || return false
    for z in [6, 7, 8, 9]
      z in m || return false
    end
  end
  true
end

function test_molecular_formula()::Bool
  smi_to_mf = Dict{String, String}("C"=>"CH4", "CO"=>"CH4O",
        "C.C"=>"CH4.CH4", "c1ccccc1"=>"C6H6")
  for (smi, mf) in smi_to_mf
    m = Molecule()
    build_from_smiles(m, smi) || return false
    mf == molecular_formula(m) || return false
  end
  true
end

function test_valence_ok()::Bool
  m = Molecule()
  build_from_smiles(m, "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F") || return false
  valence_ok(m) || return false
  build_from_smiles(m, "[CH2]-C") || return false
  valence_ok(m) && return false
  build_from_smiles(m, "CC(C)(C)(C)(C)(C)(C)(C)(C)(C)C") || return false
  valence_ok(m) && return false
  true
end

function test_atom_valence_ok()::Bool
  m = Molecule()
  build_from_smiles(m, "C=F") || return is_failure("Bad smiles")
  valence_ok(m) && return is_failure("Valence not bad", m)
  valence_ok(m[0]) || return is_failure("Atom 0 bad", m)
  valence_ok(m[1]) && return is_failure("Atom 1 good", m)
  true
end

function test_nrings()::Bool
  smiles_and_nrings = Dict{String, Int}("C"=>0, "CC"=>0, "CCC"=>0,
        "C1CC1"=>1, "C12CC1C2"=>2, "C1CC1C1CC1"=>2)
  m = Molecule()
  for (smi, nr) in smiles_and_nrings
    build_from_smiles(m, smi) || return false
    nrings(m) == nr || return false
  end
  true
end

function test_nrings_atom()::Bool
  smiles_and_nrings = Dict{String, Array{Int}}("C"=>[0], "CC"=>[0,0], "CCC"=>[0,0,0],
        "C1CC1"=>[1,1,1], "C12CC1C2"=>[2,1,2,1], "C1CC1C1CC1"=>[1,1,1,1,1,1])
  m = Molecule()
  for (smi, rings) in smiles_and_nrings
    build_from_smiles(m, smi) || return false
    for (atom, nr) in enumerate(rings)
      atom = atom - 1
      nrings(m, atom) == nr || return false
    end
  end
  true
end

function test_nrings_size()::Bool
  m = Molecule()
  for rsize in 3:10
    smi = "C1" * "C"^(rsize-1) * "1"
    build_from_smiles(m, smi) || return false
    for i in 1:natoms(m)
      nrings(m, i - 1, rsize) == 1 || return false
    end
  end
  build_from_smiles(m, "C12CC1C2") || return false
  nrings(m, 0, 3) == 2 || return false
  nrings(m, 1, 3) == 1 || return false
  nrings(m, 2, 3) == 2 || return false
  nrings(m, 3, 3) == 1 || return false
  true
end

function test_is_ring_atom()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  is_ring_atom(m, 0) == false || return false
  build_from_smiles(m, "C1CC1C")
  is_ring_atom(m, 0) == true || return false
  is_ring_atom(m, 1) == true || return false
  is_ring_atom(m, 2) == true || return false
  is_ring_atom(m, 3) == false || return false

  true
end

function test_ring_bond_count()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  ring_bond_count(m, 1) == 0 || return false
  build_from_smiles(m, "C1CC1C")
  ring_bond_count(m, 0) == 1 || return false
  ring_bond_count(m, 1) == 1 || return false
  ring_bond_count(m, 2) == 1 || return false
  ring_bond_count(m, 3) == 0 || return false

  build_from_smiles(m, "C12CC2C1") || return false
  ring_bond_count(m, 0) == 2 || return false
  ring_bond_count(m, 1) == 1 || return false
  ring_bond_count(m, 2) == 2 || return false
  ring_bond_count(m, 3) == 1 || return false

  true
end

function test_fused_system_size()::Bool
  m = Molecule();
  build_from_smiles(m, "C") || return false
  fused_system_size(m, 0) == 0 || return false
  build_from_smiles(m, "C1CC1C") || return false
  fused_system_size(m, 0) == 1 || return false
  fused_system_size(m, 3) == 0 || return false

  build_from_smiles(m, "C12CC2C1") || return false
  for i in 1:natoms(m)
    fused_system_size(m, i - 1) == 2 || return false
  end
  true
end

function test_fused_system_identifier()::Bool
  m = Molecule()
  build_from_smiles(m, "CC") || return false
  fused_system_identifier(m, 0) == fused_system_identifier(m, 1) || return false
  build_from_smiles(m, "C1CC1") || return false
  fused_system_identifier(m, 0) == fused_system_identifier(m, 1) || return false
  fused_system_identifier(m, 1) == fused_system_identifier(m, 2) || return false
  build_from_smiles(m, "C12CC2C1") || return false
  for i in 2:natoms(m)
    fused_system_identifier(m, 0) == fused_system_identifier(m, i - 1) || return false
  end
  true
end

function test_rings_with_fused_system_identifier()::Bool
  m = Molecule()
  build_from_smiles(m, "C1CC1")
  fsid = fused_system_identifier(m, 0);
  rings_with_fused_system_identifier(m, fsid) == 1 || return false
  build_from_smiles(m, "C1CC1C1CC1")
  fsid = fused_system_identifier(m, 0);
  rings_with_fused_system_identifier(m, fsid) == 1 || return false
  fsid = fused_system_identifier(m, 3);
  rings_with_fused_system_identifier(m, fsid) == 1 || return false
  build_from_smiles(m, "C12CC2C1") || return false
  fsid = fused_system_identifier(m, 3);
  rings_with_fused_system_identifier(m, fsid) == 2 || return false
  true
end

function test_in_same_ring()::Bool
  m = Molecule();
  build_from_smiles(m, "CC") || return false
  in_same_ring(m, 0, 1) && return false
  build_from_smiles(m, "C1CC1C") || return false
  in_same_ring(m, 0, 1) || return false
  in_same_ring(m, 1, 2) || return false
  in_same_ring(m, 0, 2) || return false
  ! in_same_ring(m, 0, 3) || return false
  ! in_same_ring(m, 1, 3) || return false
  ! in_same_ring(m, 2, 3) || return false

  build_from_smiles(m, "C12CC2C1") || return false
  in_same_ring(m, 0, 1) || return false
  in_same_ring(m, 0, 2) || return false
  in_same_ring(m, 0, 3) || return false
  in_same_ring(m, 1, 2) || return false
  ! in_same_ring(m, 1, 3) || return false
  in_same_ring(m, 2, 3) || return false
end

function test_in_same_aromatic_ring()::Bool
  m = Molecule()
  build_from_smiles(m, "Oc1cccnc1") || return false
  for i in 1:6
    in_same_aromatic_ring(m, 0, i - 1) && return false
  end
  for i in 1:6
    for j in i:6
      in_same_aromatic_ring(m, i, j) || return false
    end
  end
  build_from_smiles(m, "c12cncc2cccc1") || return false
  in_same_aromatic_ring(m, 1, 2) || return false
  in_same_aromatic_ring(m, 1, 3) || return false
  in_same_aromatic_ring(m, 1, 4) || return false
  ! in_same_aromatic_ring(m, 1, 5) || return false
  ! in_same_aromatic_ring(m, 5, 1) || return false
  ! in_same_aromatic_ring(m, 1, 6) || return false
  ! in_same_aromatic_ring(m, 6, 1) || return false
  ! in_same_aromatic_ring(m, 1, 7) || return false
  ! in_same_aromatic_ring(m, 1, 8) || return false
end

function test_in_same_ring_system()::Bool
  m = Molecule()
  build_from_smiles(m, "C1CC1CC1CC1") || return false
  in_same_ring_system(m, 0, 1) || return false
  in_same_ring_system(m, 0, 2) || return false
  in_same_ring_system(m, 1, 2) || return false
  in_same_ring_system(m, 4, 5) || return false
  in_same_ring_system(m, 4, 6) || return false
  in_same_ring_system(m, 5, 6) || return false
  ! in_same_ring_system(m, 2, 3) || return false
  ! in_same_ring_system(m, 3, 4) || return false
  ! in_same_ring_system(m, 2, 4) || return false
  ! in_same_ring_system(m, 1, 5) || return false
  ! in_same_ring_system(m, 0, 6) || return false
end

function test_ring_membership()::Bool
  m = Molecule();
  build_from_smiles(m, "CCC") || return false
  r = ring_membership(m)
  sum(r) == 0 || return false
  build_from_smiles(m, "C1CC1") || return false
  r = ring_membership(m)
  # https://stackoverflow.com/questions/47564825/check-if-all-the-elements-of-a-julia-array-are-equal#:~:text=The%20shortest%20way%20I%20can,%3D%3D%20arr)%20.
  # all(y->y==r[1], r) || return false
  all(y->y==1, r) || return false
  build_from_smiles(m, "C12CC2C1")
  r = ring_membership(m)
  [2, 1, 2, 1] == r || return false
  true
end

function test_rings_containing_both()::Bool
  m = Molecule();
  build_from_smiles(m, "C1CC1C") || return false
  rings_containing_both(m, 0, 1) == 1 || return false
  rings_containing_both(m, 0, 2) == 1 || return false
  rings_containing_both(m, 1, 2) == 1 || return false
  rings_containing_both(m, 0, 3) == 0 || return false
  rings_containing_both(m, 1, 3) == 0 || return false
  rings_containing_both(m, 2, 3) == 0 || return false
  build_from_smiles(m, "C12CC2C1C") || return false
  rings_containing_both(m, 0, 1) == 1 || return false
  rings_containing_both(m, 0, 2) == 2 || return false
  rings_containing_both(m, 0, 3) == 1 || return false
  rings_containing_both(m, 1, 2) == 1 || return false
  rings_containing_both(m, 2, 3) == 1 || return false
  rings_containing_both(m, 3, 4) == 0 || return false
  true
end

function test_is_part_of_fused_ring_system()::Bool
  m = Molecule()
  build_from_smiles(m, "CCC") || return false
  ! is_part_of_fused_ring_system(m, 0) || return false
  build_from_smiles(m, "C1CC1") || return false
  ! is_part_of_fused_ring_system(m, 0) || return false
  ! is_part_of_fused_ring_system(m, 1) || return false
  build_from_smiles(m, "C12CC2C1CC") || return false
  is_part_of_fused_ring_system(m, 0) || return false
  is_part_of_fused_ring_system(m, 1) || return false
  is_part_of_fused_ring_system(m, 2) || return false
  is_part_of_fused_ring_system(m, 3) || return false
  ! is_part_of_fused_ring_system(m, 4) || return false
end

function test_ring()::Bool
  m = Molecule();
  build_from_smiles(m, "C1CC1C") || return false
  r = ring(m, 0)
  l = length(r)
  length(r) == 3 || return false
  for i in 0:2
    i in r || return false
  end
  # Once I switched to using pointers for Rings, this stopped working.
  #atoms = collect(r)
  #for i in 0:2
  #  i in atoms || return false
  #end
  #3 in atoms && return false

  build_from_smiles(m, "C1CC1CC1CC1") || return false
  for i in 1:nrings(m)
    r = ring(m, i - 1)
    length(r) == 3 || return false
  end
  all([x in ring(m, 0) for x in 0:2]) || return false
  all([x in ring(m, 1) for x in 4:6]) || return false
  true
end

function same_atoms(s1, s2)::Bool
  length(s1) == length(s2) || return false
  set1 = Set{Int}()
  for i in s1
    push!(set1, i)
  end
  for i in s2
    i in set1 || return false
  end
  true
end

function test_rings()::Bool
  m = Molecule()
  build_from_smiles(m, "C1CC1C2CCC2C3CCCC3") || return false
  expected = [[0, 1, 2], [3, 4, 5, 6], [7, 8, 9, 10, 11]]
  for (i, r) in enumerate(sssr_rings(m))
    same_atoms(expected[i], r) || return false
  end
  true
end

function test_ring_containing_atom()::Bool
  m = Molecule();
  build_from_smiles(m, "C1CC1C2CCC2C3CCCC3") || return false
  for i in 0:2
    r = ring_containing_atom(m, i)
    length(r) == 3 || return false
  end
  for i in 3:6
    r = ring_containing_atom(m, i)
    length(r) == 4 || return false
  end
  for i in 7:(natoms(m) - 1)
    r = ring_containing_atom(m, i)
    length(r) == 5 || return false
  end
  true
end

# Fix for zero arrays
function test_label_atoms_by_ring_system()::Bool
  m = LillyMol.MolFromSmiles("C1CC1C2CCC2C3CCCCC3C")
  rsys = label_atoms_by_ring_system(m)
  for i in 1:3
    rsys[i] == 1 || return false
  end
  for i in 4:7
    rsys[i] == 2 || return false
  end
  for i in 8:13
    rsys[i] == 3 || return false
  end
  rsys[14] == 0 || return false

  true
end

# Fix for zero arrays
function test_label_atoms_by_ring_system_including_spiro_fused()::Bool
  m = LillyMol.MolFromSmiles("C1CC12CC2")
  rsys = label_atoms_by_ring_system_including_spiro_fused(m)
  for i in 1:natoms(m)
    rsys[i] == 1 || return false
  end
  true
end

function test_nrings_including_non_sssr_rings()::Bool
  m = LillyMol.MolFromSmiles("C12C3C4C1C5C2C3C45")
  for i in eachindex(m)
    nrings_including_non_sssr_rings(m, i) == 3 || return false
  end
  true
end

function test_non_sssr_rings()::Bool
  m = LillyMol.MolFromSmiles("C12C3C4C1C5C2C3C45")
  non_sssr_rings(m) == 1 || return false
  true
end

function test_non_sssr_ring()::Bool
  m = LillyMol.MolFromSmiles("C12C3C4C1C5C2C3C45")
  # 1 == non_sssr_rings(m) || return false
  r = non_sssr_ring(m, 0)
  length(r) == 4 || return false
  true
end

function test_is_spiro_fused()::Bool
  m = LillyMol.MolFromSmiles("C1CC12CC2")
  is_spiro_fused(m, 0) && return false
  is_spiro_fused(m, 1) && return false
  is_spiro_fused(m, 2) || return false
  is_spiro_fused(m, 3) && return false
  is_spiro_fused(m, 4) && return false
  true
end

function test_is_halogen()::Bool
  m = LillyMol.MolFromSmiles("FC(Cl)(Br)CI")
  is_halogen(m, 0) || return false
  is_halogen(m, 1) && return false
  is_halogen(m, 2) || return false
  is_halogen(m, 3) || return false
  is_halogen(m, 4) && return false
  is_halogen(m, 5) || return false
  true
end

function test_ncon_molecule()::Bool
  m = LillyMol.MolFromSmiles("CCC(C)C(C)(C)C")
  ncon(m, 0) == 1 || return false
  ncon(m, 1) == 2 || return false
  ncon(m, 2) == 3 || return false
  ncon(m, 4) == 4 || return false
  true
end

function test_nbonds_molecule()::Bool
  m = LillyMol.MolFromSmiles("CCC=CC#C")
  nbonds(m, 0) == 1 || return false
  nbonds(m, 1) == 2 || return false
  nbonds(m, 2) == 3 || return false
  nbonds(m, 4) == 4 || return false
  true
end

function test_maximum_connectivity()::Bool
  m = Molecule()
  maximum_connectivity(m) == 0 || return false
  build_from_smiles(m, "C") || return false
  maximum_connectivity(m) == 0 || return false
  build_from_smiles(m, "CC") || return false
  maximum_connectivity(m) == 1 || return false
  build_from_smiles(m, "CCC") || return false
  maximum_connectivity(m) == 2 || return false
  build_from_smiles(m, "CC(C)C") || return false
  maximum_connectivity(m) == 3 || return false
  build_from_smiles(m, "CC(C)(C)C") || return false
  maximum_connectivity(m) == 4 || return false
  true
end

function test_connections_molecule()::Bool
  m = Molecule();
  build_from_smiles(m, "C") || return false
  c = connections(m, 0)
  length(c) == 0 || return false
  build_from_smiles(m, "CC") || return false
  c = connections(m, 0)
  length(c) == 1 || return false
  c[1] == 1 || return false
  c = connections(m, 1)
  length(c) == 1 || return false
  c[1] == 0 || return false
  build_from_smiles(m, "CCC") || return false
  c = connections(m, 1)
  (0 in c && 2 in c) || return false
  build_from_smiles(m, "CC(C)(C)C") || return false
  c = connections(m, 1)
  length(c) == 4 || return false
  (0 in c && 2 in c && 3 in c && 4 in c) || return false
  unordered_elements_are_array(collect(c), [0, 2, 3, 4]) || return false
  true
end

function test_isotopically_labelled_smiles()::Bool
  m = Molecule()
  build_from_smiles(m, "CNO")
  s = isotopically_labelled_smiles(m)
  s == "C[1NH][2OH]" || return false
  true
end

function test_is_aromatic()::Bool
  m = LillyMol.MolFromSmiles("C1CCC1")
  for i in natoms(m)
    is_aromatic(m, i - 1) && return false
  end

  build_from_smiles(m, "c1ccccc1") || return false
  for i in natoms(m)
    is_aromatic(m, i - 1) || return false
  end
  true
end

function test_getindex_molecule()::Bool
  m = LillyMol.MolFromSmiles("[1CH3]NOF")
  atomic_number(m[0]) == 6 || return false
  atomic_number(m, 0) == 6 || return false
  atomic_number(m[1]) == 7 || return false
  atomic_number(m, 1) == 7 || return false
  atomic_number(m[2]) == 8 || return false
  atomic_number(m, 2) == 8 || return false
  atomic_number(m[3]) == 9 || return false
  atomic_number(m, 3) == 9 || return false
  atomic_symbol(m[0]) == "C" || return false
  atomic_symbol(m, 0) == "C" || return false
  ncon(m[0]) == 1 || return false
  ncon(m, 0) == 1 || return false
  ncon(m[1]) == 2 || return false
  ncon(m, 1) == 2 || return false
  ncon(m[2]) == 2 || return false
  ncon(m, 2) == 2 || return false
  ncon(m[3]) == 1 || return false
  ncon(m, 3) == 1 || return false
  formal_charge(m[0]) == 0 || return false
  formal_charge(m, 0) == 0 || return false
  isotope(m[0]) == 1 || return false
  isotope(m, 0) == 1 || return false
  isotope(m[1]) == 0 || return false
  isotope(m, 1) == 0 || return false
  for (ndx, atom) in enumerate(m)
    ndx = ndx - 1
    atomic_number(m[ndx]) == atomic_number(atom) || return false
    atomic_symbol(m[ndx]) == atomic_symbol(atom) || return false
    ncon(m[ndx]) == ncon(atom) || return false
    formal_charge(m[ndx]) == formal_charge(atom) || return false
    isotope(m[ndx]) == isotope(atom) || return false
  end
  true
end

function test_set_formal_charge()::Bool
  m = LillyMol.MolFromSmiles("NCC")
  set_formal_charge!(m, 0, 1)
  formal_charge(m, 0) == 1 || return false
  true
end

function test_number_formally_charged_atoms()::Bool
  m = LillyMol.MolFromSmiles("NC")
  number_formally_charged_atoms(m) == 0 || return false
  build_from_smiles(m, "[NH3+]C")
  number_formally_charged_atoms(m) == 1 || return false
  build_from_smiles(m, "C(=O)[O-]")
  number_formally_charged_atoms(m) == 1 || return false
  true
end

function test_net_formal_charge()::Bool
  m = LillyMol.MolFromSmiles("NC")
  net_formal_charge(m) == 0 || return false
  build_from_smiles(m, "[NH3+]C")
  net_formal_charge(m) == 1 || return false
  build_from_smiles(m, "C(=O)[O-]")
  net_formal_charge(m) == -1 || return false
  true
end

function test_bond_molecule()::Bool
  m = LillyMol.MolFromSmiles("CC=CC#C")
  b = bond(m, 0)
  is_single_bond(b) || return false
  is_double_bond(bond(m, 1)) || return false
  is_triple_bond(bond(m, 3)) || return false
  build_from_smiles(m, "c1ccccc1")
  compute_aromaticity_if_needed(m)
  for i in 1:nedges(m)
    is_aromatic(bond(m, i - 1)) || return false
  end
  true
end

function test_bond_between_atoms()::Bool
  m = LillyMol.MolFromSmiles("CC=C")
  b = bond_between_atoms(m, 0, 1)
  is_single_bond(b) || return false
  b = bond_between_atoms(m, 1, 2)
  is_double_bond(b) || return false
  true
end

function test_number_symmetry_classes()::Bool
  m = LillyMol.MolFromSmiles("C")
  number_symmetry_classes(m) == 1 || return false
  build_from_smiles(m, "CC") || return false
  number_symmetry_classes(m) == 1 || return false
  build_from_smiles(m, "FC(F)(F)C(C)C") || return false
  number_symmetry_classes(m) == 4 || return false
  build_from_smiles(m, "C1CC1") || return false
  number_symmetry_classes(m) == 1 || return false
  true
end

function test_symmetry_class()::Bool
  m = Molecule()
  build_from_smiles(m, "FC(F)(F)C(C)C") || return false
  symmetry_class(m, 0) == symmetry_class(m, 2) || return false
  symmetry_class(m, 0) == symmetry_class(m, 3) || return false
  symmetry_class(m, 5) == symmetry_class(m, 6) || return false
  build_from_smiles(m, "CN")
  symmetry_class(m, 0) == symmetry_class(m, 1) && return false
  true
end

function test_symmetry_equivalents()::Bool
  m = Molecule()
  build_from_smiles(m, "FC(F)(F)C(C)C") || return false
  s = symmetry_equivalents(m, 0)
  length(s) == 2 || return false
  (2 in s && 3 in s) || return false
  s = symmetry_equivalents(m, 6)
  length(s) == 1 || return false
  (6 in s) || return false
  true
end

# Fix for zero arrays
function test_symmetry_classes()::Bool
  m = Molecule()
  build_from_smiles(m, "FC(F)(F)C(C)C") || return false
  c = symmetry_classes(m)
  c[1] == c[3] || return false
  c[1] == c[4] || return false
  c[6] == c[7] || return false
  true
end

function test_attached_heteroatom_count()::Bool
  m = Molecule()
  build_from_smiles(m, "FC(F)(F)C(C)N") || return false
  attached_heteroatom_count(m, 0) == 0 || return false
  attached_heteroatom_count(m, 1) == 3 || return false
  attached_heteroatom_count(m, 2) == 0 || return false
  attached_heteroatom_count(m, 4) == 1 || return false
  attached_heteroatom_count(m, 5) == 0 || return false
  attached_heteroatom_count(m, 6) == 0 || return false
  true
end

function test_bond_length()::Bool
  m = Molecule()
  build_from_smiles(m, "C{{0,0,0}}C{{1,1,1}}") || return false
  @test bond_length(m, 0, 1) ≈ sqrt(3.0) atol=0.001
  true
end

function test_bond_angle()::Bool
  m = Molecule()
  build_from_smiles(m, "C{{0,0,0}}C{{1,0,0}}C{{1,1,1}}") || return false
  @test bond_angle(m, 0, 1, 2) ≈ (π / 2.0) atol=0.001
  true
end

function test_dihedral_angle()::Bool
  m = Molecule()
  build_from_smiles(m, "C{{0,0,0}}C{{1,0,0}}C{{1,1,1}}") || return false
  @test bond_angle(m, 0, 1, 2) ≈ (π / 2.0) atol=0.001
  true
end

function test_add_bond()::Bool
  m = Molecule()
  build_from_smiles(m, "C.C") || return false
  add_bond!(m, 0, 1, SINGLE_BOND)
  smiles(m) == "CC" || return false
  build_from_smiles(m, "C.C") || return false
  add_bond!(m, 0, 1, DOUBLE_BOND)
  smiles(m) == "C=C" || return false
  build_from_smiles(m, "C.C") || return false
  add_bond!(m, 0, 1, TRIPLE_BOND)
  smiles(m) == "C#C" || return false
  true
end

function test_are_bonded()::Bool
  m = Molecule()
  build_from_smiles(m, "C.CC")
  are_bonded(m, 0, 1) && return false
  are_bonded(m, 1, 2) || return false
  true
end

function test_bond_between_atoms()::Bool
  m = Molecule()
  build_from_smiles(m, "CCC") || return false
  b = bond_between_atoms(m, 0, 1)
  ((a1(b) == 0 && a2(b) == 1) || (a1(b) == 1 && a2(b) == 0)) || return false
  true
end

function test_bond_list()::Bool
  m = Molecule()
  build_from_smiles(m, "CC=CC#C") || return false
  for (bnum, bond) in enumerate(bond_list(m))
    if bnum == 1
      (a1(bond) == 0 && a2(bond) == 1) || return false
    elseif bnum == 2
      (a1(bond) == 1 && a2(bond) == 2) || return false
    elseif bnum == 3
      (a1(bond) == 2 && a2(bond) == 3) || return false
    elseif bnum == 4
      (a1(bond) == 3 && a2(bond) == 4) || return false
    end
  end
  for (bnum, bond) in enumerate(bond_list(m))
    if bnum == 1
      is_single_bond(bond) || return false
    elseif bnum == 2
      is_double_bond(bond) || return false
    elseif bnum == 3
      is_single_bond(bond) || return false
    elseif bnum == 4
      is_triple_bond(bond) || return false
    end
  end
  true
end

function test_add_molecule()::Bool
  m1 = LillyMol.MolFromSmiles("C")
  m2 = LillyMol.MolFromSmiles("N")
  add!(m1, m2)
  smiles(m1) == "C.N" || return false
  add!(m1, m2)
  smiles(m1) == "C.N.N" || return false
  true
end

function test_remove_atom()::Bool
  m = LillyMol.MolFromSmiles("CN")
  remove_atom!(m, 0)
  smiles(m) == "N" || return false
  true
end

function test_remove_atom2()::Bool
  m = Molecule()
  build_from_smiles(m, "CCC") || return is_failure("bad smiles", m)
  remove_atom!(m, 1)
  smiles(m) == "C.C" || return is_failure("After atom removal", m)

  add_atom!(m, 7)
  add_bond!(m, 0, 2, SINGLE_BOND)
  add_bond!(m, 1, 2, SINGLE_BOND)
  smiles(m) == "CNC" || return is_failure("After adding atom", m)
  true
end

function test_remove_atoms()::Bool
  m = LillyMol.MolFromSmiles("CNCOF")
  s = SetOfAtoms()
  add!(s, 1)
  add!(s, 3)
  remove_atoms!(m, s)
  smiles(m) == "C.C.F" || return false
  true
end

function test_remove_atoms_vector()::Bool
  m = LillyMol.MolFromSmiles("CNCOF")
  to_remove = [0, 1, 0, 1, 0]
  remove_atoms!(m, to_remove) == 2 || return false
  smiles(m) == "C.C.F" || return false
  true
end

function test_delete_fragment()::Bool
  m = LillyMol.MolFromSmiles("CCC.NNN.O")
  delete_fragment!(m, 1)
  smiles(m) == "CCC.O" || return false
  true
end

function test_delete_fragment2()::Bool
  m = Molecule()
  build_from_smiles(m, "CC.CCC.CC") || return is_failure("Bad smiles", m)
  # Fragment numbering is not guaranteed. Better to programatically determine
  # an atom in the fragment and remove the fragment containing that atom.
  delete_fragment!(m, 1)
  smiles(m) == "CC.CC" || return is_failure("remaining fragment wrong", m)

  # Remove the fragment containing a nitrogen atom

  build_from_smiles(m, "CC.CNC.CC") || return is_failure("Bad smiles", m)
  to_remove = -1
  for (ndx, atom) in enumerate(m)
    if atomic_number(atom) == 7
      to_remove = fragment_membership(m, ndx)
      break
    end
  end

  to_remove >= 0 || return is_failure("No atom to remove", m)
  delete_fragment!(m, to_remove)
  smiles(m) ==  "CC.CC" || return is_failure("Frag remove bad", m)
  true
end


function test_largest_fragment()::Bool
  m = Molecule()
  build_from_smiles(m, "C.CC.CCC") || return is_failure("Cannot build", m)
  natoms(m) == 6 || return is_failure("Bad atom count", m)
  number_fragments(m) == 3 || return is_failure("Not 3 fragments", m)
  reduce_to_largest_fragment!(m)
  smiles(m) == "CCC" || return is_failure("Not CCC largest frag", m)
  true
end
  
function test_largest_fragment_carefully()::Bool
  m = Molecule()
  build_from_smiles(m, "O=N(=O)C1=CC(=C(O)C(=C1)N(=O)=O)N(=O)=O.CCSC(=N)N CHEMBL4531203") || return is_failure("bad smiles", m)
  reduce_to_largest_fragment_carefully!(m)
  smiles(m) == "O=N(=O)C1=CC(=C(O)C(=C1)N(=O)=O)N(=O)=O" || return is_failure("Bad smiles", m)

  build_from_smiles(m, "O(C(=O)C(C1=NC(C)(C)CC2=CC=CC=C12)CC)CC.OC1=C(N(=O)=O)C=C(N(=O)=O)C=C1N(=O)=O CHEMBL1352385") ||
      return is_failure("Bad smiles", m)
  reduce_to_largest_fragment_carefully!(m)
  smiles(m) ==  "O(C(=O)C(C1=NC(C)(C)CC2=CC=CC=C12)CC)CC" || return is_failure("Wrong smiles", m)

  build_from_smiles(m, "N1=C(NC(C)C1)CC1=CC=CC=C1.C1(=CC(=CC(=C1O)N(=O)=O)N(=O)=O)N(=O)=O CHEMBL609421") || 
                return is_failure("Bad smiles", m)
  reduce_to_largest_fragment_carefully!(m)
  smiles(m) == "N1=C(NC(C)C1)CC1=CC=CC=C1" || return is_failure("Wrong smiles", m)

  true
end

function test_remove_fragment_containing_atom()::Bool
  m = LillyMol.MolFromSmiles("C.N.O.S")
  remove_fragment_containing_atom!(m, 1)
  smiles(m) == "C.O.S" || return false
  true
end

function test_remove_all()::Bool
  m = LillyMol.MolFromSmiles("OC(=O)C")
  remove_all!(m, 8) == 2 || return false
  smiles(m) == "CC" || return false
  true
end

function test_set_auto_create_new_elements()::Bool
  LillyMol.set_auto_create_new_elements(1)
  m = LillyMol.MolFromSmiles("[Th][Eq][U]IC[K][Br]O[W]NFO[Xj][Um]PSO[Ve][Rt][La][Zy][D]O[G]")
  length(m) == 24 || return false
  true
end

function test_set_atomic_symbols_can_have_arbitrary_length()::Bool
  LillyMol.set_atomic_symbols_can_have_arbitrary_length(1)
  m = LillyMol.MolFromSmiles("[Ala][Gly][Ser][Hello]")
  natoms(m) == 4 || return false
  atomic_symbol(m, 0) == "Ala" || return false
  atomic_number(m, 1) < 0 || return false
  LillyMol.set_atomic_symbols_can_have_arbitrary_length(0)
  true
end
  
function test_set_atomic_symbols_can_have_arbitrary_length2()::Bool
  m = Molecule()
  LillyMol.set_atomic_symbols_can_have_arbitrary_length(true)
  LillyMol.set_auto_create_new_elements(true)
  build_from_smiles(m, "[Ala][Glu]") || return is_failure("Bad smiles", m)
  natoms(m) == 2 || return is_failure("Not 2 atoms", m)
  atomic_symbol(m, 0) == "Ala" || return is_failure("Not Ala", m)
  atomic_symbol(m, 1) == "Glu" || return is_failure("Not Glu", m)
  occursin("Ala", m) || return is_failure("Ala not occursin", m)
  contains(m, "Glu") || return is_failure("mol not contains Glu", m) 
  LillyMol.set_auto_create_new_elements(false)
  LillyMol.set_atomic_symbols_can_have_arbitrary_length(false)
  true
end

function test_fused_to()::Bool
  m = Molecule()
  build_from_smiles(m, "C1C2CC12") || return is_failure("Bad smiles", m)
  nrings(m) == 2 || return is_failure("Not 2 rings", m)
  length(ring(m, 0)) == 3 || return is_failure("Ring 0 not 3", m)
  length(ring(m, 1)) == 3 || return is_failure("Ring 1 not 3", m)
  r0 = ring(m, 0)
  r1 = ring(m, 1)
  fused_system_identifier(r0) == fused_system_identifier(r1) || return is_failure("fsid non match", m)
  is_fused(r0) || return is_failure("r0 is_fused", m)
  is_fused(r1) || return is_failure("r1 is_fused", m)
  is_fused_to(r0, r1) || return is_failure("is_fused_to r0 r1", m)
  is_fused_to(r1, r0) || return is_failure("is_fused_to r1 r0", m)

  fragment_membership(r0) == fragment_membership(r1) || return is_failure("frag mismatch", m)

  largest_number_of_bonds_shared_with_another_ring(r0) == 1 || return is_failure("bonds shared 0", m)
  largest_number_of_bonds_shared_with_another_ring(r1) == 1 || return is_failure("bonds shared 1", m)
end

function test_remove_all_non_natural_elements()::Bool
  LillyMol.set_auto_create_new_elements(1)
  m = LillyMol.MolFromSmiles("[Xx]C[Yy]")
  remove_all_non_natural_elements!(m) == 2 || return false
  smiles(m) == "C" || return false
  LillyMol.set_auto_create_new_elements(0)
  true
end

function test_remove_explicit_hydrogens()::Bool
  m = LillyMol.MolFromSmiles("[H]C([H])([H])CC")
  valence_ok(m) || return false
  remove_explicit_hydrogens!(m) == 3 || return false
  true
end

function test_chop()::Bool
  m = LillyMol.MolFromSmiles("CN")
  chop!(m, 1)
  smiles(m) == "C" || return false
  build_from_smiles(m, "CNN") || return false
  chop!(m, 2)
  smiles(m) == "C" || return false
  true
end

function test_remove_bonds_to_atom()::Bool
  m = LillyMol.MolFromSmiles("FCN")
  remove_bonds_to_atom!(m, 1)
  smiles(m) == "F.C.N" || return false
  true
end

function test_remove_bond()::Bool
  m = LillyMol.MolFromSmiles("OCN")
  remove_bond!(m, 0)
  smiles(m) == "O.CN" || return false
  true
end

function test_remove_bond_between_atoms()::Bool
  m = LillyMol.MolFromSmiles("OCN")
  remove_bond_between_atoms!(m, 0, 1)
  smiles(m) == "O.CN" || return false
  true
end

function test_remove_all_bonds()::Bool
  m = LillyMol.MolFromSmiles("CCC")
  remove_all_bonds!(m)
  smiles(m) = "C.C.C" || return false
  true
end

function test_molecular_weight()::Bool
  m = LillyMol.MolFromSmiles("FC(Cl)(Br)CNC(=O)O")
  abs(220.4245 - molecular_weight(m)) < 0.001 || return false
  true
end

function test_molecular_weight_count_isotopes()::Bool
  m = LillyMol.MolFromSmiles("[1C]C")
  abs(molecular_weight_count_isotopes(m) - 16.03452) < 0.001 || return false
  true
end

function test_molecular_weight_ignore_isotopes()::Bool
  m = LillyMol.MolFromSmiles("[1C]C")
  abs(molecular_weight_ignore_isotopes(m) - 27.045221) < 0.001 || return false
  true
end

function test_highest_coordinate_dimensionality()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  highest_coordinate_dimensionality(m) == 0 || return false
  build_from_smiles(m, "C{{1,0,0}}")
  highest_coordinate_dimensionality(m) == 1 || return false
  build_from_smiles(m, "C{{1,1,0}}")
  highest_coordinate_dimensionality(m) == 2 || return false
  build_from_smiles(m, "C{{1,1,1}}")
  highest_coordinate_dimensionality(m) == 3 || return false
  true
end

function test_exact_mass()::Bool
  m = LillyMol.MolFromSmiles("CNOF")
  abs(exact_mass(m) - 65.027691992) < 0.001 || return false
  true
end

function test_number_fragments()::Bool
  m = Molecule()
  for i in 1:10
    smi = "C."^i
    smi *= "C"
    build_from_smiles(m, smi) || return false
    number_fragments(m) == (i + 1) || return false
  end
  true
end

function test_fragment_membership()::Bool
  m = Molecule()
  build_from_smiles(m, "C.CC.CCC")
  fragment_membership(m, 0) == 0 || return false
  fragment_membership(m, 1) == 1 || return false
  fragment_membership(m, 2) == 1 || return false
  fragment_membership(m, 3) == 2 || return false
  fragment_membership(m, 4) == 2 || return false
  fragment_membership(m, 5) == 2 || return false
  true
end
  
function test_fragment_membership2()::Bool
  m = Molecule()
  build_from_smiles(m, "C.CC.CCC.CCCC") || return is_failure("Bad smiles", m)
  fragment_membership(m) == [0, 1, 1, 2, 2, 2, 3, 3, 3, 3] || return is_failure("Bad frag membership", m)
  true
end

function test_fragment_membership_vector()::Bool
  m = Molecule()
  build_from_smiles(m, "C.CC.CCC")
  expected = [0, 1, 1, 2, 2, 2]
  fragm = Array{Int64}(undef, natoms(m))
  fragment_membership(m, fragm)
  for (i, f) in enumerate(fragm)
    expected[i] == f || return false
  end
  true
end

function test_atoms_in_fragment()::Bool
  smi = join(["C"^i for i in 1:10], '.')
  m = Molecule()
  build_from_smiles(m, smi) || return false
  for i in 0:9
    atoms_in_fragment(m, i) == (i + 1) || return false
  end
  true
end

function test_get_atoms_in_fragment()::Bool
  m = Molecule()
  smi = join(["C"^i for i in 1:10], '.')
  build_from_smiles(m, smi)
  for i in 0:9
    s = get_atoms_in_fragment(m, i)
    length(s) == (i + 1) || return false
  end
  true
end

function test_largest_fragment2()::Bool
  m = Molecule()
  smi = join(["C"^i for i in 1:10], '.')
  build_from_smiles(m, smi)
  largest_fragment(m) == 9 || return false
  true
end

function test_identify_spinach()::Bool
  m = LillyMol.MolFromSmiles("CC1C(O)C1N")
  spch = identify_spinach(m)
  spch[1] == 1 || return false
  spch[2] == 0 || return false
  spch[3] == 0 || return false
  spch[4] == 1 || return false
  spch[5] == 0 || return false
  spch[6] == 1 || return false
  true
end

function test_rings_in_fragment()::Bool
  m = LillyMol.MolFromSmiles("CCC.C1CC1.C1CC1C2CC2")
  rings_in_fragment(m, 0) == 0 || return false
  rings_in_fragment(m, 1) == 1 || return false
  rings_in_fragment(m, 2) == 2 || return false
  true
end

# Cannot figure out how to return a std::vector<Molecule>.
# This will change once I figure that out.
function test_create_components_not_invoked()::Bool
  m = LillyMol.MolFromSmiles("CCC.C1CC1.C1CC1C2CC2")
  nf = number_fragments(m)
  nf == 3 || return false
  frags = [Molecule() for i in 1:nf]
  create_components(m, frags)
  length(frags) == 3 || return false
  true
end

function test_create_components()::Bool
  m = LillyMol.MolFromSmiles("CCC.C1CC1.C1CC1C2CC2")
  nf = number_fragments(m)
  nf == 3 || return false
  components = Components()
  create_components(m, components)
  length(components) == 3 || return is_failure("Not 3 components", m)
  for (ndx, c) in enumerate(components)
    number_fragments(c) == 1 || return is_failure("Component not 1 frag", m)

    if ndx == 1
      smiles(c) == "CCC" || return is_failure("1 not CCC", m)
      nrings(c) == 0 || return is_failure("1 has rings", m)
    elseif ndx == 2
      smiles(c) == "C1CC1" || return is_failure("2 not C1CC1", m)
      nrings(c) == 1 || return is_failure("2 not 1 ring", m)
    elseif ndx == 3
      smiles(c) == "C1CC1C1CC1" || return is_failure("3 not C1CC1C2CC2", m)
      nrings(c) == 2 || return is_failure("3 not 2 rings", m)
    else
      return is_failure("Invalid index")
    end
  end
  for ring in rings(components[3])
    length(ring) == 3 || return is_failure("Not 3 membered", m)
  end
  true
end

function test_create_subset()::Bool
  m = LillyMol.MolFromSmiles("Cc1ccccc1N")
  subset = [0, 1, 1, 1, 1, 1, 1, 0]
  length(subset) == natoms(m) || return false
  s = create_subset(m, subset)
  unique_smiles(s) == "c1ccccc1" || return false
  true
end

function test_create_subset_set_of_atoms()::Bool
  m = LillyMol.MolFromSmiles("Cc1ccccc1N")
  # Initializer list does not work.
  # subset = SetOfAtoms(1,2,3,4,5,6)
  subset = SetOfAtoms()
  for i in 1:6
    LillyMol.push!(subset, i)
  end
  s = create_subset(m, subset)
  unique_smiles(s) == "c1ccccc1" || return false
end

function test_random_smiles()::Bool
  m = Molecule()
  build_from_smiles(m, "CC(=O)OC1=CC=CC=C1C(=O)O aspirin") || return is_failure("bad smiles", m)
  initial_smiles = unique_smiles(m)
  seen = Set()
  for i in 1:10
    smiles = random_smiles(m)
    smiles in seen && continue

    push!(seen, smiles)
    m2 = Molecule()
    build_from_smiles(m2, smiles) || return is_failure("invalid random smiles $(smiles)", m)
    unique_smiles(m2) == initial_smiles || return is_failure("Smiles mismatch", m)
  end
  true
end
  
function test_smiles_starting_atom()::Bool
  m = Molecule()
  build_from_smiles(m, "OCN") || return is_failure("Bad smiles", m)
  smiles(m) == "OCN" || return is_failure("Wrong smiles", m)
  smiles_starting_with_atom(m, 1) == "C(O)N" || return is_failure("Incorrect", m)
  true
end

function test_isotopically_labelled_smiles2()::Bool
  m = Molecule()
  build_from_smiles(m, "OCN") || return is_failure("Bad smiles", m)
  isotopically_labelled_smiles(m) == "O[1CH2][2NH2]" || return is_failure("wrong smiles", m)
  true
end

function test_symmetry()::Bool
  m = Molecule()
  build_from_smiles(m, "NC1=CC=C(C=C1)C(F)(F)F CHEMBL1162294") || return is_failure("Bad smiles", m)
  number_symmetry_classes(m) == 7 || return is_failure("Not 7 classes", m)
  symmetry_class(m, 2) == symmetry_class(m, 6) || return is_failure("2 6 not equivalent", m)
  symmetry_class(m, 3) == symmetry_class(m, 5) || return is_failure(" 3 5 not equivalent", m)

  symmetry_class(m, 8) == symmetry_class(m, 9) || return is_failure("8 9 not equivalent", m)
  symmetry_class(m, 9) == symmetry_class(m, 10) || return is_failure("9 10 not equivalent", m)
  symmetry_equivalents(m, 8) == [9, 10] || return is_failure("8 is not [9, 10]", m)
  isempty(symmetry_equivalents(m, 0)) || return is_failure("0 not empty", m)
  isempty(symmetry_equivalents(m, 4)) || return is_failure("4 not empty", m)
  isempty(symmetry_equivalents(m, 7)) || return is_failure("7 not empty", m)
  true
end


function test_reduce_to_largest_fragment()::Bool
  m = LillyMol.MolFromSmiles("CC")
  reduce_to_largest_fragment!(m)
  smiles(m) == "CC" || return false
  build_from_smiles(m, "C.C")
  reduce_to_largest_fragment!(m)
  smiles(m) == "C" || return false
  build_from_smiles(m, "C.N")
  reduce_to_largest_fragment!(m)
  smiles(m) == "C" || return false
  build_from_smiles(m, "N.C")
  reduce_to_largest_fragment!(m)
  smiles(m) == "N" || return false
  build_from_smiles(m, "C.CCCC")
  reduce_to_largest_fragment!(m)
  smiles(m) == "CCCC" || return false
  build_from_smiles(m, "CCCC.C")
  reduce_to_largest_fragment!(m)
  smiles(m) == "CCCC" || return false
  smi = join(Random.shuffle(["C"^i for i in 1:10]), '.')
  build_from_smiles(m, smi)
  reduce_to_largest_fragment!(m)
  smiles(m) == "C"^10 || return false
  true
end

function test_reduce_to_largest_organic_fragment()::Bool
  m = LillyMol.MolFromSmiles("[Na].C")
  reduce_to_largest_organic_fragment!(m)
  smiles(m) == "C" || return false
  true
end

function test_reduce_to_largest_fragment_carefully()::Bool
  m = LillyMol.MolFromSmiles("O=N(=O)c1cc(N(=O)=O)cc(N(=O)=O)c1O."* "OCN"^5)
  reduce_to_largest_fragment_carefully!(m)
  smiles(m) == "OCN"^5 || return false
  true
end
  
function test_fragment_related()::Bool
  m = Molecule()
  build_from_smiles(m, "C1CC1C1CC1") || return is_failure("cannot build", m)
  natoms(m) == 6 || return is_failure("Not 6 atoms", m)
  nrings(m) == 2 || return is_failure("Not 2 rings", m)
  number_fragments(m) == 1 || return is_failure("Not 1 frgament", m)
  in_same_ring_system(m, 0, 1) || return is_failure("0 and 1 not in same ring", m)
  in_same_ring_system(m, 1, 2) || return is_failure("1 and 2 not in same ring", m)

  in_same_ring_system(m, 2, 3) && return is_failure("2 and 3 in same ring", m)
  in_same_ring_system(m, 3, 4) || return is_failure("3 and 4 not in same ring", m)
  in_same_ring_system(m, 4, 5) || return is_failure("4 and 5 not in same ring", m)

  remove_bond_between_atoms!(m, 2, 3)
  natoms(m) == 6 || return is_failure("rmbond not 6 atoms", m)
  nrings(m) == 2 || return is_failure("rmbond not 2 rings", m)
  number_fragments(m) == 2 || return is_failure("rmbond not 2 frags", m)
  in_same_ring(m, 0, 1) || return is_failure("rmbond 0 and 1 not in same ring", m)
  in_same_ring(m, 1, 2) || return is_failure("rmbond 1 and 2 not in same ring", m)
  in_same_ring(m, 2, 3) && return is_failure("rmbond 2 and 3 in same ring", m)
  in_same_ring(m, 3, 4) || return is_failure("rmbond 3 and 4 not in same ring", m)
  in_same_ring(m, 4, 5) || return is_failure("rmbond 4 and 5 not in same ring", m)

  fused_system_identifier(m, 2) != fused_system_identifier(m, 3) || return is_failure("fsid error", m)
  atoms_in_fragment(m, 0) == 3 || return is_failure("not 3 atoms in frag 0", m);
  atoms_in_fragment(m, 1) == 3 || return is_failure("not 3 atoms in frag 1", m);

  fragment_membership(m, 0) == fragment_membership(m, 1) || return is_failure("frag membership 0 1", m)
  fragment_membership(m, 1) == fragment_membership(m, 2) || return is_failure("frag membership 1 2", m)
  fragment_membership(m, 3) == fragment_membership(m, 4) || return is_failure("frag membership 3 4", m)
  fragment_membership(m, 4) == fragment_membership(m, 5) || return is_failure("frag membership 4 5", m)
  # The actual numbers should not be depended on, this test is likely unstable
  label_atoms_by_ring_system(m) == [1, 1, 1, 2, 2, 2] || return is_failure("label_atoms_by_ring_system", m)

  # Restore to where we started
  add_bond!(m, 2, 3, SINGLE_BOND)
  number_fragments(m) == 1 || return is_failure("Rejoin not 1 frag", m)
  smiles(m) == "C1CC1C1CC1" || return is_failure("rejoin smiles", m)
  number_ring_systems(m) == 2 || return is_failure("rejoin number_ring_systems", m)
  label_atoms_by_ring_system(m) == [1, 1, 1, 2, 2, 2] || return is_failure("label_atoms_by_ring_system rejoin", m)
  true
end

function test_atoms_in_fragment2()::Bool
  m = Molecule()
  build_from_smiles(m, "CCCC.c1ccccc1") || return is_failure("Invalid smiles", m)
  atoms_in_fragment(m, 0) == 4 || return is_failure("Not 4 atoms in frag 0", m)
  atoms_in_fragment(m, 1) == 6 || return is_failure("Not 6 atoms in frag 1", m)
end


function test_organic_only()::Bool
  m = LillyMol.MolFromSmiles("IC(Cl)(Br)C(F)NCOCSP")
  organic_only(m) || return false
  build_from_smiles(m, "[B]")
  organic_only(m) && return false
  build_from_smiles(m, "[Si]")
  organic_only(m) && return false
  build_from_smiles(m, "[Se]")
  organic_only(m) && return false
  true
end

function test_contains_non_periodic_table_elements()::Bool
  LillyMol.set_auto_create_new_elements(1)
  m = Molecule()
  build_from_smiles(m, "[He][Ll]O[W][Rl][D]") || return false
  contains_non_periodic_table_elements(m) || return false
  LillyMol.set_auto_create_new_elements(0)
  true
end

function test_longest_path()::Bool
  smi = "C"^64
  m = Molecule()
  build_from_smiles(m, smi) || return false
  longest_path(m) == (length(smi) - 1) || return false
  true
end

function test_atoms_between()::Bool
  m = Molecule()
  build_from_smiles(m, "CC1CC(C)C1") || return false
  btw = atoms_between(m, 0, 4)
  length(btw) == 3 || return false
  btw[1] == 1 || return false
  btw[2] == 2 || return false
  btw[3] == 3 || return false
  true
end

function test_implicit_hydrogens()::Bool
  m = LillyMol.MolFromSmiles("CCC(C)C(C)(C)C")
  implicit_hydrogens(m, 0) == 3 || return false
  implicit_hydrogens(m, 1) == 2 || return false
  implicit_hydrogens(m, 2) == 1 || return false
  implicit_hydrogens(m, 4) == 0 || return false
  build_from_smiles(m, "CC(=O)O") || return false
  implicit_hydrogens(m, 0) == 3 || return false
  implicit_hydrogens(m, 1) == 0 || return false
  implicit_hydrogens(m, 2) == 0 || return false
  implicit_hydrogens(m, 3) == 1 || return false
  build_from_smiles(m, "NCNCN(C)C#N") || return false
  implicit_hydrogens(m, 0) == 2 || return false
  implicit_hydrogens(m, 2) == 1 || return false
  implicit_hydrogens(m, 4) == 0 || return false
  implicit_hydrogens(m, 7) == 0 || return false
  build_from_smiles(m, "C[N+](C)C") || return false
  valence_ok(m) && return false
  implicit_hydrogens(m, 1) == 0 || return false
  build_from_smiles(m, "C[NH+](C)C") || return false
  valence_ok(m) || return false
  implicit_hydrogens(m, 1) == 1 || return false
  build_from_smiles(m, "SC") || return false
  implicit_hydrogens(m, 0) == 1 || return false
  true
end

function test_explicit_hydrogens()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  explicit_hydrogens(m, 0) == 0 || return false
  build_from_smiles(m, "[H]C") || return false
  explicit_hydrogens(m, 1) == 1 || return false
  build_from_smiles(m, "[H]C[H]") || return false
  explicit_hydrogens(m, 1) == 2 || return false
  build_from_smiles(m, "[H]C([H])[H]") || return false
  explicit_hydrogens(m, 1) == 3 || return false
  build_from_smiles(m, "[H]C([H])([H])[H]") || return false
  explicit_hydrogens(m, 1) == 4 || return false
  true
end

function test_hcount()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  hcount(m, 0) == 4 || return false
  build_from_smiles(m, "[H]C")
  hcount(m, 1) == 4 || return false
  build_from_smiles(m, "[H]C[H]")
  hcount(m, 1) == 4 || return false
  build_from_smiles(m, "[H]C([H])[H]")
  hcount(m, 1) == 4 || return false
  build_from_smiles(m, "[H]C([H])([H])[H]")
  hcount(m, 1) == 4 || return false
  true
end

function test_make_implicit_hydrogens_explicit()::Bool
  m = Molecule()
  build_from_smiles(m, "CC")
  make_implicit_hydrogens_explicit!(m)
  smiles(m) == "C([H])([H])([H])C([H])([H])[H]" || return false
  true
end

function test_move_hydrogens_to_end_of_connection_table()::Bool
  m = Molecule()
  build_from_smiles(m, "CN")
  atomic_number(m, 0) == 6 || return false
  atomic_number(m, 1) == 7 || return false
  move_hydrogens_to_end_of_connection_table!(m, 6)
  atomic_number(m, 0) == 7 || return false
  atomic_number(m, 1) == 6 || return false
  smiles(m) == "NC" || return false
  true
end

function test_pi_electrons()::Bool
  m = Molecule()
  build_from_smiles(m, "OCNC=NC#N")
  pi_electrons(m, 0) == 2 || return false
  pi_electrons(m, 1) == 0 || return false
  pi_electrons(m, 2) == 2 || return false
  pi_electrons(m, 3) == 1 || return false
  pi_electrons(m, 4) == 1 || return false
  pi_electrons(m, 5) == 2 || return false
  pi_electrons(m, 6) == 2 || return false
  true
end

function test_lone_pair_count()::Bool
  m = Molecule()
  build_from_smiles(m, "OCNC=NC#N")
  lone_pair_count(m, 0) == 2 || return false
  lone_pair_count(m, 1) == 0 || return false
  lone_pair_count(m, 2) == 1 || return false
  lone_pair_count(m, 3) == 0 || return false
  lone_pair_count(m, 4) == 1 || return false
  lone_pair_count(m, 5) == 0 || return false
  lone_pair_count(m, 6) == 1 || return false
  true
end

function test_aromatic_atom_count()::Bool
  m = LillyMol.MolFromSmiles("Cc1nccc1")
  aromatic_atom_count(m) == 5 || return false
  true
end

function test_aromatic_ring_count()::Bool
  m = LillyMol.MolFromSmiles("c12ccccc2ccc(c1)c1ccccc1")
  aromatic_ring_count(m) == 3 || return false
  true
end

function test_saturated()::Bool
  m = LillyMol.MolFromSmiles("CC=Cc1ccccc1")
  saturated(m, 0) || return false
  saturated(m, 1) && return false
  saturated(m, 2) && return false
  saturated(m, 3) && return false
  true
end

function test_number_chiral_centres()::Bool
  m = Molecule()
  build_from_smiles(m, "C[C@H](F)N") || return false
  number_chiral_centres(m) == 1 || return false
  true
end

function test_chiral_centre_at_atom()::Bool
  m = LillyMol.MolFromSmiles("C[C@H](F)N")
  c = chiral_centre_at_atom(m, 1)
  top_front(c) == 0 || return false
  left_down(c) == 2 || return false
  right_down(c) == 3 || return false
  # TODO:ianwatson figure out lone pairs and implicit hydrogens
  true
end

function test_chiral_centre_in_molecule_not_indexed_by_atom_number()::Bool
  m = LillyMol.MolFromSmiles("C[C@H](F)N[C@@H](C)CC")
  c = chiral_centre_in_molecule_not_indexed_by_atom_number(m, 0)
  centre(c) == 1 || return false
  c = chiral_centre_in_molecule_not_indexed_by_atom_number(m, 1)
  centre(c) == 4 || return false
  true
end

function test_remove_chiral_centre_at_atom()::Bool
  m = LillyMol.MolFromSmiles("C[C@H](F)N[C@@H](C)CC")
  remove_chiral_centre_at_atom!(m, 4)
  smiles(m) == "C[C@H](F)NC(C)CC" || return false
  true
end

function test_remove_all_chiral_centres()::Bool
  m = LillyMol.MolFromSmiles("C[C@H](F)N[C@@H](C)CC")
  remove_all_chiral_centres!(m)
  smiles(m) == "CC(F)NC(C)CC" || return false
  true
end

function test_invert_chirality_on_atom()::Bool
  m = LillyMol.MolFromSmiles("C[C@H](F)N[C@@H](C)CC")
  invert_chirality_on_atom!(m, 4)
  smiles(m) == "C[C@H](F)N[C@H](C)CC" || return false
  true
end

function test_smarts_equivalent_for_atom()::Bool
  m = LillyMol.MolFromSmiles("Cc1ncccc1O")
  smarts_equivalent_for_atom(m, 0) == "[CD1H3v1]" || return false
  smarts_equivalent_for_atom(m, 1) == "[cD3H0v4;r6]" || return false
  true
end

function test_smarts()::Bool
  m = LillyMol.MolFromSmiles("Cn1cnc2n(C)c(=O)n(C)c(=O)c12")
  smarts(m) == "[C][n]1[c][n][c]2[c]1[c](=[O])[n]([C])[c](=[O])[n]2[C]" || return false
  true
end

function test_atom_map_number()::Bool
  m = Molecule()
  build_from_smiles(m, "[CH3:1]C[CH3:3]")
  atom_map_number(m, 0) == 1 || return false
  atom_map_number(m, 1) == 0 || return false
  atom_map_number(m, 2) == 3 || return false
  true
end
  
function test_atom_map_numbers()::Bool
  m = Molecule()
  build_from_smiles(m, "C[CH2:1][CH2:2]C") || return is_failure("Bad smiles", m)
  valence_ok(m) || return is_failure("valence", m)
  atom_map_number(m, 0) == 0 || return is_failure("Atom map 0", m)
  atom_map_number(m, 1) == 1 || return is_failure("Atom map 1", m)
  atom_map_number(m, 2) == 2 || return is_failure("Atom map 2", m)
  atom_map_number(m, 3) == 0 || return is_failure("Atom map 3", m)

  set_atom_map_number!(m, 3, 3)
  smiles(m) == "C[CH2:1][CH2:2][CH3:3]" || return is_failure("Smiles no good", m)

  atom_with_atom_map_number(m, 3) == 3 || return is_failure("not amap 3", m)

  reset_atom_map_numbers!(m)
  smiles(m) ==  "CCCC" || return is_failure("Atom map not reverted", m)
  true
end

function test_sort_atoms()::Bool
  m = Molecule()
  build_from_smiles(m, "CO.N") || return is_failure("Bad smiles", m)
  atomic_number(m, 0) == 6 || return is_failure("not carbon", m)
  atomic_number(m, 1) == 8 || return is_failure("not oxygen", m)
  atomic_number(m, 2) == 7 || return is_failure("not nitrogen", m)
  order = [1, 2, 0]
  sort_atoms!(m, order, 1)
  atomic_number(m, 0) == 8 || return is_failure("Not oxygen", m)
  atomic_number(m, 1) == 6 || return is_failure("Not carbon", m)
  atomic_number(m, 2) == 7 || return is_failure("Not nitrogen", m)
  smiles(m) == "OC.N" || return is_failure("Not sorted", m)
  true
end


function test_set_atom_map_number()::Bool
  m = LillyMol.MolFromSmiles("CCC")
  set_atom_map_number!(m, 1, 8)
  smiles(m) == "C[CH2:8]C" || return false
  atom_map_number(m, 1) == 8 || return false
  set_atom_map_number!(m, 1, 0)
  smiles(m) == "CCC" || return false
  true
end

function test_set_include_atom_map_with_smiles()::Bool
  m = Molecule();
  build_from_smiles(m, "[CH3:7]C") || return false
  set_include_atom_map_with_smiles(false)
  atom_map_number(m, 0) == 7 || return false
  # This is kind of unfortunate, perhaps this could be fixed/changed.
  smiles(m) == "[CH3]C" || return false
  set_include_atom_map_with_smiles(true)
  true
end

function test_atom_with_atom_map_number()::Bool
  m = Molecule()
  build_from_smiles(m, "C[CH2:4]C")
  atom_with_atom_map_number(m, 4) == 1 || return false
  atom_with_atom_map_number(m, 7) < 0 || return false
  true
end

function test_reset_all_atom_map_numbers()::Bool
  m = LillyMol.MolFromSmiles("C[CH2:5]C")
  reset_all_atom_map_numbers!(m)
  smiles(m) == "CCC" || return false
  true
end

function test_unset_unnecessary_implicit_hydrogens_known_values()::Bool
  m = Molecule()
  build_from_smiles(m, "[OH][CH2][NH2]")
  unset_unnecessary_implicit_hydrogens_known_values!(m)
  smiles(m) == "OCN" || return false
  true
end

function test_discern_chirality_from_3d_structure()::Bool
  m = Molecule()
  build_from_smiles(m, "C{{0,1,1}}C{{0,0,0}}(N{{0,1,-1}})(O{{-1,-1,0}})C{{1,-1,0}}C") || return false
  smiles(m) == "CC(N)(O)CC" || return false
  discern_chirality_from_3d_structure(m)
  smiles(m) == "C[C@](N)(O)CC" || return false
  setz!(m, 0, -1.0)
  setz!(m, 2, 1.0)
  discern_chirality_from_3d_structure(m)
  smiles(m) == "C[C@@](N)(O)CC" || return false
  true
end

function test_bonds_between()::Bool
  m = LillyMol.MolFromSmiles("Oc1ccc(C)cc1")
  bonds_between(m, 0, 1) == 1 || return false
  bonds_between(m, 0, 2) == 2 || return false
  bonds_between(m, 0, 3) == 3 || return false
  bonds_between(m, 0, 4) == 4 || return false
  bonds_between(m, 0, 5) == 5 || return false
  bonds_between(m, 0, 6) == 3 || return false
  bonds_between(m, 0, 7) == 2 || return false
  bonds_between(m, 1, 2) == 1 || return false
  bonds_between(m, 1, 3) == 2 || return false
  bonds_between(m, 1, 4) == 3 || return false
  bonds_between(m, 1, 5) == 4 || return false
  bonds_between(m, 1, 6) == 2 || return false
  bonds_between(m, 1, 7) == 1 || return false
  true
end

function test_returns_vector()::Bool
  v = LillyMol.returns_vector(1, 2)
  v == [1, 2] || return false
  true
end

function test_standardise()
  standardise = ChemicalStandardisation()
  activate_all(standardise)
  m = Molecule()
  build_from_smiles(m, "CC[N+](=O)[O-]")
  process(standardise, m)
  return smiles(m) == "CCN(=O)=O"
end

function test_aspirin()::Bool
  m = Molecule()
  build_from_smiles(m, "CC(=O)OC1=CC=CC=C1C(=O)O aspirin") || return is_failure("Cannot build aspirin", m)
  molecular_formula(m) == "C9H8O4" || return is_failure("Bad mf", m)
  isapprox(amw(m), 180.16, atol=0.01) || return is_failure("Bad amw", m)
  isapprox(exact_mass(m), 180.04225873, atol=0.000001) || return is_failure("Bad exact_mass", m)
  natoms(m) == 13 || return is_failure("not 13 atoms", m)
  nedges(m) == 13 || return is_failure("not 13 edges", m)
  nrings(m) == 1 || return is_failure("not 1 ring", m)
  aromatic_ring_count(m) == 1 || return is_failure("not 1 aromatic ring", m)
  aromatic_atom_count(m) == 6 || return is_failure("not 6 aromatic atoms", m)
  connections(m, 0) == [1] || return is_failure("Wrong connections", m)
  number_fragments(m) == 1 || return is_failure("not 1 fragment", m)
  connections(m, 1) == [0, 2, 3] || return is_failure("Wrong connections 1", m)
  organic_only(m) || return is_failure("Not organic", m);
  true
end

function test_xlogp()::Bool
  m = Molecule()
  build_from_smiles(m, "CC(=O)OC1=CC=CC=C1C(=O)O aspirin") || return is_failure("Bad smiles")
  (result, flag) = xlogp(m)
  # This does not work:
  # TypeError: non-boolean (CxxWrap.CxxWrapCore.CxxBool) used in boolean context
  #flag || return is_failure("xlogp failed", m)
  isapprox(result, 1.426, atol=0.001) || return is_failure("xlogp wrong", m)
  true
end

function test_alogp()::Bool
  smiles = [
    "S=C(NN=CC1=CC=C2OCOC2=C1)NCC(=C)C CHEMBL1989706",
    "S1C(=C(O)N2N=C(C)N=C12)C(N1CCOCC1)C1=CC=CO1 CHEMBL1578764",
    "O=C(OC)C1=CC(=CC(=C1)CN)NC(=O)C1=CC(=CC=C1)OC CHEMBL4904536",
    "C1=C(C)N=C(N=C1NC1=NNC(=C1)C1CC1)NCC1=NC=CN1 CHEMBL3964210",
    "S(=O)(=O)(N(C)CC(=O)NC1=CC(=CC=C1)OC)C1=CC=C2NC(=O)CCC2=C1 CHEMBL1513029",
    "CCS(=O)(=O)NC1=CC=C(OC2=CN=CN=C2)C(=C1)C1=CN(C)C(=O)C2=C1C=CN2 CHEMBL3967185",
    "C1(=CC=C2C(=C1)N=C1C(=CC(=O)C3=C1N=CC=C3)O2)C(=O)NCCCCC(N)C(=O)O CHEMBL4780009",
    "FC1=CC=C(C(=O)C2=C(NC(=O)CC3=COC4=C3C=CC(=C4C)C)C3=C(O2)C=CC=C3)C=C1 CHEMBL1436782",
    "S1C(=C(C2=C1CCCC2)C(=O)OCC)NC(=O)CSC1=CN(C2=C1C=CC=C2)CCNC(=O)C1=C(F)C=CC=C1F CHEMBL1438572",
    "C(#C)C(O)C=CCCCCCCC#CC(O)C#CCCCCC=CCCCCC=CCCCCCCCCCCCCCC=CC#C CHEMBL505112"
  ]

end
  
function test_cubane()::Bool
  m = Molecule()
  build_from_smiles(m, "C12C3C4C1C5C2C3C45") || return is_failure("Bad smiles", m)
  nrings(m) == 5 || return is_failure("Nrings not 5", m)
  non_sssr_rings(m) == 1 || return is_failure("Must be 1 non sssr rings", m)
  for atnum in eachindex(m)
    ring_bond_count(m, atnum) == 3 || return is_failure("rbc not 3", m)
  end

  for ring in sssr_rings(m)
    length(ring) == 4 || return is_failure("Ring not 4", m)
  end

  ring_membership(m) == [3, 3, 3, 3, 2, 2, 2, 2] || return is_failure("ring membership", m)

  fused_system_identifier(m, 0) == fused_system_identifier(m, 1) || return is_failure("fused_system_identifier", m)
  true
end
  
function test_distance_matrix()::Bool
  m = Molecule()
  build_from_smiles(m, "C1CCC1") || return is_failure("Bad smiles", m)
  bonds_between(m, 0, 1) == 1 || return is_failure("Bonds btw 1", m)
  bonds_between(m, 0, 2) == 2 || return is_failure("Bonds btw 2", m)

  build_from_smiles(m, "N(CC1=CC=C(OCCCC2=CC=CC=C2)C=C1)(CC1=CC=C(OCCCC2=CC=CC=C2)C=C1)CCCCN CHEMBL349114") || return is_failure("Bad smiles", m)
  fragment_membership(m, 30) == fragment_membership(m, 39) || return is_failure("30 30 frag", m)
  bonds_between(m, 30, 39) == 18 || return is_failure("30 39 dist", m)
  longest_path(m) == 26 || return is_failure("longest path", m)
  most_distant_pair(m) == (13, 30) || return is_failure("most distanct pair", m)
  bonds_between(m, 13, 30) == 26 || return is_failure("13, 30", m)
  true
end

#      C
#    / | \
# C-C  |  C-C
#    \ | /
#      C
#
function test_ring_related()::Bool
  m = Molecule()
  build_from_smiles(m, "CC1C2C(C)C12") || return is_failure("Bad smiles", m)
  natoms(m) == 6 || return is_failure("not 6 atoms", m)
  natoms(m, 6) == 6 || return is_failure("not 6 carbon atoms", m)

  ncon(m, 0) == 1 || return is_failure("ncon 0", m)
  ncon(m, 1) == 3 || return is_failure("ncon 1", m)
  ncon(m, 2) == 3 || return is_failure("ncon 2", m)
  ncon(m, 3) == 3 || return is_failure("ncon 3", m)
  ncon(m, 4) == 1 || return is_failure("ncon 4", m)
  ncon(m, 5) == 3 || return is_failure("ncon 5", m)

  nrings(m) == 2 || return is_failure("not 2 rings", m)
  non_sssr_rings(m) == 0 || return is_failure("non sssr rings", m)
  nrings(m, 0) == 0 || return is_failure("nrings 0", m)
  nrings(m, 1) == 1 || return is_failure("nrings 1", m)
  nrings(m, 2) == 2 || return is_failure("nrings 2", m)
  nrings(m, 3) == 1 || return is_failure("nrings 3", m)
  nrings(m, 4) == 0 || return is_failure("nrings 4", m)
  nrings(m, 5) == 2 || return is_failure("nrings 5", m)

  ring_bond_count(m, 0) == 0 || return is_failure("ring_bond_count 0", m)
  ring_bond_count(m, 1) == 2 || return is_failure("ring_bond_count 1", m)
  ring_bond_count(m, 2) == 3 || return is_failure("ring_bond_count 2", m)
  ring_bond_count(m, 3) == 2 || return is_failure("ring_bond_count 3", m)
  ring_bond_count(m, 4) == 0 || return is_failure("ring_bond_count 4", m)
  ring_bond_count(m, 5) == 3 || return is_failure("ring_bond_count 5", m)

  for i in eachindex(m)
    attached_heteroatom_count(m, i) == 0 || return is_failure("Attached heteroatom count", m)
  end

  is_ring_atom(m, 0) && return is_failure("0 is ring atom", m)
  is_ring_atom(m, 1) || return is_failure("1 not ring atom", m)

  in_ring_of_given_size(m, 0, 3) && return is_failure("0 in ring", m)
  in_ring_of_given_size(m, 1, 3) || return is_failure("1 not in 3", m)
  in_ring_of_given_size(m, 1, 4) && return is_failure("1 not in 4", m)
  in_ring_of_given_size(m, 2, 3) || return is_failure("2 not in 3", m)
  in_ring_of_given_size(m, 3, 3) || return is_failure("3 not in 3", m)
  in_ring_of_given_size(m, 5, 3) || return is_failure("5 not in 3", m)

  in_same_ring(m, 1, 2) || return is_failure("1 2 not in same ring", m)
  in_same_ring(m, 1, 5) || return is_failure("1 5 not in same ring", m)
  in_same_ring(m, 3, 2) || return is_failure("3 2 not in same ring", m)
  in_same_ring(m, 3, 5) || return is_failure("3 5 not in same ring", m)
  in_same_ring(m, 1, 3) && return is_failure("1 3 in same ring", m)

  fused_system_identifier(m, 1) == fused_system_identifier(m, 2) || return is_failure("fsid mismatch 1, 2", m)
  fused_system_identifier(m, 1) == fused_system_identifier(m, 3) || return is_failure("fsid mismatch 1, 3", m)
  fused_system_identifier(m, 1) == fused_system_identifier(m, 5) || return is_failure("fsid mismatch 1, 5", m)
  fused_system_identifier(m, 1) == fused_system_identifier(m, 5) || return is_failure("fsid mismatch 1, 5", m)

  fused_system_size(m, 1) == 2 || return is_failure("fused_system_size 1", m)
  fused_system_size(m, 0) == 0 || return is_failure("fused_system_size 0", m)

  nrings(m) == 2 || return is_failure("nrings 2", m)
  # The order of the rings should not be relied upon, so this test is
  # possibly fragile.
  # cannot get this to work yet.
  # ring(m, 0) == [1, 2, 5] || return is_failure("ring0", m)
  # ring(m, 1) == [2, 3, 5] || return is_failure("ring1", m)
  for r in rings(m)
    length(r) == 3 || return is_failure("Ring not 3", m)
  end

  largest_ring_size(m) == 3 || return is_failure("largest_ring_size", m)
  number_ring_systems(m) == 1 || return is_failure("number_ring_systems". m)
  is_spiro_fused(m, 1) && return is_failure("Spiro 1", m)
  is_spiro_fused(m, 2) && return is_failure("Spiro 2", m)
  return true
end
  

function test_atom_iterator()::Bool
  m = Molecule()
  build_from_smiles(m, "C1CC1") || return is_failure("Bad smiles", m)
  valence_ok(m) || return is_failure("Valence not ok", m)
  # Iterate through all atoms in the molecule.
  for (ndx, atom) in enumerate(m)
    atomic_number(m, ndx - 1) == 6 || return is_failure("Not carbon", m)
    ring_bond_count(m, ndx - 1) == 2 || return is_failure("Not 2 ring bonds", m)
  end

  build_from_smiles(m, "CC(N)(O)F") || return is_failure("Bad smiles", m)
  # Loop through bonds connected to atom 1. Check with `m` that
  # the atom is bonded to 1.
  for bond in m[1] 
    atom = other(bond, 1)
    are_bonded(m, atom, 1) || return is_failure("Not bonded to 1", m)
  end
  true
end


function test_iterate_bonds()::Bool
  m = Molecule()
  build_from_smiles(m, "O=C1N(C(=O)N2[C@@]1([C@@H]1[C@H]([C@H]2C2=CC=CC(=C2)C)C(=O)N(C1=O)CC1=CC=CC=C1)CC)C1=CC=CC=C1 CHEMBL1434237") || return is_failure("Bad smiles", m)
  # Count number of aromatic bonds
  compute_aromaticity_if_needed(m)
  aromatic_bonds = 0
  for bond in bonds(m)
    if is_aromatic(bond)
      aromatic_bonds += 1
    end
  end
  aromatic_bonds == 18 || return is_failure("Not 18 aromatic bond", m)
  return true
end

function test_index_bond_list()::Bool
  m = Molecule();
  build_from_smiles(m, "CC=CC#C") || return is_failure("Bad smiles")
  for (ndx, bond) in enumerate(bonds(m))
    if ndx == 1
       is_single_bond(bond) || return is_failure("1 not single", m)
    elseif ndx == 2
       is_double_bond(bond) || return is_failure("2 not double", m)
    elseif ndx == 3
       is_single_bond(bond) || return is_failure("2 not single", m)
    elseif ndx == 4
       is_triple_bond(bond) || return is_failure("3 not triple", m)
    else
      return is_failure("Index out of range", m)
    end
  end
  true
end

function test_loop_bond_list()::Bool
  m = Molecule()
  build_from_smiles(m, "CC=CC#C")
  blist = bonds(m)
  for i in 1:nedges(m)
    if i == 1
       is_single_bond(blist[i]) || return is_failure("0 not single", m)
    elseif i == 2
       is_double_bond(blist[i]) || return is_failure("1 not double", m)
    elseif i == 3
       is_single_bond(blist[i]) || return is_failure("2 not single", m)
    elseif i == 4
       is_triple_bond(blist[i]) || return is_failure("3 not triple", m)
    else
      return is_failure("Index out of range", m)
    end
  end
  true
end

function build_benzene()::Bool
  m = Molecule()
  for i in 1:6
    add!(m, 6)
  end

  natoms(m) == 6 || return is_failure("Not 6 atoms", m)
  natoms(m, 6) == 6 || return is_failure("Not 6 carbon", m)

  add_bond!(m, 0, 1, SINGLE_BOND)
  add_bond!(m, 1, 2, DOUBLE_BOND)
  add_bond!(m, 2, 3, SINGLE_BOND)
  add_bond!(m, 3, 4, DOUBLE_BOND)
  add_bond!(m, 4, 5, SINGLE_BOND)
  add_bond!(m, 5, 0, DOUBLE_BOND)

  aromatic_ring_count(m) == 1 || return is_failure("Not 1 aromatic ring", m)
  return true
end

function test_atom_iterator_and_valence()::Bool
  m = Molecule()
  build_from_smiles(m, "[1CH3]C(=[2CH2])#[3CH]") || return is_failure("Bad smiles", m)
  valence_ok(m) && return is_failure("ok valence", m)

  for bond in m[1]
    o = other(bond, 1)
    if is_single_bond(bond)
        isotope(m, o) == 1 || return is_failure("Not isotope 1", m)
    elseif is_double_bond(bond)
      isotope(m, o) == 2 || return is_failure("Not isotope 2", m)
    elseif is_triple_bond(bond)
      isotope(m, o) == 3 || return is_triple_bond("Not isotope 3", m)
    end
  end
  return true
end

function test_index_atom()::Bool
  m = Molecule()
  # build_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C caffeine") || return is_failure("Bad smiles")
  build_from_smiles(m, "CC(N)(O)C") || return is_failure("Bad smiles")
  length(m[0]) == 1  || return is_failure("Atom 0 not 1 connection", m)
  length(m[1]) == 4  || return is_failure("Atom 1 not 4 connection", m)
  for i in 2:4
    length(m[i]) == 1  || return is_failure("Atom $(i) not 1 connection", m)
  end
  # Loop through atoms atttached to atom 1
  for i in 1:4
    bond = m[1][i]
    is_single_bond(bond) || return is_failure("Not single bond", m)

    o = other(bond, 1)
    if i == 1
      o == 0 || return is_failure("0 not 0", m)
      atomic_number(m, o) == 6 || return is_failure("0 not carbon", m)
    elseif i == 2
      o == 2 || return is_failure("1 not 2", m)
      atomic_number(m, o) == 7 || return is_failure("1 not nitrogen", m)
    elseif i == 3
      o == 3 || return is_failure("2 not 3", m)
      atomic_number(m, o) == 8 || return is_failure("2 not oxygen", m)
    elseif i == 4
      o == 4 || return is_failure("3 not 4", m)
      atomic_number(m, o) == 6 || return is_failure("3 not carbon", m)
    end
  end

  true
end

function test_imidazole()::Bool
  m = Molecule()
  build_from_smiles(m, "N1C=NC2=CC=CC=C12") || return is_failure("Bad smiles", m)
  compute_aromaticity_if_needed(m)
  for ring in rings(m)
    is_aromatic(ring) || return is_failure("benzimidazole ring not aromatic", m)
  end
  true
end

# Do a 'substructure search'. Generally do not do this, do an
# actual substructure search, but this might be illustrative.
# Identify a O=a atom pair.
# There are many ways this could be done, and which will be most
# efficient will depend on the molecules being examined.
function test_find_exocyclic_bond()::Bool
  m = Molecule()
  build_from_smiles(m, "N1NC(=O)C=C1 CHEMBL4227850") || return is_failure("Bad smiles", m)
  aromatic_ring_count(m) == 1 || return is_failure("Not 1 aromatic ring", m)

  # Start search with singly bonded O atoms
  found_match = 0
  for (ndx, atom) in enumerate(m)
    ndx = ndx - 1
    ncon(atom) == 1 || continue
    atomic_number(atom) == 8 || continue

    # Bond to the first neighbour
    bond = atom[1]
    is_double_bond(bond) || continue
    o = other(bond, ndx)
    if is_aromatic(m, o) 
      found_match  += 1
    end
  end

  found_match > 0 || return is_failure("No aromatic found", m)

  # Start search by looking at atoms in aromatic rings
  found_match = 0
  compute_aromaticity_if_needed(m)
  for ring in sssr_rings(m)
    is_aromatic(ring) || continue
    for atom_number in ring
      atom = m[atom_number]
      ncon(atom) == 2 && continue

      for bond in atom
        is_double_bond(bond) || continue
        o = other(bond, atom_number)
        ncon(m, o) == 1 || continue
        if atomic_number(m, o) == 8
          found_match += 1
          break
        end
      end
    end
  end
  found_match > 0 || return is_failure("No aromat found II", m)

# Start search with 3 connected aromatic atoms
  found_match = 0
  for (ndx, atom) in enumerate(m)
    ndx = ndx - 1
    ncon(atom) == 3 || continue

    is_aromatic(m, ndx) || continue
    for bond in atom
      is_double_bond(bond) || continue
      o = other(bond, ndx)
      ncon(m, o) == 1 || continue
      atomic_number(m, o) == 8 || continue
      found_match += 1
    end
  end
  found_match > 0 || return is_failure("No aromatic found III", m)

  # On 50k random Chembl molecules the times for the methods are (seconds)
  # 1 5.18
  # 2 9.14
  # 3 8.47
  # In this case, the singly bonded oxygen is the best place to start
  # since they are comparatively rare.
  true
end

#function test_gather_rings()::Bool
#  m = Molecule()
#  build_from_smiles(m, "C1=CC=C2C(=C1)NC=N2")
#  rings = RingAtoms()
#  gather_rings(m, rings)
#  length(rings) == 2 || return is_failure("Not 2 rings", m)
#
#  rings[1] == [3, 8, 7, 6, 4] || return is_failure("Bad ring 1", m)
#  rings[2] == [0, 5, 4, 3, 2, 1] || return is_failure("Bad ring 2", m)
#
#  true
#end

function test_gather_rings()::Bool
  m = Molecule()
  build_from_smiles(m, "C1=CC=C2C(=C1)NC=N2")
  rings = RingInformation()
  gather_rings(m, rings, RIP_NONE)
  length(rings) == 2 || return is_failure("Not 2 rings", m)

  rings[1] == [3, 8, 7, 6, 4] || return is_failure("Bad ring 1", m)
  rings[2] == [0, 5, 4, 3, 2, 1] || return is_failure("Bad ring 2", m)

  true
end

function iterate_ring_information()::Bool
  m = Molecule()
  build_from_smiles(m, "C1=CC=C2C(=C1)NC=N2")
  rings = RingInformation()
  gather_rings(m, rings, RIP_NONE)
  length(rings) == 2 || return is_failure("Not 2 rings", m)

  for (ndx, ring) in enumerate(rings)
    if ndx == 1
      ring == [3, 8, 7, 6, 4] || return is_failure("Bad ring 1", m)
    elseif ndx == 2
      ring == [0, 5, 4, 3, 2, 1] || return is_failure("Bad ring 2", m)
    else
      return is_failure("Bad index")
    end
  end

  true
end

function test_scaffold()::Bool
  m = Molecule()
  build_from_smiles(m, "C") || return is_failure("Bad smiles")
  to_scaffold!(m)
  smiles(m) == "C" || return is_failure("Scaffold not C", m)

  build_from_smiles(m, "O1N=C(C(=O)N2CCCC2)C=C1COC1=CC=C2N=CC=CC2=C1 CHEMBL1589003") || return is_failure("Bad smiles", m)
  to_scaffold!(m)
  smiles(m) == "O1N=C(CN2CCCC2)C=C1COC1=CC=C2N=CC=CC2=C1" || return is_failure("Wrong scaffold", m)

  build_from_smiles(m, "O=C(N(C1=CC=C(C)C=C1)CC(=O)NCCOC)CCC(=O)NC1=CC=CC=N1 CHEMBL1576099") || return is_failure("Bad smiles", m)
  to_scaffold!(m)
  smiles(m) == "C(NC1=CC=CC=C1)CCCNC1=CC=CC=N1" || return is_failure("Scaffold incorrect", m)

  build_from_smiles(m, "O=C1N(C(=O)C2=C1C(=CC=C2)N(=O)=O)CC(=O)N1CC2=CC=CC=C2CC1 CHEMBL2134451") || return is_failure("Bad smiles", m)
  to_scaffold!(m)
  smiles(m) == "C1N(CC2=C1C=CC=C2)CCN1CC2=CC=CC=C2CC1" || return is_failure("Scaffold not good", m)

  build_from_smiles(m, "O=C(C1=CC=CN1CC(=O)NCC1N(CCC1)CC)C1=CC=CC=C1C CHEMBL1404612") || return is_failure("Bad Smiles", m)
  to_scaffold!(m)
  smiles(m) == "C(C1=CC=CN1CCNCC1NCCC1)C1=CC=CC=C1" || return is_failure("scaffold not ok", m)

  build_from_smiles(m, "O=C(N1[C@H](C(=O)NC2C3=CC=CC=C3CCC2)CCC1)[C@@H](NC(=O)[C@H](C)NC)CC(=O)O CHEMBL1570483") || return is_failure("Bad smiles", m)
  to_scaffold!(m)
  smiles(m) == "N1[C@H](CNC2C3=CC=CC=C3CCC2)CCC1" || return is_failure("scaffold not formed ok", m)

  return true
end

function test_coords()::Bool
  m = Molecule()
  build_from_smiles(m, "C{{0,0,0}}C{{1,1,1}}C{{2,0,0}}") || return is_failure("Bad smiles", m)
  highest_coordinate_dimensionality(m) == 3 || return is_failure("Not 3d", m)
  isapprox(x(m, 0), 0.0, atol=1.0e-05) || return is_failure("X not 0", m)
  isapprox(y(m, 0), 0.0, atol=1.0e-05) || return is_failure("Y not 0", m)
  isapprox(z(m, 0), 0.0, atol=1.0e-05) || return is_failure("Z not 0", m)

  isapprox(x(m, 1), 1.0, atol=1.0e-05) || return is_failure("X not 1", m)
  isapprox(y(m, 1), 1.0, atol=1.0e-05) || return is_failure("Y not 1", m)
  isapprox(z(m, 1), 1.0, atol=1.0e-05) || return is_failure("Z not 1", m)

  isapprox(distance_between_atoms(m, 0, 1), sqrt(3.0), atol=1.e0-5) || return is_failure("0 1 wroing distance", m)
  isapprox(distance_between_atoms(m, 0, 2), 2.0, atol=1.e0-5) || return is_failure("0 2 wroing distance", m)

  isapprox(bond_angle(m, 0, 1, 2), 1.230959, atol=1.0e-05) || return is_failure("Bad bond angle", m)
  return true
end
  
function test_getxyz()::Bool
  m = Molecule()
  build_from_smiles(m, "C{{0,0,0}}C{{1,1,1}}C{{2,0,0}}") || return is_failure("Bad smiles", m)
  coords = get_coordinates(m)
  size(coords) == (natoms(m), 3) || return is_failure("Wrong size", m)

  atol = 1.0e-05
  isapprox(coords[1, 1], 0.0, atol=atol) || return is_failure("1,1 not 0", m)
  isapprox(coords[1, 2], 0.0, atol=atol) || return is_failure("1,2 not 0", m)
  isapprox(coords[1, 3], 0.0, atol=atol) || return is_failure("1,3 not 0", m)

  isapprox(coords[2, 1], 1.0, atol=atol) || return is_failure("2,1 not 1", m)
  isapprox(coords[2, 2], 1.0, atol=atol) || return is_failure("2,2 not 1", m)
  isapprox(coords[2, 3], 1.0, atol=atol) || return is_failure("2,3 not 1", m)

  isapprox(coords[3, 1], 2.0, atol=atol) || return is_failure("3,1 not 2", m)
  isapprox(coords[3, 2], 0.0, atol=atol) || return is_failure("3,2 not 0", m)
  isapprox(coords[3, 3], 0.0, atol=atol) || return is_failure("3,3 not 0", m)

  return true
end

function test_setxyz()::Bool
  m = Molecule()
  build_from_smiles(m, "C"^5) || return is_failure("Bad smiles", m)
  matoms = natoms(m)
  coords = Matrix{Float32}(undef, matoms, 3)
  for i in 1:matoms
    for j in 1:3
      coords[i,j] = 4 * i + j
    end
  end

  set_xyz!(m, coords)

  atol = 1.0f-05
  for i in 1:matoms
    atnum = i - 1;
    for j in 1:3
      isapprox(convert(Float32, x(m, atnum)), convert(Float32, 4 * i + 1), atol=atol) || return is_failure("Bad x", m)
      isapprox(convert(Float32, y(m, atnum)), convert(Float32, 4 * i + 2), atol=atol) || return is_failure("Bad y", m)
      isapprox(convert(Float32, z(m, atnum)), convert(Float32, 4 * i + 3), atol=atol) || return is_failure("Bad z", m)
    end
  end

  true
end

# Cannot get this to work.
# MethodError: no method matching Array{Float32, 3}(::Vector{Float32})
# presumably something to do with the multi-dimensional array.
function test_dihedral_scan()::Bool
  m = Molecule()
  build_from_smiles(m, "C{{-2,1,0}}C{{-1,0,0}}C{{0,0,0}}C{{1,1,0}}") || return is_failure("Bad smiles", m)
  bump_check = 0.0 
  angle = 45.0
  coords = dihedral_scan(m, 1, 2, angle, bump_check)
# length(coords) == 7 || return is_failure("Not 7", m)

  expected = [-45.0, -90.0, -135.0, 180, 135, 90, 45]
  atol = 1.0f-04
# for i in 1:length(coords)
#    set_coordinates!(m, coords[:,:,i])
#    found = signed_dihedral_angle(m, 0, 1, 2, 3)
#    isapprox(found * 180.0 / 3.14159265, expected[i], atol=atol)
# end
  true
end

function test_rule_of_five()::Bool
  m = Molecule()
  build_from_smiles(m, "COC1=CC=C(C=C1)N2C3=C(CCN(C3=O)C4=CC=C(C=C4)N5CCCCC5=O)C(=N2)C(=O)N Eliquis") || return is_failure("Bad smiles", m)
  donor = 0
  acceptor = 0

  # no consideration of charged atoms.
  for (ndx, atom) in enumerate(m)
    z = atomic_number(atom)
    # Intercept the most common case first
    z == 6 && continue

    h = hcount(m, ndx - 1)
    if z == 7
      if h == 0
        acceptor += 1
      else
        donor += h
      end
    elseif z == 8
      if h == 0
        acceptor += 1
      else
        donor += 1
      end
    end
  end

  donor ==  2 || return is_failure("Not 2 donors", m)
  acceptor == 8 || return is_failure("Not 8 acceptors", m)
  true
end

CHIRAL_SMILES1 = [
"C(O)[C@@H](N)C CHEMBL1229871",
"O1[C@H](C1)CF CHEMBL501668",
"O1[C@H](C1)CBr CHEMBL504705",
"O1[C@H](C1)CCl CHEMBL448626",
"BrC[C@@H](C)O CHEMBL446288",
"C(Cl)[C@@H](C)Cl CHEMBL373466",
"SC[C@@H](N)C CHEMBL37279",
"OC[C@H](N)CC CHEMBL3184640",
"C1(=O)[C@H](N)CO1 CHEMBL2219717",
"ClC[C@H](O)CO CHEMBL1794186",
]
CHIRAL_SMILES2 = [
"O[C@H]1[C@H](N)CCC1 CHEMBL2374489",
"O[C@H]1[C@H](N)CCC1 CHEMBL2375114",
"N1C[C@H](O)[C@@H](O)C1 CHEMBL2335511",
"C1[C@@H](O)[C@@H](O)CN1 CHEMBL2207396",
"C1[C@@H](N)[C@H]1C(=O)O CHEMBL403157",
"C1NC[C@@H](O)[C@@H]1O CHEMBL396701",
"C1C[C@@H](O)[C@@H](O)C1 CHEMBL399324",
"C1=NC[C@@H](O)[C@H]1O CHEMBL389969",
"O1[C@H](C)[C@@H]1C(=O)O CHEMBL370643",
"O1C[C@@H](O)[C@H](O)C1 CHEMBL350524",
]

function test_number_chiral_smiles1()::Bool
  CHIRAL_SMILES1 = [
    "C(O)[C@@H](N)C CHEMBL1229871",
    "O1[C@H](C1)CF CHEMBL501668",
    "O1[C@H](C1)CBr CHEMBL504705",
    "O1[C@H](C1)CCl CHEMBL448626",
    "BrC[C@@H](C)O CHEMBL446288",
    "C(Cl)[C@@H](C)Cl CHEMBL373466",
    "SC[C@@H](N)C CHEMBL37279",
    "OC[C@H](N)CC CHEMBL3184640",
    "C1(=O)[C@H](N)CO1 CHEMBL2219717",
    "ClC[C@H](O)CO CHEMBL1794186",
  ]
  mols = [LillyMol.MolFromSmiles(s) for s in CHIRAL_SMILES1]
  for m in mols
    number_chiral_centres(m) == 1 || return is_failure("Not 1 chiral centre", m)
  end
  true
end

function test_number_chiral_smiles2()::Bool
  CHIRAL_SMILES2 = [
    "O[C@H]1[C@H](N)CCC1 CHEMBL2374489",
    "O[C@H]1[C@H](N)CCC1 CHEMBL2375114",
    "N1C[C@H](O)[C@@H](O)C1 CHEMBL2335511",
    "C1[C@@H](O)[C@@H](O)CN1 CHEMBL2207396",
    "C1[C@@H](N)[C@H]1C(=O)O CHEMBL403157",
    "C1NC[C@@H](O)[C@@H]1O CHEMBL396701",
    "C1C[C@@H](O)[C@@H](O)C1 CHEMBL399324",
    "C1=NC[C@@H](O)[C@H]1O CHEMBL389969",
    "O1[C@H](C)[C@@H]1C(=O)O CHEMBL370643",
    "O1C[C@@H](O)[C@H](O)C1 CHEMBL350524",
  ]
  mols = [LillyMol.MolFromSmiles(s) for s in CHIRAL_SMILES2]
  for m in mols
    number_chiral_centres(m) == 2 || return is_failure("Not 2 chiral centre", m)
  end
  true
end

function test_chiral_implicit_hydrogen()::Bool
  m = Molecule()
  build_from_smiles(m, "C[C@H](N)F") || return is_failure("Bad smiles", m)
  number_chiral_centres(m) == 1 || return is_failure("not 1 chiral centres", m)
  # self.assertIsNone(m.chiral_centre_at_atom(0))
  c = chiral_centre_at_atom(m, 1)
  # self.assertIsNotNone(c)
  # Which atom gets assigned to which position is unpredictable.
  top_front(c) == 0 || return is_failure("Top front", m)
  left_down(c) == 2 || return is_failure("Left down", m)
  right_down(c) == 3 || return is_failure("Right down", m)
  # TODO:ianwatson implement something sensible for this.
  is_chiral_implicit_hydrogen(top_back(c)) || return is_failure("Not implicit H", m)

  # All smiles variants must be identical
  usmi = unique_smiles(m)

  for i in 1:10
    m2 = LillyMol.MolFromSmiles(random_smiles(m))
    usmi == unique_smiles(m2) || return is_failure("Usmi mismatch", m)
  end
  return true
end

function test_invert_chirality()::Bool
  m = Molecule()
  build_from_smiles(m, "C[C@H](N)F") || return is_failure("bad smiles", m)
  usmi = unique_smiles(m)
  invert_chirality_on_atom!(m, 1)
  usmi == unique_smiles(m) && return is_failure("Usmi unchanged", m)

  remove_chiral_centre_at_atom!(m, 1)
  usmi == unique_smiles(m) && return is_failure("Usmi same after removal", m)
  occursin("@", unique_smiles(m)) && return is_failure("@ in usmi", m)
  return true
end
  
function test_iterate_chiral_centres()::Bool
  m = Molecule();
  build_from_smiles(m, "O[C@H]1[C@@H](O)C[C@@H](N)[C@H]1O CHEMBL268037") || return is_failure("Bad smiles", m)
  number_chiral_centres(m) == 4 || return is_failure("Not 4 chiral centres", m)
  atoms = SetOfAtoms()

  for c in chiral_centres(m)
    push!(atoms, centre(c))
  end
  atoms == [1, 2, 5, 7] || return is_failure("Not right atoms $(atoms)", m)
  true
end

function test_index_chiral_centres()::Bool
  m = Molecule();
  build_from_smiles(m, "O[C@H]1[C@@H](O)C[C@@H](N)[C@H]1O CHEMBL268037") || return is_failure("Bad smiles", m)
  number_chiral_centres(m) == 4 || return is_failure("Not 4 chiral centres", m)
  atoms = SetOfAtoms()

  for (ndx, c) in enumerate(chiral_centres(m))
    if ndx == 1
      centre(c) == 1 || return is_failure("1 not 1", m)
    elseif ndx == 2
      centre(c) == 2 || return is_failure("2 not 2", m)
    elseif ndx == 3
      centre(c) == 5 || return is_failure("3 not 5", m)
    elseif ndx == 4
      centre(c) == 7 || return is_failure("4 not 1", m)
    else
      return is_failure("Invalid index $(ndx)", m)
    end
  end
  true
end


function test_charge()::Bool
  m = Molecule()
  build_from_smiles(m, "CC(=O)[O-]") || return is_failure("Bad smiles", m)
  net_formal_charge(m) == -1 || return is_failure("net_formal_charge wrong", m)
  number_formal_charges(m) == 1 || return is_failure("Not 1 formal charge", m)
  has_formal_charges(m) || return is_failure("has_formal_charges", m)

  m2 = Molecule()
  build_from_smiles(m2, "CC(=O)O") || return is_failure("Bad smiles", m2)
  formal_charge(m2, 3) ==  0 || return is_failure("O has charge", m2)
  net_formal_charge(m2) ==  0 || return is_failure("neutral has charge", m2)
  set_formal_charge!(m2, 3, -1)
  formal_charge(m2, 3) == -1 || return is_failure("set_formal_charge failed", m2)
  net_formal_charge(m2) == -1 || return is_failure("Did not get -1 net", m2)

  return true
end

function write_smiles_tempfile(prefix::String, smiles::String)::String
  tmpdir = Base.Filesystem.mktempdir()
  fname = joinpath(tmpdir, "$(prefix).smi")
  open(fname, "w") do file
    write(file, smiles)
  end

  fname
end

function test_read_smiles()::Bool
  m = Molecule()

  smiles = "C methane\nCC ethane\nCCC propane\nC1CC1 cyclopropane\nCCCC butane\nc1ccccc1 benzene\n"

  fname = write_smiles_tempfile("abc", smiles)

  nsmiles = length(split(smiles, "\n")) - 1

  reader = LillyMol.MoleculeReader(SMI, fname)

  molecules_remaining(reader) == nsmiles || return is_failure("Wrong count remaining", m)

  my_count = 0
  while next_molecule(reader, m)
    my_count += 1
  end
  my_count == molecules_read(reader) || return is_failure("Incorrect count", m)
  molecules_read(reader) == nsmiles || return is_failure("Wrong molecule count", m)
  true
end

function test_read_smiles_with_errors_fail()::Bool
  m = Molecule()

  set_display_smiles_interpretation_error_messages(false)

  smiles = "C methane\nCC ethane\nCCCQ propane\nC1CC1 cyclopropane\nCCCC butane\nc1ccccc1 benzene\n"

  fname = write_smiles_tempfile("abc", smiles)

  nsmiles = length(split(smiles, "\n")) - 1

  reader = LillyMol.MoleculeReader(SMI, fname)

  molecules_remaining(reader) == nsmiles || return is_failure("Wrong count remaining", m)

  my_count = 0
  while next_molecule(reader, m)
    my_count += 1
  end
  my_count == molecules_read(reader) || return is_failure("Incorrect count", m)
  molecules_read(reader) == 2 || return is_failure("Wrong molecule count", m)

  set_display_smiles_interpretation_error_messages(true)

  true
end

function test_read_smiles_with_errors_skip()::Bool
  m = Molecule()

  set_display_smiles_interpretation_error_messages(false)

  smiles = "C methane\nCC ethane\nCCCQ propane\nC1CC1 cyclopropane\nCCCC butane\nc1ccccc1 benzene\n"

  fname = write_smiles_tempfile("abc", smiles)

  nsmiles = length(split(smiles, "\n")) - 1

  reader = LillyMol.MoleculeReader(SMI, fname)

  molecules_remaining(reader) == nsmiles || return is_failure("Wrong count remaining", m)
  set_connection_table_errors_allowed(reader)

  my_count = 0
  while next_molecule(reader, m)
    my_count += 1
  end
  my_count == molecules_read(reader) || return is_failure("Incorrect count", m)
  molecules_read(reader) == (nsmiles - 1) || return is_failure("Wrong molecule count", m)
  connection_table_errors_encountered(reader) == 1 || return is_failure("Wrong error count", m)

  set_display_smiles_interpretation_error_messages(true)

  true
end

# Substructure Related
function test_query_from_smarts_ok()::Bool
  q = SubstructureQuery()
  build_from_smarts(q, "CC") || return is_failure("Bad smarts")
  true
end

function test_query_from_smarts_bad()::Bool
  q = SubstructureQuery()
  build_from_smarts(q, "J") && return is_failure("bad smarts ok")
  true
end

function test_query_from_textproto()::Bool
  proto_contents = raw"
query {
  smarts: \"Fc\"
}
"
  fname = write_smiles_tempfile("txtp", proto_contents)
  q = SubstructureQuery()
  read_proto(q, fname) || return is_failure("Cannot read $(fname)")

  m = Molecule()
  build_from_smiles(m, "Fc1cc(F)ccc1") || return is_failure("Bad smiles")
  matches(q, m) || return is_failure("query does not match")
  substructure_search(q, m) == 2 || return is_failure("Wrong match count", m)

  set_find_unique_embeddings_only!(q, true)
  substructure_search(q, m) == 2 || return is_failure("Wrong unique match count", m)

  set_perceive_symmetry_equivalent_matches!(q, false)
  substructure_search(q, m) == 1 || return is_failure("Wrong symmetric match count", m)
  true
end
  
function test_carbon_self_search()::Bool
  m = Molecule()
  build_from_smiles(m, "C methane") || return is_failure("bad smiles")
  q = SubstructureQuery()
  build_from_smarts(q, "C") || return is_failure("Cannot parse C smarts")
  substructure_search(q, m) == 1 || return is_failure("Not 1 match to C", m)
  q in m || return is_failure("C not in C", m)
  true
end

function test_ethane_c()::Bool
  m = Molecule()
  build_from_smiles(m, "CC ethane") || return is_failure("Bad smiles")
  q = SubstructureQuery()
  build_from_smarts(q, "C") || return is_failure("Bad smarts")
  substructure_search(q, m) == 2 || return is_failure("Not 2 matched", m)
  set_perceive_symmetry_equivalent_matches!(q, false)
  substructure_search(q, m) == 1 || return is_failure("Symmetry", m)
  true
end
  
function test_ethane_cc_embeddings_do_not_overlap()::Bool
  m = Molecule()
  build_from_smiles(m, "CC ethane") || return is_failure("Bad smiles")
  q = SubstructureQuery()
  build_from_smarts(q, "CC") || return is_failure("Bad smarts")
  substructure_search(q, m) == 2 || return is_failure("Not to C in CC", m)
  set_embeddings_can_overlap!(q, false)
  substructure_search(q, m) == 1 || return is_failure("Overlap", m)
  true
end
  
function test_ethane_cc_find_unique_embeddings_only()::Bool
  m = Molecule()
  build_from_smiles(m, "CC ethane") || return is_failure("Bad smiles")
  q = SubstructureQuery()
  build_from_smarts(q, "CC") || return is_failure("Bad smarts")
  substructure_search(q, m) == 2 || return is_failure("Not 2 matches", m)
  set_find_unique_embeddings_only!(q, true)
  substructure_search(q, m) == 1 || return is_failure("Bad count unique", m)
  true
end

  
function test_ethane_cc_find_one_embedding_per_atom()::Bool
  m = Molecule()
  build_from_smiles(m, "CC ethane") || return is_failure("bad smiles")
  q = SubstructureQuery()
  build_from_smarts(q, "CC") || return is_failure("Bad smarts")
  substructure_search(q, m) == 2 || return is_failure("Not 2 matches", m)
  set_find_one_embedding_per_atom!(q, true)
  substructure_search(q, m) == 2 || return is_failure("not 2 after", m)
  true
end

  

function test_results_returned()::Bool
  m = Molecule()
  build_from_smiles(m, "CC ethane") || return is_failure("Bad smiles")
  q = SubstructureQuery()
  build_from_smarts(q, "CC") || return is_failure("Bad smarts")
  sresults = SubstructureResults()
  substructure_search(q, m, sresults) == 2 || return is_failure("Not 2 matches", m)

  number_embeddings(sresults) == 2 || return is_failure("Not 2 embeddings", m)
  # Count the number of times each atoms is matched.
  matched = [0, 0]

  for embedding in embeddings(sresults)
    length(embedding) == 2 || return is_failure("Embedding not 2 atoms", m)
    for a in embedding
      matched[a + 1] += 1
    end
  end
  matched == [2, 2] || return is_failure("Not t matches each", m)
  true
end

function iterate_substructure_results()::Bool
  m = Molecule()
  build_from_smiles(m, "CCC") || return is_failure("Bad smiles")
  q = SubstructureQuery()
  build_from_smarts(q, "CC") || return is_failure("Bad smarts")
  sresults = SubstructureResults()
  substructure_search(q, m, sresults) == 4 || return is_failure("Not 4 embeddings", m)

  matches = embeddings(sresults)
  for i in eachindex(matches)
    e = matches[i]
    length(e) == 2 || return is_failure("Not 2 atoms in match", m)
  end

  true
end

function test_sresults_set_vector_all()::Bool
  m = Molecule()
  build_from_smiles(m, "Oc1c(O)c(O)c(O)c(O)c1O") || return is_failure("Bad smiles")
  q = SubstructureQuery()
  build_from_smarts(q, "Oc") || return is_failure("Bad smarts")
  sresults = SubstructureResults()
  substructure_search(q, m, sresults) == 6 || return is_failure("Not 6 matches", m)
  v = each_embedding_set_vector(sresults, natoms(m), 2)
  v == [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2] || return is_failure("Vector not match", m)

  true
end

function test_sresults_set_vector_partial()::Bool
  m = Molecule()
  build_from_smiles(m, "CNC") || return is_failure("Bad smiles")
  q = SubstructureQuery()
  build_from_smarts(q, "N") || return is_failure("Bad smarts")
  sresults = SubstructureResults()
  substructure_search(q, m, sresults) == 1 || return is_failure("Not 1 matches", m)
  v = each_embedding_set_vector(sresults, natoms(m), 3)
  v == [0, 3, 0] || return is_failure("Vector not match", m)

  true
end
  

function test_matched_atoms_returned()::Bool
  m = Molecule()
  build_from_smiles(m, "Oc1ccccc1") || return is_failure("Bad smiles")
  q = SubstructureQuery()
  build_from_smarts(q, "Occ") || return is_failure("Bad smarts")
  sresults = SubstructureResults()
  substructure_search(q, m, sresults)
  for embedding in embeddings(sresults)
    set_isotopes!(m, embedding, 1)
  end
  aromatic_smiles(m) == "[1OH][1c]1[1cH]ccc[1cH]1" || return is_failure("Not isotopically labelled", m)
  true
end

function test_substructure_search_as_vector()::Bool
  m = Molecule()

  build_from_smiles(m, "OC(=O)CC(=O)O") || return is_failure("Bad smiles")

  q = SubstructureQuery()
  build_from_smarts(q, "[OD1]-C=O") || return is_failure("Bad smarts")

  results = substructure_search_as_vector(q, m)
  size(results) == (2, 3) || return is_failure("Wrong shape", m)
  results[1,:] == [0, 1, 2] || return is_failure("First embedding", m)
  results[2,:] == [6, 4, 5] || return is_failure("First embedding", m)
  true
end

function test_reaction1()::Bool
  rxn_text = raw"
scaffold {
  smarts: \"[Pb]\"
  change_element {
    atom: 0
    element: \"Au\"
  }
}
"
  fname = write_smiles_tempfile("rxn1", rxn_text)
  rxn = Reaction()
  read_textproto(rxn, fname) || return is_failure("Cannot read reaction")

  m = Molecule()
  build_from_smiles(m, "[Pb]") || return is_failure("Bad smiles")
  in_place_transformation(rxn, m) == 1 || return is_failure("No reaction", m)
  smiles(m) == "[Au]" || return is_failure("Not gold", m)
  true
end

function test_reaction2()::Bool
  rxn_text = raw"
scaffold {
  smarts: \"Cl-c\"
  remove_atom: 0
}
"
  fname = write_smiles_tempfile("rxn2", rxn_text)
  rxn = Reaction()
  read_textproto(rxn, fname) || return is_failure("Cannot read reaction")

  m = Molecule()
  build_from_smiles(m, "Clc1ccc(Cl)cc1") || return is_failure("Bad smiles")
  in_place_transformation(rxn, m) == 1 || return is_failure("No reaction", m)
  aromatic_smiles(m) == "c1ccccc1" || return is_failure("Not removed", m)
  true
end

function test_reaction_single_reagent()::Bool
  rxn_text = raw"
scaffold {
  id: 0
  smarts: \"Cl-c\"
  remove_atom: 0
}
sidechain {
  id: 1
  reagent: \"O\"
  smarts: \"O\"
  join {
    a1: 1
    a2: 0
    btype: SS_SINGLE_BOND
  }
}
"
  fname = write_smiles_tempfile("rxn3", rxn_text)
  rxn = Reaction()
  read_textproto(rxn, fname) || return is_failure("Cannot read reaction")

  m = Molecule()
  build_from_smiles(m, "Clc1ccc(Cl)cc1") || return is_failure("Bad smiles")
  in_place_transformation(rxn, m) == 1 || return is_failure("No reaction", m)
  unique_smiles(m) == "Oc1ccc(O)cc1" || return is_failure("Not reacted", m)
  true
end

function test_reaction_iterator1()::Bool
  rxn_text = raw"
scaffold {
  id: 0
  smarts: \"[OH]-C(=O)-[#6]\"
  remove_atom: 0
}
sidechain {
  id: 1
  smarts: \"[NH2]-[CX4]\"
  join {
    a1: 1
    a2: 0
    btype: SS_SINGLE_BOND
  }
}
"
  fname = write_smiles_tempfile("rxn4", rxn_text)
  rxn = Reaction()
  read_textproto(rxn, fname) || return is_failure("Cannot read reaction")

  scaffold = Molecule()
  build_from_smiles(scaffold, "c1cc(C(=O)O)ccc1CC(=O)O benzoic acid") || return is_failure("Bad smiles")

  amine_smiles = [
    "CC(C)N ID:EN300-18989",
    "CCCN ID:EN300-20505",
    "NCC=C ID:EN300-19995",
    "NCC#C ID:EN300-21409",
    "NCCO ID:EN300-19392",
    "NC1CC1 ID:EN300-21353",
    "CC(C)(C)N ID:EN300-16766",
    "CCCCN ID:EN300-19577",
    "CCC(C)N ID:EN300-19017"
  ]

  amines = [LillyMol.MolFromSmiles(s) for s in amine_smiles]

  smc = SidechainMatchConditions()
  [add_sidechain(rxn, 0, m, smc) for m in amines]
   
  sresults = SubstructureResults()
  substructure_search(rxn, scaffold, sresults) == 2 || return is_failure("Not 2 acids", scaffold)

  seen = Set{String}()
  for embedding in embeddings(sresults)
    iter = ReactionIterator()
    initialise(iter, rxn)
    # iter = ReactionIterator(rxn)

    while active(iter)
      product = Molecule()
      perform_reaction(rxn, scaffold, embedding, iter, product) || return is_failure("Cannot perform reaction", scaffold) || break
      usmi = unique_smiles(product)
      usmi in seen && return is_failure("Duplicate product", product)
      push!(seen, usmi)
      increment!(iter)
    end
  end

  length(seen) == 18 || return is_failure("Not 18 products", scaffold)

  true
end


boobar()
@test test_empty_molecule()
@test test_methane()
@test test_set_of_atoms()
@test test_set_of_atoms_equals()
@test test_build_from_smiles()
@test test_copy_constructor()
@test test_hydrogen_related()
@test test_copy_operation()
@test test_copy_constructor_with_name()
@test test_isotopes()
@test test_set_name()
@test test_natoms()
@test test_natoms_atomic_number()
@test test_natoms_element()
@test test_length_molecule()
@test test_atomic_number()
@test test_atomic_number_in()
@test test_enumerate_atom()
@test test_nedges()
@test test_molecular_formula()
@test test_valence_ok()
@test test_atom_valence_ok()
@test test_standardise()
@test test_nrings()
@test test_nrings_atom()
@test test_nrings_size()
@test test_is_ring_atom()
@test test_fused_system_size()
@test test_fused_system_identifier()
#@test test_fused_to()
@test test_rings_with_fused_system_identifier()
@test test_in_same_ring()
@test test_in_same_aromatic_ring()
@test test_in_same_ring_system()
@test test_ring_membership()
@test test_rings_containing_both()
@test test_is_part_of_fused_ring_system()
@test test_ring()
@test test_rings()
@test test_ring_containing_atom()
@test test_ring_related()
@test test_label_atoms_by_ring_system()
@test test_label_atoms_by_ring_system_including_spiro_fused()
@test test_nrings_including_non_sssr_rings()
@test test_non_sssr_rings()
@test test_non_sssr_ring()
@test test_is_spiro_fused()
@test test_is_halogen()
@test test_ncon_molecule()
@test test_nbonds_molecule()
@test test_maximum_connectivity()
@test test_connections_molecule()
@test test_isotopically_labelled_smiles()
@test test_isotopically_labelled_smiles2()
@test test_is_aromatic()
@test test_getindex_molecule()
@test test_number_formally_charged_atoms()
@test test_net_formal_charge()
@test test_bond_molecule()
@test test_bond_between_atoms()
@test test_number_symmetry_classes()
@test test_symmetry_class()
@test test_symmetry_equivalents() broken=true
@test test_symmetry_classes()
@test test_symmetry()
@test test_attached_heteroatom_count()
@test test_bond_length()
@test test_bond_angle()
#@test test_dihedral_angle()
#@test test_signed_dihedral_angle()
#@test test_set_bond_length()
#@test test_set_bond_angle()
#@test test_set_dihedral_angle()
#@test test_xyx_molecule()
#@test test_distance_between_atoms()
#@test tests_longest_intra_molecular_distance()
@test test_add_bond()
@test test_set_bond_type_between_atoms()
@test test_set_atomic_number()
@test test_are_bonded()
@test test_bond_between_atoms()
@test test_bond_list()
@test test_add_molecule()
@test test_remove_atom()
@test test_remove_atom2()
@test test_remove_atoms()
@test test_remove_atoms_vector()
@test test_delete_fragment()
@test test_delete_fragment2()
@test test_remove_fragment_containing_atom()
@test test_remove_all()
@test test_set_auto_create_new_elements()
@test test_set_atomic_symbols_can_have_arbitrary_length()
@test test_set_atomic_symbols_can_have_arbitrary_length2()
@test test_remove_all_non_natural_elements()
@test test_remove_explicit_hydrogens()
@test test_chop()
@test test_remove_bonds_to_atom()
@test test_remove_bond()
@test test_remove_bond_between_atoms()
@test test_remove_all_bonds()
@test test_molecular_weight()
@test test_molecular_weight_count_isotopes()
@test test_molecular_weight_ignore_isotopes()
@test test_highest_coordinate_dimensionality()
@test test_exact_mass()
@test test_number_fragments()
@test test_fragment_membership()
@test test_fragment_membership2()
@test test_fragment_membership_vector()
@test test_atoms_in_fragment()
@test test_get_atoms_in_fragment()
@test test_largest_fragment()
@test test_largest_fragment_carefully()
@test test_largest_fragment2()
@test test_identify_spinach()
@test test_rings_in_fragment()
@test test_create_components()
@test test_returns_vector()
@test test_create_subset()
@test test_create_subset_set_of_atoms()
@test test_reduce_to_largest_fragment()
@test test_reduce_to_largest_organic_fragment()
@test test_reduce_to_largest_fragment_carefully()
@test test_fragment_related()
@test test_atoms_in_fragment2()
@test test_organic_only()
@test test_contains_non_periodic_table_elements()
@test test_longest_path()
@test test_bonds_between()
@test test_distance_matrix()
@test test_atoms_between()
@test test_implicit_hydrogens()
@test test_explicit_hydrogens()
@test test_hcount()
@test test_make_implicit_hydrogens_explicit()
@test test_move_hydrogens_to_end_of_connection_table()
@test test_pi_electrons()
@test test_lone_pair_count()
@test test_aromatic_atom_count()
@test test_aromatic_ring_count()
@test test_saturated()
@test test_number_chiral_centres()
@test test_chiral_centre_at_atom()
@test test_chiral_centre_in_molecule_not_indexed_by_atom_number()
@test test_remove_chiral_centre_at_atom()
@test test_remove_all_chiral_centres()
@test test_invert_chirality_on_atom()
@test test_smarts_equivalent_for_atom()
@test test_smarts()
@test test_atom_map_number()
@test test_atom_map_numbers()
@test test_set_atom_map_number()
@test test_set_include_atom_map_with_smiles()
@test test_atom_with_atom_map_number()
@test test_reset_all_atom_map_numbers()
@test test_unset_unnecessary_implicit_hydrogens_known_values()
@test test_discern_chirality_from_3d_structure()
@test test_set_formal_charge()
@test test_sort_atoms()
@test test_random_smiles()
@test test_smiles_starting_atom()
@test test_atom_iterator()
@test test_iterate_bonds()
@test test_index_bond_list()
@test test_loop_bond_list()
@test build_benzene()
@test test_imidazole()
@test test_find_exocyclic_bond()
@test test_atom_iterator_and_valence()
@test test_index_atom()
@test test_scaffold()
@test test_gather_rings()
@test iterate_ring_information()
@test test_coords()
@test test_getxyz()
@test test_setxyz()
@test test_dihedral_scan() skip=true
@test test_rule_of_five()
@test test_number_chiral_smiles1()
@test test_number_chiral_smiles2()
@test test_chiral_implicit_hydrogen()
@test test_invert_chirality()
@test test_iterate_chiral_centres()
@test test_index_chiral_centres()
@test test_charge()
@test test_xlogp()
@test test_aspirin()
@test test_cubane()

@test test_read_smiles()
@test test_read_smiles_with_errors_fail()
@test test_read_smiles_with_errors_skip()

# Substructure related
@test test_query_from_smarts_ok()
@test test_query_from_smarts_bad()
@test test_query_from_textproto()
@test test_carbon_self_search()
@test test_ethane_c()
@test test_ethane_cc_embeddings_do_not_overlap()
@test test_ethane_cc_find_unique_embeddings_only()
@test test_ethane_cc_find_one_embedding_per_atom()
@test test_results_returned()
@test iterate_substructure_results()
@test test_sresults_set_vector_all()
@test test_sresults_set_vector_partial()
@test test_matched_atoms_returned()
@test test_substructure_search_as_vector()

@test test_reaction1()
@test test_reaction2()
@test test_reaction_single_reagent()
@test test_reaction_iterator1()

