# Pharmacophore_2d

## Objective
Given a set of active molecules, generate a set of substructure queries that describe the
pharmacophoric features in that set of molecules.

Generally the tool will be driven by a set of queries which define the pharmacophore. As those
pharmacophoric atoms are identified in the molecules, the query is constructed to preserve
the relative separation of those atoms. The atomic characteristics of the matched atoms can
also be specified, so there is considerable flexibility on how precisely the atom matching is
done.

For example, if 3-aminopropanol, `NCCCO`, were an active molecule, we might define the pharmacophoric
elements to be a primary amine and an alcohol. In the starting molecules the pharmacophores are
separated by 3 bonds.

An exact match to this query might look like
```
name: "1aminopropanol"
query {
  respect_initial_atom_numbering: true
  query_atom {
    id: 0
    atom_properties {
      atomic_number: 7
    }
  }
  query_atom {
    id: 4
    atom_properties {
      atomic_number: 8
    }
  }
  separated_atoms {
    a1: 4
    a2: 0
    bonds_between: 4
  }
}
```
This query will match any N...O grouping. But note that there are no constraints on
what kind of Nitrogen or Oxygen atom (alpiphatic only). This can be controlled by the textproto config
file that controls execution. If we use
```
functional_group: "SMARTS:[ND1H2]"
functional_group: "SMARTS:[OD1H]"
all_functional_groups_must_match: true

atomic_property: ATOMIC_NUMBER
atomic_property: HCOUNT
```
as the configuration file, the resulting query file will be
```
name: "1aminopropanol"
query {
  respect_initial_atom_numbering: true
  query_atom {
    id: 0
    atom_properties {
      atomic_number: 7
      hcount: 2
    }
  }
  query_atom {
    id: 4
    atom_properties {
      atomic_number: 8
      hcount: 1
    }
  }
  separated_atoms {
    a1: 4
    a2: 0
    bonds_between: 4
  }
}
```
where we see that the Hydrogen count has been added to the atomic properties
that must be matched. Note that this information did *not* derive from the smarts.
The smarts matched the atoms, and the hydrogen count is derived from the matched
atoms. If would probably not make sense to specify a highly specific smarts, and then
only use `atomic_number` as the query property.

Atomic properties from the enum in the proto definition can be used.
```
enum AtomicProperty {
  UNSPECIFIED = 0;
  ATOMIC_NUMBER = 1;
  NCON = 2;
  AP_AROMATIC = 3;
  RING_BOND_COUNT = 4;  
  HAS_PI_ELECTRON = 5;  // Not implemented, do not use
  PI_ELECTRON_COUNT = 6;  // Not implemented, do not use
  UNSATURATION = 7;
  ISOTOPE = 8;
  RING_SIZE = 9;
  SPINACH = 10;
  FUSED_SYSTEM_SIZE = 11;
  HCOUNT = 12;
}
```
Any number of these properties can be combined. All properties will be transferred
from the matched atom to the query file, and so all properties must be matched.

## Fuzzy matching.
In the above example, the match and the `separated_atoms` directives are exact. We can
obtain less precise matches by turning on directives in the configuraiton file.

```
functional_group: "SMARTS:[ND1H2]"
functional_group: "SMARTS:[OD1H]"
all_functional_groups_must_match: true

min_separation: 3
max_separation: 8
delta_longer: 4
ncon_becomes_min_ncon: true
hcount_becomes_min: true

atomic_property: ATOMIC_NUMBER
atomic_property: HCOUNT
atomic_property: NCON
```
generates the query file
```
name: "1aminopropanol"
query {
  respect_initial_atom_numbering: true
  query_atom {
    id: 0
    atom_properties {
      atomic_number: 7
      min_ncon: 1
      min_hcount: 2
    }
  }
  query_atom {
    id: 4
    atom_properties {
      atomic_number: 8
      min_ncon: 1
      min_hcount: 1
    }
  }
  separated_atoms {
    a1: 4
    a2: 0
    max_bonds_between: 8
  }
}
```
where now some of the atomic properties are specified as minimum rather than specific
values. And bonds_between is now max_bonds_between, because of the `delta_longer` setting.

The proto definition enables setting of several attributes that enable conversion of
specific values in the molecules to min/max values in the query proto.

As a check, the query generated from a molecule should match the starting molecule.

A more complex example might involve a pharmacophore defined by an aromatic nitrogen
atom and an amide group. In addition, this shows how the number of rotatable bonds
between pharmacophoric atoms can be specified.
```
functional_group: "SMARTS:[/IWfss2nr5]"
functional_group: "SMARTS:O=C-[NR0]"
all_functional_groups_must_match: true

min_separation: 5
max_separation: 10

atomic_property: ATOMIC_NUMBER
atomic_property: AP_AROMATIC
atomic_property: SPINACH
atomic_property: FUSED_SYSTEM_SIZE

extra_rotbond: 3
```
which might generate a query file such as
```
name: "CHEMBL1471636"
query {
  respect_initial_atom_numbering: true
  query_atom {
    id: 7
    atom_properties {
      atomic_number: 7
      aromatic: true
      fused_system_size: 2
      match_spinach_only: 0
    }
  }
  query_atom {
    id: 10
    atom_properties {
      atomic_number: 7
      match_spinach_only: 1
    }
  }
  query_atom {
    id: 11
    atom_properties {
      atomic_number: 6
      match_spinach_only: 1
    }
    query_bond {
      btype: SS_SINGLE_BOND
      other_end: 10
    }
  }
  query_atom {
    id: 12
    atom_properties {
      atomic_number: 8
      match_spinach_only: 1
    }
    query_bond {
      btype: SS_DOUBLE_BOND
      other_end: 11
    }
  }
  separated_atoms {
    a1: 12
    a2: 7
    bonds_between: 5
    max_rotbond: 5
  }
}
```
### Example
If we were to look for pharmacaphore equivalents for Dasatinib, that input file
might look like
```
functional_group: "SMARTS:[OD1]CCN1CCNCC1"
functional_group: "SMARTS:O=C-Nc1c(Cl)cccc1C"
delta_shorter: 2
delta_longer: 7
atomic_property: ATOMIC_NUMBER
atomic_property: NCON
atomic_property: AP_AROMATIC
atomic_property: RING_BOND_COUNT
atomic_property: UNSATURATION
atomic_property: RING_SIZE
```

