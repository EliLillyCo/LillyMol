# Substructure Queries as Protos

LillyMol has two ways to describe substructure queries beyond command line
SMARTS:

- the original LillyMol query-file format, which still works;
- protocol buffer query files, which are the preferred format for new work.

All new substructure features are exposed only through protocol buffers. The reason
for having query files at all is that some chemically useful ideas are awkward,
fragile, or effectively impossible to express as plain SMARTS. Protocol buffers
make those ideas explicit and maintainable.

The schema is the authority. If this document and the schema disagree, trust
[substructure.proto](/src/Molecule_Lib/substructure.proto). The main message is
`SubstructureSearch::SubstructureQuery`, usually written as textproto.

For background on why LillyMol uses protocol buffers, see
[why_protos.md](why_protos.md). For SMARTS syntax and LillyMol SMARTS
extensions, see [substructure.md](substructure.md).

## Quick Start

A minimal textproto query wraps a SMARTS in a `query` block:

```textproto
name: "aromatic nitrogen"
query {
  smarts: "n"
}
```

Run it with `tsubstructure`:

```shell
tsubstructure -q PROTO:aromatic_nitrogen.textproto molecules.smi
```
Recently LillyMol query reading was updated so that the `PROTO:` prefix is no longer
always needed. The first few lines of the file are examined, and if the file looks
like textproto form, it is interpreted as textproto. For certainty and to avoid
any ambiguity, use `PROTO:`. For convenience, omit `PROTO:`.

A file containing many query-file names can be supplied with `PROTOFILE`:

```shell
ls *.textproto > all_queries
tsubstructure -q PROTOFILE:all_queries molecules.smi
```

For very large query collections, serialized
`SubstructureSearch::SubstructureQuery` protos can be read from TFDataRecord
files:

```shell
tsubstructure -q TFPROTO:queries.tfrecord molecules.smi
```

A textproto can also be supplied directly on the command line, although this is
usually only useful for small tests:

```shell
tsubstructure -q 'proto:query { smarts: "n" }' molecules.smi
```

## Drafting A Query With Codex

Tools like Codex or Claude can be useful for drafting proto query files,
especially when the desired query can first be stated in ordinary chemical
language. Treat the generated query as a starting point, not as an authority.
Always inspect the SMARTS or proto fields, then test the query against small
positive and negative examples.

### Simple Example: Phenolic Oxygen

**User**

Write a LillyMol substructure proto query that matches phenolic oxygen.

**Codex**

A phenolic oxygen is an oxygen directly attached to an aromatic carbon. A minimal
query is:

```textproto
name: "phenolic oxygen"
query {
  smarts: "[OD1]-c"
}
```

`[OD1]-c` matches a singly connected oxygen attached to an aromatic carbon. If
the query should require a hydrogen on the oxygen, use `[OD1H1]-c`.

Save the query as `phenolic_oxygen.textproto` and test it with a few small
examples:

```shell
cat > phenol_examples.smi <<'EOF'
Oc1ccccc1 phenol
COc1ccccc1 anisole
CCO ethanol
O=C(O)c1ccccc1 benzoic_acid
EOF

tsubstructure -q PROTO:phenolic_oxygen.textproto phenol_examples.smi
```

The expected match is `phenol`. `anisole` should not match because its oxygen has
two heavy-atom connections, `ethanol` should not match because the carbon is not
aromatic, and `benzoic_acid` should not match because the oxygen is attached to a
carbonyl carbon rather than an aromatic carbon.

To help visualize the matched atoms in a more complex query, write the matched
molecules and place isotopic labels on the matched atoms.

```shell
tsubstructure -m matches -j 1 -q PROTO:phenolic_oxygen ...
```

### More Complex Example: Aniline With An Electron-Withdrawing Group

**User**

I want a query for anilines where the aromatic ring also has an
electron-withdrawing group. The group should be on the same isolated aromatic
ring as the aniline nitrogen. Examples of electron-withdrawing groups are nitro,
cyano, carbonyl, fluoro, chloro, bromo, and iodo.

**Codex**

The local aniline part is simple SMARTS:

```smarts
[ND1H2]-c
```

The harder part is making sure the electron-withdrawing group is on the same
ring. If we tried to do this as ordinary SMARTS, we would need to encode ring
position and each allowed substituent pattern directly. A narrow para-nitro
example is readable:

```smarts
[ND1H2]-c1ccc([N+](=O)[O-])cc1
```

But that only covers one substituent at one ring position. A broader SMARTS-only
version quickly becomes difficult to maintain because ortho, meta, and para
positions, and all allowed electron-withdrawing groups, must be enumerated.

A LillyMol SMARTS expression could use the LillyMol `||` extension to combine
alternatives, but `||` is not standard SMARTS. For a portable SMARTS-only query,
those alternatives would need to be expanded explicitly. In a LillyMol proto
query, we can instead separate the concepts:

- find an isolated aromatic ring with one of the electron-withdrawing
  environments;
- assign a global id to that ring;
- require the aniline nitrogen to be attached to an atom on that same ring.

```textproto
name: "aniline same isolated aromatic ring as EWG"
query {
  ring_specifier {
    aromatic: true
    fused: 0
    base {
      environment: "a-[N+](=O)[O-]||a-C#N||a-C(=O)[#6,O,N]||a-[F,Cl,Br,I]"
      set_global_id: 1
    }
  }

  smarts: "[ND1H2]-[/IWgid1a]"
}
```

`environment` describes what must be attached somewhere to the ring. The
`set_global_id: 1` directive labels the atoms in the ring that satisfied the ring
specifier. `/IWgid1` in the SMARTS then makes the aniline match use an atom from
that same ring.

**User**

Does the environment become part of the returned match? I only want the aniline
N and the ring atom in the embedding.

**Codex**

No. In this query, the ring environment is a requirement on the ring, but the
environment atoms are not returned as matched atoms. The returned embedding comes
from the final SMARTS:

```textproto
smarts: "[ND1H2]-[/IWgid1a]"
```

This is useful for reactions or matched-atom reporting because the proto can use
rich context to decide whether a molecule matches, while still returning only the
atoms specified by the main SMARTS.

Test it with simple positive and negative examples:

```shell
cat > aniline_ewg_examples.smi <<'EOF'
Nc1ccc([N+](=O)[O-])cc1 p_nitroaniline
Nc1ccccc1 aniline
Nc1ccc(C#N)cc1 p_cyanoaniline
Nc1ccccc1c1ccc([N+](=O)[O-])cc1 different_ring_same_molecule
CCN ethylamine
EOF

tsubstructure -q PROTO:aniline_same_ring_ewg.textproto aniline_ewg_examples.smi
```

The expected matches are `p_nitroaniline` and `p_cyanoaniline`. Plain `aniline`
lacks the electron-withdrawing environment, `different_ring_same_molecule` contains
an aniline and a nitro-substituted ring but not on the same ring, and
`ethylamine` lacks the aromatic ring.

### Complex Example: Meta Substitution On A Constrained Ring

**User**

I need a query for a six-membered aromatic ring where fluorine and hydroxy are
meta to each other. A simple SMARTS for the positional relationship might be
`F-c:c:c-[OH]`, but I also want the ring to have no other halogen substituent.
How would I express that as a proto query?

**Codex**

Start with the part that really is local SMARTS. We want to return the fluorine,
the fluorinated ring atom, the hydroxy oxygen, and the hydroxy-bearing ring atom:

```smarts
F-c.[OH]-c
```

The two ring atoms are matched atoms 1 and 3 in that SMARTS. On a six-membered
aromatic ring, meta substitution means those two ring atoms are separated by two
bonds around the ring. We can express that with `separated_atoms` instead of
making the ring path part of the SMARTS:

```textproto
name: "meta fluoro phenol"
query {
  smarts: "F-c.[OH]-c"
  separated_atoms {
    a1: 1
    a2: 3
    bonds_between: 2
  }
}
```

This is already easier to maintain than spelling out every atom in the ring, but
it still does not constrain the ring itself. To add ring requirements, define the
ring independently and give it a global id. Then require the two matched aromatic
atoms to be on that ring.

```textproto
name: "meta fluoro phenol, only two ring substituents"
query {
  ring_specifier {
    ring_size: 6
    aromatic: true
    fused: 0
    base {
      ncon: 2
      set_global_id: 1
    }
  }

  smarts: "F-[/IWgid1a].[OH]-[/IWgid1a]"

  separated_atoms {
    a1: 1
    a2: 3
    bonds_between: 2
  }
}
```

Here the ring specifier says there must be an isolated six-membered aromatic
ring with exactly two exocyclic connections. The SMARTS then requires both
substituted ring atoms to come from that same ring. Because `ncon: 2` allows only
two connections from the ring to atoms outside the ring, this excludes any
additional halogen substituent. It also excludes any other extra substituent,
such as methyl or nitro. If other non-halogen substituents should be allowed,
the query would need a more specific negative condition rather than the simple
`ncon: 2` constraint.

**User**

So `F-c:c:c-[OH]` is still a useful way to explain the idea, but the proto query
is separating the positional relationship from the ring constraints?

**Codex**

Yes. The SMARTS names the two substituents and their attachment atoms. The
`separated_atoms` block says those attachment atoms must be meta. The
`ring_specifier` block says what kind of ring they must be attached to. That
division is the main advantage of a query file: each chemical idea is visible as
its own condition.

Test it with:

```shell
cat > meta_fluoro_phenol_examples.smi <<'EOF'
Oc1cccc(F)c1 m_fluorophenol
Oc1ccc(F)cc1 p_fluorophenol
Oc1c(F)cccc1 o_fluorophenol
Oc1cccc(F)c1Cl extra_chloro
Oc1cccc(F)c1C extra_methyl
EOF

tsubstructure -q PROTO:meta_fluoro_phenol.textproto meta_fluoro_phenol_examples.smi
```

The expected match is `m_fluorophenol`. The ortho and para examples fail the
bond-separation constraint. `extra_chloro` and `extra_methyl` fail because the
ring has more than two exocyclic connections.

## Mental Model

A substructure proto has two levels.

`SubstructureQuery` is the outer container. It has a `name`, optional comments,
one or more `query` blocks, and optional logical operators connecting those
query blocks.

`SingleSubstructureQuery` is each `query` block. This is where the actual SMARTS,
ring specifications, environments, hit-count requirements, and other matching
constraints live.

For example:

```textproto
name: "carbon or nitrogen"
query {
  smarts: "C"
}
query {
  smarts: "N"
}
```

With two `query` blocks and no `logexp`, the default relationship is OR. This is
like the SMARTS expression `C||N`.

To require both components, use `logexp`:

```textproto
name: "carbon and nitrogen"
query {
  smarts: "C"
}
query {
  smarts: "N"
}
logexp: SS_AND
```

With three components, provide two operators:

```textproto
query {
  smarts: "C"
}
query {
  smarts: "N"
}
query {
  smarts: "F"
}
logexp: [SS_AND, SS_OR]
```

The `smarts` field inside one `query` block must not contain top-level `&&` or
`||` operators. Put each component in its own `query` block and connect them with
`logexp`.

Avoid making a single huge composite query when you really have many independent
queries. Separate query files are easier to test, easier to name, and usually
easier to maintain.

## When To Use A Proto

Use plain SMARTS when the query is naturally expressed as atom and bond
relationships:

```shell
tsubstructure -s 'Nc1ccccc1 aniline' molecules.smi
```

Use a proto when you need one or more of these:

| Need | Proto construct |
| ---- | --------------- |
| Whole-molecule prefilter before atom matching | `required_molecular_properties` |
| Ring constraints not tied to a specific SMARTS atom | `ring_specifier` |
| Ring-system constraints | `ring_system_specifier` |
| Link a ring/ring-system prefilter to atoms in the returned embedding | `set_global_id` plus `/IWgidN` |
| Require context next to matched atoms | `environment` |
| Reject context next to matched atoms | `environment_no_match` |
| Constrain a substituent attached to a ring or matched atom | `substituent` |
| Require a matched atom to satisfy one of several rooted SMARTS | `matched_atom_must_be` |
| Require or reject a motif within a bond range of matched atoms | `nearby_atoms` |
| Constrain the path or region between matched atoms | `link_atoms`, `separated_atoms`, `region`, `inter_ring_region` |
| Need a complex down the bond directive | `down_the_bond` |

## Whole-Molecule Prefilters

`required_molecular_properties` checks simple molecular properties before atom
matching starts. It is mainly an expressive and sometimes efficient way to say
"do not even try this query unless the molecule has the required global shape".

```textproto
name: "oxygen rich aromatic molecules"
query {
  required_molecular_properties {
    min_natoms: 20
    min_aromatic_rings: 2
    elements_needed {
      atomic_number: 8
      min_hits_needed: 2
    }
    elements_needed {
      atomic_symbol: "P"
      hits_needed: 0
    }
  }
  smarts: "o"
}
```

These prefilters are independent of the later atom match. In the example above,
the oxygen atoms counted by `elements_needed` are not necessarily the atom matched
by `smarts: "o"`.

Commonly used whole-molecule properties include:

- `natoms`
- `nrings`
- `heteroatoms_in_molecule`
- `ring_systems`
- `aromatic_rings`
- `aromatic_atoms`
- `non_aromatic_rings`
- `number_isotopic_atoms`
- `number_fragments`
- `atoms_in_spinach`
- `inter_ring_atoms`
- `net_formal_charge`
- `any_net_formal_charge`
- `elements_needed`
- `required_bond`

Most count properties support exact, minimum, and maximum forms such as
`nrings`, `min_nrings`, and `max_nrings`. See the schema for the complete list.

### Required Bonds

`required_bond` checks whether a bond exists anywhere in the molecule. It does
not identify a bond that must participate in the later atom embedding.

```textproto
query {
  required_molecular_properties {
    required_bond {
      atomic_number_1: 6
      btype: SS_SINGLE_BOND
      atomic_number_2: 8
    }
  }
  smarts: "N"
}
```

This says "only search molecules that contain a carbon-oxygen single bond, and
then look for nitrogen".

## Ring And Ring-System Specifications

Ring and ring-system specifiers are among the most useful proto features. They
allow queries such as:

- an isolated aromatic six-membered ring;
- a fused 6-6 aromatic ring system;
- a ring system with one heteroatom and two substituent groups;
- a ring with a particular environment attached somewhere;
- logical combinations of ring or ring-system environment requirements;
- the total connectivity to a ring or ring system;

A ring specifier searches for a ring satisfying the conditions. A ring-system
specifier searches for a fused or isolated ring system satisfying the conditions.
By themselves, these constructs are preconditions. They do not automatically say
that subsequent SMARTS within a query block match must occur in that same ring.

```textproto
name: "aniline and an EWG somewhere"
query {
  ring_specifier {
    aromatic: true
    fused: 0
    base {
      environment: "a-[F,Cl]||a-N(=O)=O||a-C=O"
    }
  }
  smarts: "[ar6]-[ND1H2]"
}
```

This requires an unfused aromatic ring with an electron-withdrawing group, and it
also requires an aniline-like nitrogen. But those two conditions are independent:
the aniline nitrogen might be attached to a different ring.

## Global IDs

Use global ids when atoms found by a ring, ring-system, substituent, or other
non-SMARTS construct must be the same atoms used by the SMARTS match.

```textproto
name: "aniline EWG on the same ring"
query {
  ring_specifier {
    aromatic: true
    fused: 0
    base {
      environment: "a-[F,Cl]||a-N(=O)=O||a-C=O"
      set_global_id: 5
    }
  }
  smarts: "[/IWgid5a]-[ND1H2]"
}
```

`set_global_id: 5` labels all atoms in the ring that satisfied the ring
specifier. `/IWgid5` in the SMARTS then requires the matched aromatic atom to be
one of those labelled atoms. Same for a ring system specifier.

A negative global id in SMARTS, such as `/IWgid-3`, means do NOT match any atom
having global id 3. The form `/IWgid!3` is also recognised.

This is the single most important rule for proto ring queries: if a ring or
ring-system condition must be tied to the returned atom embedding, use a global
id.

## Ring Base Conditions

`ring_specifier` and `ring_system_specifier` both contain a `base` message. The
base describes properties common to rings and ring systems, such as environment,
substitution, and global-id assignment.

```textproto
name: "disubstituted aniline EWG"
query {
  ring_specifier {
    aromatic: true
    fused: 0
    base {
      environment: "a-[F,Cl]||a-N(=O)=O||a-C=O"
      ncon: 2
      set_global_id: 5
    }
  }
  smarts: "[/IWgid5a]-[ND1H2]"
}
```

This matches an unfused aromatic ring with exactly two connections to atoms
outside the ring, one of which is the aniline nitrogen matched by the SMARTS.

See the `SubstructureRingBase` message in the proto definition file
[substructure.proto](/src/Molecule_Lib/substructure.proto) for information
about what matching attributes are available in the base ring specification.

## Ring Systems

Use `ring_system_specifier` when the chemically meaningful object is the whole
fused system rather than one SSSR ring.

```textproto
name: "fused 6-6 aromatic system"
query {
  ring_system_specifier {
    rings_in_system: 2
    aromatic_ring_count: 2
    ring_size_requirement {
      ring_size: 6
      count: 2
    }
    base {
      set_global_id: 1
    }
  }
  smarts: "[/IWgid1]-O-[CH3]"
}
```
with an attached methoxy. Note that the attached methoxy could have been specified as
an environment for the ring system, but then the matched atoms would not have been
available via the query match.

Useful ring-system and ring-system `base` fields include:

- `rings_in_system`, `min_rings_in_system`, `max_rings_in_system`
- `ring_size_requirement`
- `aromatic_ring_count`
- `non_aromatic_ring_count`
- `degree_of_fusion`
- `atoms_in_system`
- `number_spinach_groups`
- `number_non_spinach_groups`
- `atoms_in_spinach_group`
- `length_of_spinach_group`
- `distance_to_another_ring`
- `strongly_fused_ring_count`
- `ring_systems_extend_across_spiro_fusions`
- `environment_sets_global_id`

By default, atoms matched by a ring or ring-system environment do not set any
associated global id; only the atoms that are part of the ring or ring system are
assigned global ids. That can be changed with `environment_sets_global_id: true`.
Note however that if the environment is an OR condition, only the atoms associated
with the first matching condition will be marked. If it is an AND condition, the
environment atoms associated with the matching conditions are marked.

Setting `number_non_spinach_groups: 1` identifies a terminal ring system,
setting a larger number implies an "internal" ring.

## Substituents

A `substituent` describes a group attached to a ring, ring system, or matched
atom. It can constrain size, ring count, length from the attachment point,
heteroatom count, unsaturation, and required or disqualifying SMARTS.

```textproto
name: "fused aromatic system with amide substituent"
query {
  ring_system_specifier {
    rings_in_system: 2
    aromatic_ring_count: 2
    ring_size_requirement {
      ring_size: 6
      count: 2
    }
    base {
      set_global_id: 1
      max_ncon: 4
      substituent {
        min_natoms: 3
        max_natoms: 12
        max_nrings: 1
        hits_needed: 1
        max_length: 8
        required_smarts: "[ND2H]-C(=O)[#6]"
        disqualifying_smarts: "[Cl,Br,I,OH]||O=N=O"
        set_global_id: 2
      }
    }
  }

  smarts: "[/IWgid1]!@[/IWgid2D2]"
}
```

The ring system atoms receive global id 1. The substituent receives global id 2.
The final SMARTS requires a non-ring bond between the ring system and the matched
substituent atom.

`substituent_no_match` is the corresponding rejection form.

## Environments
The concept of environment was introduced in with ring and ring systems
specifications: atoms that sprout from the atoms in the ring system.
The environment concept described here relates to the atoms matched by
the overall smarts of the query.

An environment is a small query attached to an already matched atom. The
environment atoms are used to accept or reject the embedding, but they are not
returned as matched atoms.

Use `environment` for required context and `environment_no_match` for rejecting
context.

```textproto
query {
  smarts: "[ND<4G0T0]"
  environment_no_match {
    attachment {
      attachment_point: 0
      btype: SS_SINGLE_BOND
    }
    smarts: "a"
    smarts: "[G>0]"
    smarts: "[+]"
    smarts: "C-[+]"
    smarts: "C-C-[+]"
    smarts: "C-[O,F,P,S,Cl,Br,I]"
  }
}
```

By default, environment atoms must match atoms that are not already part of the
main embedding. That can be changed with
`environment_must_match_unmatched_atoms`. This is different from the environment
concept within a ring or ring system. In that case, the first matched atom in
the environment must be part of the ring. Here the first atom of the environment
must be an unmatched atom attached to a matched atom.

`attachment_point` is repeated, so one environment can be allowed at multiple
matched atoms:

```textproto
query {
  smarts: "[ND1H2][/IWfsid1ar6]:a:a:a"
  environment {
    attachment {
      attachment_point: [2, 4]
      btype: SS_SINGLE_BOND
    }
    smarts: "[Cl,Br,I]"
  }
}
```

If multiple environments are specified, they must all match unless `and_id` or
`or_id` is used. Those fields exist for more complex logic, but they are harder
to read than separate well-named query files.

## Matched-Atom Constraints

`matched_atom_must_be` applies one or more rooted SMARTS to a specific matched
atom. The first atom in each SMARTS is matched against the already matched atom.
This is often easier to read than recursive SMARTS.

```textproto
query {
  smarts: "[Cx0]"
  matched_atom_must_be {
    atom: 0
    smarts: "[CH3]"
    smarts: "[CD2]-[OH]"
    smarts: "[CD2]-C(F)(F)F"
    smarts: "[CD2]-O-[CD1]"
  }
}
```

Positive SMARTS entries are ORed: at least one must match. Entries beginning
with `!` are rejections: if any rejection matches, the embedding fails.

```textproto
query {
  smarts: "[ND1H2]"
  matched_atom_must_be {
    atom: 0
    smarts: "!N-C=O"
    smarts: "!N-S=O"
    smarts: "!N-C=S"
    smarts: "!N-a"
    smarts: "!N-[!#6]"
    smarts: "!N...{<3}N"
  }
}
```

## Nearby Atoms

`nearby_atoms` describes motifs that must, or must not, occur within a bond
range of matched atoms. It is useful for less specific lead optimisation ideas
such as "a donor within three bonds of this ring system".

```textproto
query {
  ring_system_specifier {
    rings_in_system: 2
    aromatic_ring_count: 2
    base {
      min_heteroatom_count: 1
      set_global_id: 1
      environment: "[o,nD2H0]&&[cD3]-!@*"
    }
  }

  smarts: "[/IWgid1]"

  nearby_atoms {
    smarts: "[OH,SH]"
    max_bonds_between: 3
    hits_needed: 1
  }
  nearby_atoms {
    smarts: "[ND1H2T0]-[G0]"
    max_bonds_between: 2
    rejection: true
  }
}
```

This matches a two-aromatic-ring system with at least one heteroatom, requires an
OH or SH within three bonds of the matched ring-system atom, and rejects a primary
amine within two bonds.

By default, nearby-atom matches cannot overlap the main matched atoms. Use
`can_overlap_matched_atoms` if overlap is intended. Use `matched_atom` to limit
distance checks to particular matched atoms; otherwise distances are checked
against any matched atom.

## Atom Queries Without SMARTS

Most proto queries use SMARTS, but `query_atom` can describe atoms and bonds
structurally with proto fields. This is verbose, but useful for generated queries
or where explicit field constraints are clearer than SMARTS.

```textproto
query {
  query_atom {
    id: 0
    atom_properties {
      atomic_number: 7
      ncon: 2
      hcount: 1
    }
    single_bond: 1
  }
  query_atom {
    id: 1
    atom_properties {
      atomic_number: 6
    }
  }
}
```

Each `SubstructureAtom` must have a unique `id`. Bonds can be specified with
fields such as `single_bond`, `double_bond`, `aromatic_bond`, `bond`, or with
`bond_smarts`.

Most hand-written queries should use SMARTS unless the proto form makes the
intent substantially clearer.

## Match-Count And Embedding Controls

Common controls within a `query` block include:

| Field | Effect |
| ----- | ------ |
| `hits_needed`, `min_hits_needed`, `max_hits_needed` | Require a number of embeddings. |
| `max_matches_to_find` | Stop searching after this many embeddings. |
| `one_embedding_per_start_atom` | Keep at most one embedding for each starting atom. |
| `unique_embeddings_only` | Discard duplicate embeddings. |
| `embeddings_do_not_overlap` | Discard embeddings that overlap earlier embeddings. |
| `save_matched_atoms` | Store matched atoms in the result. Defaults to true. |
| `perceive_symmetric_equivalents` | Include symmetric equivalents. Defaults to true. |
| `compress_embeddings` | Merge all embeddings into one embedding. |
| `sort_by_preference_value` | Sort rather than discard low-preference matches. |

Be careful with options that discard embeddings based on discovery order,
especially `embeddings_do_not_overlap` and distance-between-hit constraints.
They can introduce atom-ordering dependence.

## Distances, Paths, And Regions

### No Matched Atoms Between

The SMARTS `...{}` construct can be written as `no_matched_atoms_between` when
multiple such constraints would make the SMARTS unreadable.

```textproto
query {
  smarts: "[OH]-C(=O).c-C(=O)-[OH]"
  no_matched_atoms_between {
    a1: 1
    a2: 3
    qualifier: "{>6}"
  }
}
```

### Link Atoms

`link_atoms` constrains shortest-path distance between two matched atoms. Unlike
`no_matched_atoms_between`, link atoms can pass through matched atoms.

```textproto
query {
  smarts: "[NH>0].O=c"
  link_atoms {
    a1: 0
    a2: 1
    max_distance: 8
  }
}
```

### Separated Atoms

`separated_atoms` also constrains bond separation between matched atoms and can
constrain the number of rotatable bonds on the shortest path.

```textproto
query {
  smarts: "N.c"
  separated_atoms {
    a1: 0
    a2: 1
    min_bonds_between: 3
    max_bonds_between: 7
    max_rotbond: 4
  }
}
```

### Down The Bond

The SMARTS down-the-bond construct can be represented as `down_the_bond`:

```textproto
query {
  smarts: "[Cx>0]-N"
  down_the_bond {
    a1: 0
    a2: 1
    max_natoms: 3
  }
}
```

There is conceptual overlap between `down_the_bond` and `substituent`. Prefer the
construct that makes the query easiest to review.

### Region

A `region` is defined by two matched atoms. The region consists of the contiguous
atoms between those two atoms and can be constrained by size, ring count, and
branching. If both defining atoms are in a ring, the region usually cannot be
unambiguously defined.

```textproto
query {
  smarts: "[D3R]-!@*...{<5;0[R]}*-!@[D3R]"
  region {
    atom: [1, 2]
    max_natoms: 12
    nrings: 0
  }
}
```

This identifies a linker-like region of at most 12 atoms between two rings, with
no rings in the region.

A more permissive replacement-linker query might be:

```textproto
query {
  smarts: "a-!@*...{3-5}*-!@a"
  unique_embeddings_only: true
  region {
    atom: [1, 2]
    max_nrings: 1
    max_atoms_not_on_shortest_path: 2
  }
}
```

Currently the most important implemented region attributes are atom count, ring
count, and atoms not on the shortest path. Check the schema and tests before
relying on less common fields.

### Inter-Ring Regions

`inter_ring_region` describes connected scaffold atoms between rings. This is
usually a better fit than `region` when the concept is specifically "atoms
between rings".

It can constrain the number of atoms, ring connections, and distances between
ring atoms, and can require or reject SMARTS within the region.

## Writing Maintainable Query Files

Prefer small, named query files over one giant query. Good names matter because
many LillyMol tools report the query name as the reason for a match.

Practical conventions:

- include a `name` field in every top-level query;
- keep one chemically meaningful idea per file;
- use global ids whenever a precondition must be tied to the matched embedding;
- use environments and matched-atom constraints to avoid unreadable recursive SMARTS;
- test each query against a few positive and negative example molecules;
- keep comments in the textproto when a field is present for a non-obvious reason.

## Troubleshooting

If a query matches too much, first check whether a ring, ring-system, required
molecular property, or substituent condition is independent of the SMARTS match.
If it must identify the same atoms, add `set_global_id` and match those atoms with
`/IWgidN`.

If a query does not match, simplify it until the SMARTS alone works, then add one
proto construct at a time. Whole-molecule filters and ring specifiers are
preconditions; environment and matched-atom constraints act after an embedding is
found.

If multiple embeddings appear because of symmetry, consider
`unique_embeddings_only`, `one_embedding_per_start_atom`, or more specific SMARTS.
Be aware that some embedding filters are order dependent.

If the query is slow, try cheap whole-molecule preconditions such as
`min_nrings`, `elements_needed`, or `required_bond`. This is not always a large
speed win, but it can avoid obviously impossible atom matching.

## Schema Reference Map

The most useful schema messages are:

| Message | Purpose |
| ------- | ------- |
| `SubstructureQuery` | Top-level container and logical expression across query components. |
| `SingleSubstructureQuery` | One query component: SMARTS, hit controls, environments, ring specs, etc. |
| `RequiredMolecularProperties` | Whole-molecule prefilters. |
| `SubstructureRingSpecification` | One ring satisfying constraints. |
| `SubstructureRingSystemSpecification` | One ring system satisfying constraints. |
| `SubstructureRingBase` | Shared ring/ring-system base constraints such as environment and substituent. |
| `Substituent` | Constraints on a substituent group. |
| `SubstructureEnvironment` | Required or rejected context around matched atoms. |
| `MatchedAtomMatch` | Rooted SMARTS constraints on a matched atom. |
| `NearbyAtoms` | Motifs within a bond range of matched atoms. |
| `SeparatedAtoms` | Bond and rotatable-bond separation between matched atoms. |
| `Region` | Connected region defined by matched atoms. |
| `InterRingAtoms` | Connected scaffold region between rings. |
| `SubstructureAtom` and `SubstructureAtomSpecifier` | Proto-native atom and bond query construction. |

The schema contains many fields not described here. That does not mean they are
unsupported, only that they are less commonly used or more specialized. For exact
field names, numeric types, defaults, and implementation comments, inspect
[substructure.proto](/src/Molecule_Lib/substructure.proto).

## Conclusion

Protocol buffer queries are the right LillyMol format when a chemical idea is
larger than a local atom/bond SMARTS pattern. They are especially useful for
ring-system logic, whole-molecule preconditions, environment exclusions,
substituent constraints, and readable decomposition of complex matching logic.
