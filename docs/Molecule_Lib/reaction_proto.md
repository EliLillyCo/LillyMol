# Reactions as Protos

LillyMol uses Protocol Buffers to describe reactions, specifically
text format Protocol Buffers. Some of the rationale is described
in [why protos](why_protos.md)

## Reactions
Tools like [trxn](/docs/Molecule_Tools/trxn.md) can use text protos to
describe a reaction. As described in that document, LillyMol takes
a very different approach to reactions than
[smirks](https://www.daylight.com/dayhtml/doc/theory/theory.smirks.html).
All actions in a reaction are explicit, rather than implied.

The `trxn` documentation contains many examples of reaction descriptions
in proto form. That document was designed to provide examples of many of
the functionalities avaialble. This document will focus more on
the messages available within the proto definition. Suggest reading
the `trxn` documentation first.

## Query Files
One of the advantages of LillyMol is that complex queries used for a
substructure search can also be used in a reaction unchanged. When using
`smirks` it is common for a smarts to be used for filtering a list
of reagents, and then that smarts must be copied and modified for
use in a `smirks` reaction specification. This imposes a maintenance
burden and risk of divergence between the two.

An example of re-using a query file in a reaction
```
scaffold {
  query_file: "/path/to/file.qry"
}
```
works because that file is also a textproto. That file works with
`tsubstructure`
```
tsubstructure -q PROTO:/path/to/file.qry ...
```
which means the exact same conditions used for filtering a set of
reagents can also be used in the reaction, ensuring concordance. Experience
suggests that the hardest part of enumeration is usually reagent
filtering, and complex queries can be needed.

In fact the preferred strategy may be to use the `-j` option of `tsubstructure` with 
one or more complex queries, to place isotopic labels on the atoms of interest,
and to then have a very simple reaction that just joins isotopes, and perhaps
removes some atoms. If enumerating a large library, that may also be much more
efficient since there will be less work to be done at enumeration time.

It is important to note that directives that refer to other files
are sensitive to `${foo}` type shell variable expansion - if `foo` is set.
```
scaffold {
  query_file: "${QUERIES}/scaffold.qry"
  ...
}
```
which can aid when copying files around. In fact most of LillyMol
uses common file opening tools, so this construct will work in many
different LillyMol contexts.

## Match Conditions
`trxn` contains a variety of command line options that control how
reactions are performed. There are three messages in the reaction
proto definition that correspond to these options. There is a
`MatchConditions` message that applies to both the scaffold and
the sidechains. As always, the functioning of that message is
not described here, please refer to the proto definition.

The messages `ScaffoldMatchConditions` and `SidechainMatchConditions`
both have a `MatchConditions` message - there is no inheritance in
protos.

## MatchedAtomInComponent
This is one of the most important and complex messages in a reaction proto.
The concept is explained in detail in [trxn](/docs/Molecule_Tools/trxn.md).

## StereoChemistry
Reactions can be used to set chirality via the `StereoCenter` message. As can
be seen from the message, a stereo centre can involve any matched atom in
the scaffold and/or one or more sidechains.

Unfortunately in the transition to proto form, stereochemistry became
more complicated than previously. While there is great flexibility
it is needlessly complex in the simple cases. For example
```
scaffold {
  smarts: "[Nx0]-[Cx0D3H](-[CD1])-[CD2]"
}
reaction_stereo_center {
  center {
    atom {
      component: 0
      atom: 1
    }
  }
  top_front {
    atom {
      component: 0
      atom: 0
    }
  }
  top_back {
    atom {
      component: 0
      atom: 2
    }
  }
  left_down {
    atom {
      component_and_atom: "0.3"
    }
  }
  right_down {
    implicit_hydrogen: true
  }
}
```
involved just atoms in the same reaction component. We should have the ability
to do something simple like
```
scaffold {
  smarts: ...
  stereo_center {
    centre: 1
    top_front: 0
    top_back: 2
    left_down: 3
    right_down: IMPLICIT_HYDROGEN
  }
}
```
with all atom numbers being local to the component. Or perhaps a variant of the
above that only takes the string description of a matched atom `component:atom`.

While stereochemistry does work, it is difficult to work with.
