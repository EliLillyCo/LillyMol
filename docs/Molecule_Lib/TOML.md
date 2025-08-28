# TOML
LillyMol also supports TOML and JSON as well as textproto formats.

This is done via the C++ library [toml++](https://marzer.github.io/tomlplusplus/).
tomlplusplus can read TOML and write as JSON. JSON can then be parsed by the
`text_format` utilities in Protocol Buffers. 

When using substructure queries, if the file name specified is prefixed with 'PROTO:'
then the following logic applies.

. If the file name ends in ".json", it is parsed from JSON to Proto.
. If the file name ends in ".toml", it is parsed as TOML, converted to JSON and then to Proto
. Otherwise the file is interpreted as Protocol Buffer textproto form.

Obviously these extra translation steps are expensive, but for the most common
cases of specifying a query or reaction file, it simply does not matter.
If on the other hand you are reading a large number of files, or parsing a
large number of dataitems, then clearly this overhead could become significant.

As an example here is a query from the Lilly Medchem Rules expressed in TOML form.
```
# RuleClass 
# SubClass ester

[query]
  name = "ester"
  numeric_value = 35.0
  smarts = "[OD1]=[CD3R0T2]-[OD2]"
  [query.environment_no_match]
    smarts  = "c"
    [query.environment_no_match.attachment]
      substructure_bond = "- 2"
```
json {
    "query" : {
        "environment_no_match" : {
            "attachment" : {
                "substructure_bond" : "- 2"
            },
            "smarts" : "c"
        },
        "name" : "ester",
        "numeric_value" : 35.0,
        "smarts" : "[OD1]=[CD3R0T2]-[OD2]"
    }
}
```
And as textproto it would be
```
# RuleClass
# SubClass ester

query {
  name: "ester"
  numeric_value: 35
  smarts: "[OD1]=[CD3R0T2]-[OD2]"
  environment_no_match {
    smarts: "c"
    attachment {
      substructure_bond: "- 2"
    }
  }
}
```
A simple acid+amine reaction as TOML might look like
```
[scaffold]
  id = 0
  smarts = "[OH]-C=O"
  remove_atom = [0]
[sidechain]
  id = 1
  smarts = "[ND1H2]-[CX4T1]"
  [sidechain.join]
  a1 = 1
  a2 = 0
  btype = "SS_SINGLE_BOND"
```
As JSON
```
json {
    "scaffold" : {
        "id" : 0,
        "remove_atom" : [
            0
        ],
        "smarts" : "[OH]-C=O"
    },
    "sidechain" : {
        "id" : 1,
        "join" : {
            "a1" : 1,
            "a2" : 0,
            "btype" : "SS_SINGLE_BOND"
        },
        "smarts" : "[ND1H2]-[CX4T1]"
    }
}
```
and as textproto
```
scaffold {
  id: 0
  smarts: "[OH]-C=O"
  remove_atom: [0]
}
sidechain {
  id: 1
  smarts: "[ND1H2]-[CX4T1]"
  join {
    a1: 1
    a2: 0
    btype: SS_SINGLE_BOND
  }
}
```
