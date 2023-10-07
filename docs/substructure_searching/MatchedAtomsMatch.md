# MatchedAtomsMatch

This message allows easy specification of different possibilities
for an atom match.

## Smarts
Smarts for complex atom specifications can be very hard to understand.
When confronted with a complex smarts, I have often found myself using
an editor to tokenize the smarts onto separate lines in order to understand
what is going on.

This is partly why some have referred to smarts as a write only language.

Online smarts visualization tools are sometimes helpful, other times not.

## QueryAtom
The proto representation of a query does accommodate alternate forms
for an atom, but it is verbose. For example a query for an atom that
was either a primary amine, or a secondary amine might be
```
```
