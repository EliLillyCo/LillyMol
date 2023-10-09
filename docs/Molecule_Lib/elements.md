# Elements

LillyMol contains a periodic table data structure, and each
element is an Element object. Each Atom has a pointer to
the Element that defines that Atom. Every atom must have an
Element.

LillyMol is unusual in that the universe of elements is not
restricted to the existing periodic table. Any two letter
combination can be an element, being added to the periodic
table as needed. Once that functionality was in place, we
realised that there was nothing special about two character
elements, so the concept was extended to any string.

HOWTO:

Most LillyMol tools have a `-E` option which enables customisation
of the elements. The usage message is
```
  -E autocreate  automatically create new elements when encountered
  -E PTABLE=file use 'file' for element data
  -E O:Xx        element Xx is considered Organic
  -E anylength   elements can be of arbitrary length
  -E nsqb=xx     element <xx> needs square brackets
  -E <symbol>    automatically create \<symbol\>
```

### -E autocreate
Automatically create new elements when encountered. This applies to
two letter combinations only.
All combinations of the form `[A-Z][a-z]` are recognised as elements.
This also enables such things as `R1` as a valid element. But non
alphanumeric characters are not permitted.

### -e PTABLE=file
Use 'file' for element data. This is retained for historical reasons
only. No current use cases are anticipated.

### -E O:Xx
Element Xx is considered Organic. A common use might be `-E O:Se`.

### -E anylength
Elements can be of arbitrary length. So things like `Ala` and `Phe`
become valid elements. The restriction to alphanumeric characters
is also lifted. Using this implies `-E autocreate` so you do not
need to specify both options.

If you need to perform substructure searches against these elements,
you will need to use the #{symbol} smarts extension. See the file
`Molecule_Tools/introduction.cc` for examples of how to do substructure
searches over peptides.

### -E nsqb=xx
Element <xx> needs square brackets. Use case is unclear, almost all
elements need square brackets. Mostly this gets used with Hydrogen.
We have encountered explicit Hydrogens written as either `H` or `[H]`.

### -E <symbol>
Automatically create \<symbol\>. Generally it is easier to just use
`autocreate` and allow LillyMol to create whatever otherwise unrecognised
elements it encounters. But for more control, the set of acceptable
extra elements can be specified via the `-E` option.

