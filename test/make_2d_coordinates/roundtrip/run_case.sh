#!/bin/bash
# Args from run_all_test.rb: executable, indir, outdir, then TestCase args.
# The generated coordinates are deliberately NOT compared. They come out of a
# floating point minimisation and differ in the last digit between compilers,
# platforms and optimisation levels, so a golden molfile would be a permanently
# fragile test. What must hold everywhere is that the molfiles are readable and
# describe the molecules that went in, so the structures are round tripped back
# through fileconv and only the smiles are compared.
exe="$1"; datadir="$4"
bindir=$(dirname "${exe}")

"${exe}" -S coords "${datadir}/input.smi" || exit $?
"${bindir}/fileconv" -i sdf -o smi -S roundtrip coords.sdf || exit $?
