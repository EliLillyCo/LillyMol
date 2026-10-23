#!/bin/bash
# Args from run_all_test.rb: executable, indir, outdir, then TestCase args.
#
# A molfile has no field for cis/trans - E/Z is implicit in which side of the
# double bond each substituent is drawn on - so the only way to check it is to
# lay a molecule out and read the configuration back out of the coordinates.
# That is what fileconv -i dctb does, and it is the whole chain end to end.
#
# The coordinates themselves are not compared, for the reason the roundtrip case
# explains. What is compared here is a sign rather than a number: which side of
# the double bond axis a substituent came out on. That is far from zero in any
# usable layout, so it does not vary between platforms the way the coordinates
# do. The direction markers LillyMol writes may be the other pair of the two
# that mean the same thing (C\C=C\C rather than C/C=C/C), which is why the
# golden file holds whatever it writes rather than the input smiles.
#
# -z is the opt out, and is here because it is what documents what the layout
# engine does on its own: it draws every double bond trans, so with -z the cis
# structures come back as their trans isomers. Anything else appearing under -z
# would mean the constraint is not the thing making the difference.
exe="$1"; datadir="$4"
bindir=$(dirname "${exe}")

"${exe}" -S honoured "${datadir}/cis_trans.smi" || exit $?
"${exe}" -z -S ignored "${datadir}/cis_trans.smi" || exit $?

# -i dctb is what asks fileconv to derive cis/trans from the depiction. Without
# it every one of these comes back with no double bond stereochemistry at all.
"${bindir}/fileconv" -i sdf -i dctb -o smi -S from_honoured honoured.sdf || exit $?
"${bindir}/fileconv" -i sdf -i dctb -o smi -S from_ignored ignored.sdf || exit $?
"${bindir}/fileconv" -i sdf -o smi -S no_dctb honoured.sdf || exit $?

{
  echo 'laid out to match the molecule:'
  cat from_honoured.smi
  echo 'with -z, the raw layout:'
  cat from_ignored.smi
  echo 'read back without -i dctb, so no cis/trans is derived at all:'
  cat no_dctb.smi
} > cis_trans
