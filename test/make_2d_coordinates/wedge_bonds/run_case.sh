#!/bin/bash
# Args from run_all_test.rb: executable, indir, outdir, then TestCase args.
#
# Chirality only reaches other software through wedge bonds. LillyMol writes it
# in the MDL atom parity field as well, but other toolkits ignore parity on a 2D
# record, so the bond block is what matters here.
#
# Two things are checked for every wedge, both read straight out of the molfile
# so that no other tool's interpretation gets in the way:
#
#   - which bond carries it, printed in the order the bond block stores it. The
#     first atom of a wedge is its narrow end and has to be the stereocentre;
#     written the other way round the wedge is silently discarded by readers
#     that take the first atom as the centre, RDKit among them. The parity field
#     tells us which atoms LillyMol considers chiral, so the two halves of the
#     record are cross checked against each other and 'NOT-A-CENTRE' would show
#     up here.
#   - which way it points, so that the two enantiomers below have to come out
#     as each other's opposite.
#
# The coordinates themselves are not compared - see the roundtrip case for why.
# Note that the sense of a wedge is only meaningful together with the layout it
# is drawn on: were some platform to lay these molecules out mirrored, every
# 'up' and 'down' here would swap while the depictions stayed correct. That has
# not been seen, and the sense is what a reader ultimately consumes, so it is
# compared directly rather than relative to the first record.
exe="$1"; datadir="$4"

# Reports the wedges of each record: the bond, in the order the bond block
# stores it, and which way it points.
report_wedges='
function report() {
  if (name == "") {
    return
  }
  if (nw == 0) {
    printf "%s: no wedges\n", name
  } else {
    line = name ":"
    for (i = 1; i <= nw; i++) {
      # The narrow end of a wedge is its first atom, and it must be an atom
      # LillyMol marked chiral in the parity field.
      what = (parity[wa1[i]] != 0) ? "" : " NOT-A-CENTRE"
      line = line sprintf(" %d-%d %s%s", wa1[i], wa2[i], sense[wf[i]], what)
    }
    print line
  }
  name = ""
}

BEGIN {
  sense[1] = "up"
  sense[6] = "down"
  sense[4] = "either"
}

/^\$\$\$\$/ { report(); n = 0; next }

{ n++ }

n == 1 { name = $0; nw = 0; delete parity; next }

n == 4 {
  natoms = substr($0, 1, 3) + 0
  nbonds = substr($0, 4, 3) + 0
  next
}

# Atom block. Columns 40-42 are the atom stereo parity.
n >= 5 && n <= 4 + natoms {
  parity[n - 4] = substr($0, 40, 3) + 0
  next
}

# Bond block. Columns 10-12 are the wedge/hash flag.
n >= 5 + natoms && n <= 4 + natoms + nbonds {
  flag = substr($0, 10, 3) + 0
  if (flag != 0) {
    nw++
    wa1[nw] = substr($0, 1, 3) + 0
    wa2[nw] = substr($0, 4, 3) + 0
    wf[nw] = flag
  }
  next
}
'

"${exe}" -S wedged "${datadir}/stereo.smi" || exit $?

# -w must turn the whole thing off, leaving the parity field as the only record
# of chirality, which is what the tool did before wedges were implemented.
"${exe}" -w -S nowedge "${datadir}/stereo.smi" || exit $?

{
  awk "${report_wedges}" wedged.sdf
  echo 'with -w:'
  awk "${report_wedges}" nowedge.sdf
} > wedges
