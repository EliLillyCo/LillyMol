#!/usr/bin/env bash
# Checks the recipes in SKILL.md against a few inline molecules with known answers.
# Run it again after upgrading LillyMol to see whether any documented behaviour changed.
#
#   LILLYMOL_HOME=/path/to/LillyMol PYTHON=/path/to/matching/python ./test_examples.sh
#
# Needs tsubstructure on the PATH. The python checks need run_python.sh to work
# (set PYTHON to the interpreter the bindings were built for), and are skipped otherwise.

set -u
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
fail=0

cat > "$tmp/t.smi" <<'EOF'
c1ccccc1-c1ccccc1 biphenyl
c1ccccc1 benzene
CS(=O)(=O)N sulfonamide
CC#N acetonitrile
N#CCc1ccccc1 phenylacetonitrile
CS(=O)(=O)C dimethylsulfone
CCO ethanol
EOF

# usage: check "description" expected_count tsubstructure_args...
check() {
  local desc=$1 want=$2; shift 2
  local got
  got=$(tsubstructure "$@" "$tmp/t.smi" 2>&1 | sed -n 's/.* \([0-9]*\) molecules match.*/\1/p')
  if [[ "$got" == "$want" ]]; then
    printf 'ok    %s\n' "$desc"
  else
    printf 'FAIL  %s: expected %s, got %s\n' "$desc" "$want" "${got:-nothing}"; fail=1
  fi
}

check "count nitriles"                          2 -s 'C#N'
check "count sulfonyl"                          2 -s '[SX4](=O)=O'
check "any of two queries"                      4 -s 'C#N' -s '[SX4](=O)=O'
check "all queries must match (-M mmaq)"       1 -M mmaq -s 'C#N' -s 'c'
check "-b keeps the any-match result"         4 -b -s 'C#N' -s 'c'
check "-B depends on query order (1)"          2 -B -s 'C#N' -s 'c'
check "-B depends on query order (2)"          3 -B -s 'c' -s 'C#N'
printf 'C#N nitrile\n[SX4](=O)=O sulfonyl\n' > "$tmp/q.smt"
check "smarts file (-q S:)"                     4 -q S:"$tmp/q.smt"
check "a-a also matches Kekule bonds in rings"  3 -s 'a-a'
check "a-!@a is a biaryl"                       1 -s 'a-!@a'
check "-M nokekule makes a-a a biaryl"          1 -M nokekule -s 'a-a'

# -m / -n write files named <stem>.smi
tsubstructure -s 'C#N' -m "$tmp/hit" -n "$tmp/miss" "$tmp/t.smi" >/dev/null 2>&1
hits=$(wc -l < "$tmp/hit.smi"); misses=$(wc -l < "$tmp/miss.smi")
if [[ "$hits" == 2 && "$misses" == 5 ]]; then echo "ok    -m and -n"; else echo "FAIL  -m/-n: $hits hits, $misses misses"; fail=1; fi

# embeddings versus unique matches on a symmetric query
e=$(tsubstructure -s '[SX4](=O)=O s' -m - -m QDT "$tmp/t.smi" 2>/dev/null | grep dimethylsulfone | grep -o "([0-9]* matches" | head -1)
u=$(tsubstructure -u -s '[SX4](=O)=O s' -m - -m QDT "$tmp/t.smi" 2>/dev/null | grep dimethylsulfone | grep -o "([0-9]* matches" | head -1)
if [[ "$e" == "(2 matches" && "$u" == "(1 matches" ]]; then echo "ok    embeddings (2) versus -u (1)"; else echo "FAIL  -u: '$e' '$u'"; fail=1; fi

# stdin
n=$(cat "$tmp/t.smi" | tsubstructure -s 'C#N' - 2>&1 | sed -n 's/.* \([0-9]*\) molecules match.*/\1/p')
if [[ "$n" == 2 ]]; then echo "ok    stdin"; else echo "FAIL  stdin: $n"; fail=1; fi

# python bindings
if [[ -n "${LILLYMOL_HOME:-}" && -x "${LILLYMOL_HOME}/run_python.sh" ]]; then
  cat > "$tmp/t.py" <<EOF
from lillymol import *
from lillymol_tsubstructure import *
ts = TSubstructure()
ts.add_query_from_smarts('a-!@a biaryl')
ts.add_query_from_smarts('[SX4](=O)=O sulfonyl')
mols = []
with MolReaderContext('$tmp/t.smi') as reader:
  for m in reader:
    mols.append(m)
counts = ts.num_matches(mols)
ok = (len(mols) == 7 and sum(ts.substructure_search(mols)) == 3
      and sum(1 for c in counts if c[0] > 0) == 1 and sum(1 for c in counts if c[1] > 0) == 2)
t2 = TSubstructure(); t2.add_query_from_smarts('a-a')
ok = ok and sum(t2.substructure_search(mols)) == 3
if hasattr(__import__('lillymol'), 'set_aromatic_bonds_lose_kekule_identity'):
  set_aromatic_bonds_lose_kekule_identity(1)
  ok = ok and sum(t2.substructure_search(mols)) == 1
  set_aromatic_bonds_lose_kekule_identity(0)
print('PYOK' if ok else 'PYFAIL')
EOF
  out=$("${LILLYMOL_HOME}/run_python.sh" "$tmp/t.py" 2>&1)
  if echo "$out" | grep -q PYOK; then echo "ok    python TSubstructure";
  elif echo "$out" | grep -q "ABI mismatch"; then echo "SKIP  python (ABI mismatch: set PYTHON to the interpreter the bindings were built for)"
  else echo "FAIL  python: $out" | tail -3; fail=1; fi
else
  echo "SKIP  python (LILLYMOL_HOME or run_python.sh not found)"
fi

exit $fail
