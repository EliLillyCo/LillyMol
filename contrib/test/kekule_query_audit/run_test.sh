#!/usr/bin/env bash
# Tests contrib/bin/kekule_query_audit.py on queries with known behaviour.
#
#   LILLYMOL_HOME=/path/to/LillyMol ./run_test.sh
#
# tsubstructure must be on the PATH, or found via LILLYMOL_HOME.

here=$(cd $(dirname $0) && pwd)
audit=${here}/../../bin/kekule_query_audit.py
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
fail=0

check() {  # description, condition result (0 = pass)
  if [[ $2 -eq 0 ]]; then echo "ok    $1"; else echo "FAIL  $1"; fail=1; fi
}

mkdir -p $tmp/queries
cat > $tmp/molecules.smi <<'MOLS'
c1ccccc1 benzene
c1ccccc1-c1ccccc1 biphenyl
CCO ethanol
c1ccc2ccccc2c1 naphthalene
c1ccc(cc1)-c1ccccn1 phenylpyridine
MOLS
printf 'c1ccccc1-c1ccccc1 biphenyl\n' > $tmp/biphenyl.smi

# a-a means different things with and without the Kekule forms of aromatic bonds
printf '(0 Query\n  (A C smarts "a-a")\n)\n'   > $tmp/queries/dependent.qry
printf '(0 Query\n  (A C smarts "a-!@a")\n)\n' > $tmp/queries/independent.qry
printf 'name: "dependent"\nquery {\n  smarts: "a-a"\n}\n' > $tmp/queries/dependent_proto.txtproto
printf '(0 Query\n  (A C smarts "C(")\n)\n'   > $tmp/queries/broken.qry
printf 'dependent.qry\n' > $tmp/queries/a_list
# not a substructure query: should be ignored when the directory is searched
printf 'name: "not a query"\n' > $tmp/queries/config.textproto
# a textproto query with the .textproto suffix, which is used when it is named
printf 'name: "named"\nquery {\n  smarts: "a-a"\n}\n' > $tmp/named.textproto

run() { python3 $audit "$@" 2>$tmp/stderr; }

out=$(run --jobs 2 $tmp/molecules.smi $tmp/queries)
row() { echo "$out" | grep -F "$1" | head -1; }

# default mode matches benzene, biphenyl, naphthalene and phenylpyridine, -M nokekule only the last two
[[ "$(row dependent.qry)" == *$'\tCOUNT_DIFF\t4\t2\t-2'* ]]; check "a-a depends on the mode" $?
[[ "$(row dependent_proto.txtproto)" == *$'\tCOUNT_DIFF\t4\t2\t-2'* ]]; check "textproto query depends on the mode" $?
[[ "$(row broken.qry)" == *$'\tERROR\t'* ]]; check "a query that cannot be read is an error" $?
echo "$out" | grep -q "independent.qry"; [[ $? -ne 0 ]]; check "a-!@a is not listed" $?
grep -q "1 depend\|2 depend on the mode" $tmp/stderr; check "summary counts the queries that depend on the mode" $?

out=$(run --all --jobs 2 $tmp/molecules.smi $tmp/queries)
[[ "$(echo "$out" | grep -F independent.qry)" == *$'\tsame\t2\t2\t0'* ]]; check "--all lists a-!@a as the same in both modes" $?

out=$(run --used-in --jobs 2 $tmp/molecules.smi $tmp/queries)
[[ "$(echo "$out" | grep -F '/dependent.qry')" == *a_list* ]]; check "--used-in names the list that has the query" $?

out=$(run --fkekule --all --jobs 2 $tmp/molecules.smi $tmp/queries/dependent.qry)
[[ "$(echo "$out" | head -1)" == *matches_fkekule* ]]; check "--fkekule adds a column" $?

# Same molecule is found in both modes, but different atoms are labelled
out=$(run --labels --jobs 2 $tmp/biphenyl.smi $tmp/queries/dependent.qry)
[[ "$(echo "$out" | grep -F dependent.qry)" == *$'\tLABEL_DIFF\t1\t1\t0'* ]]; check "--labels finds a change in the labelled atoms" $?
out=$(run --jobs 2 $tmp/biphenyl.smi $tmp/queries/dependent.qry)
[[ -z "$(echo "$out" | grep -F dependent.qry)" ]]; check "counts alone do not see that change" $?

run --fail-on-change --jobs 2 $tmp/molecules.smi $tmp/queries >/dev/null; [[ $? -eq 1 ]]; check "--fail-on-change gives status 1" $?
run --fail-on-change --jobs 2 $tmp/molecules.smi $tmp/queries/independent.qry >/dev/null; [[ $? -eq 0 ]]; check "--fail-on-change gives status 0 when nothing depends on the mode" $?

run --head 2 --jobs 2 $tmp/molecules.smi $tmp/queries/dependent.qry >/dev/null
grep -q "2 molecules" $tmp/stderr; check "--head limits the molecules" $?
run --every 2 --jobs 2 $tmp/molecules.smi $tmp/queries/dependent.qry >/dev/null
grep -q "3 molecules" $tmp/stderr; check "--every samples the molecules" $?

out=$(run --all --jobs 2 $tmp/molecules.smi $tmp/queries)
echo "$out" | grep -q "config.textproto"; [[ $? -ne 0 ]]; check "a .textproto file in a directory is not taken as a query" $?

out=$(run --jobs 2 $tmp/molecules.smi $tmp/named.textproto)
[[ "$(echo "$out" | grep -F named.textproto)" == *$'\tCOUNT_DIFF\t4\t2\t-2'* ]]; check "a .textproto file named on the command line is used" $?

exit $fail
