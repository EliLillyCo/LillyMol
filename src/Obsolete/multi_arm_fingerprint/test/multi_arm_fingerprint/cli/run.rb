#!/usr/bin/env ruby
require 'open3'

exe, indir = ARGV
input = File.join(indir, 'input.smi')
def run(exe, args, stdin = '')
  stdout, stderr, status = Open3.capture3(exe, *args, stdin_data: stdin)
  abort "#{args.inspect}: #{stderr}" unless status.success?
  [stdout, stderr]
end
def fingerprints(text)
  text.lines.grep(/^NCMARM</)
end
base, = run(exe, [input])
abort 'Wrong record count' unless base.lines.count { |line| line == "|\n" } == 3
abort 'Missing empty fingerprint' unless fingerprints(base)[1] == "NCMARM<>\n"
stdin, = run(exe, ['-'], File.read(input))
abort 'stdin differs' unless stdin == base
selected, = run(exe, ['-z', 'i', '-s', '[CD3]', input])
abort 'Selecting all roots differs' unless selected == base
duplicated, = run(exe, ['-z', 'i', '-s', '[CD3]', '-s', '[CD3]', input])
abort 'Duplicate roots counted twice' unless duplicated == base
missing, = run(exe, ['-z', 'i', '-s', '[#7]', input])
abort 'No-match selection should retain empty records' unless
  fingerprints(missing) == ["NCMARM<>\n"] * 3
# Both endpoints of this symmetric query are roots. Unique-embedding or symmetry
# suppression would select just one endpoint and change the counted fingerprint.
query_selected, = run(exe, ['-z', 'i', '-s', '[CD3]-[CD3]', input])
abort 'Symmetric query lost a root' unless fingerprints(query_selected)[2] == fingerprints(base)[2]
file_selected, = run(exe, ['-z', 'i', '-q', "PROTO:#{File.join(indir, 'root.textproto')}", input])
abort '-q and -s differ' unless file_selected == query_selected
filtered, = run(exe, ['-f', '-'], base.gsub(/^NCMARM<.*\n/, ''))
abort 'Filter fingerprint differs' unless fingerprints(filtered) == fingerprints(base)
abort 'Filter changed existing records' unless
  filtered.gsub(/^NCMARM<.*\n/, '') == base.gsub(/^NCMARM<.*\n/, '')
all_roots, = run(exe, ['-z', 'i', '-s', '[CD3]-C', input])
abort 'Multiple embeddings counted a root twice' unless all_roots == base
custom, = run(exe, ['-J', 'NCCUSTOM', input])
abort 'Custom tag ignored' unless custom.lines.grep(/^NCCUSTOM</).size == 3
_, report = run(exe, ['-v', '-R', '3', input])
abort 'Missing statistics' unless report.include?('centres') && report.include?('combinations')
_, bit_report = run(exe, ['-v', '-R', '3', '-'], "CC(C)C branched\nCCCC linear\n")
abort 'Incorrect bit-count summary' unless
  bit_report.include?('Distinct bits per molecule: min 0 max 1 average 0.5')
_, empty_report = run(exe, ['-v', '-'], '')
abort 'Empty input should not report undefined bit statistics' if
  empty_report.include?('Distinct bits per molecule:')
stdout, stderr, status = Open3.capture3(exe, '-R', '2', input)
abort 'Invalid bounds accepted' if status.success?
_, error, status = Open3.capture3(exe, '-s', '[#7]', input)
abort 'Query nonmatch should fail by default' if status.success?
abort 'Missing no-match diagnostic' unless error.include?('no query matches')
_, unmatched_report = run(exe, ['-v', '-z', 'i', '-s', '[#7]', input])
abort 'Wrong no-match count' unless unmatched_report.include?('3 molecules did not match any queries')
# A match at an ineligible root is not a query nonmatch.
linear, linear_report = run(exe, ['-v', '-s', 'C', '-'], "CCCC linear\n")
abort 'Ineligible matched roots should produce an empty fingerprint' unless
  fingerprints(linear) == ["NCMARM<>\n"] &&
  linear_report.include?('0 molecules did not match any queries')
# Only the first matched atom selects a root: terminal carbons are ineligible,
# even though the second query atom matches the branching carbon.
first_only, = run(exe, ['-s', '[CD1]-[CD3]', '-'], "CC(C)C branched\n")
abort 'Selected atoms beyond the first query atom' unless fingerprints(first_only) == ["NCMARM<>\n"]
no_fp = base.gsub(/^NCMARM<.*\n/, '')
no_match_filter, filter_report = run(exe, ['-v', '-f', '-z', 'i', '-s', '[#7]', '-'], no_fp)
abort 'Filter mode lost nonmatching records' unless no_match_filter.gsub(/^NCMARM<.*\n/, '') == no_fp
abort 'Filter mode should emit empty fingerprints for nonmatches' unless
  fingerprints(no_match_filter) == ["NCMARM<>\n"] * 3 &&
  filter_report.include?('3 molecules did not match any queries')
_, _, status = Open3.capture3(exe, '-z', 'unknown', input)
abort 'Unknown -z qualifier accepted' if status.success?
sync, sync_report = run(exe, ['-S', '-L', '3', '-R', '12', '-v', '-'], "CCC(CC)CC balanced\n")
abort 'Wrong synchronised combination count' unless sync_report.include?('2 combinations')
normal, normal_report = run(exe, ['-L', '3', '-R', '12', '-v', '-'], "CCC(CC)CC balanced\n")
abort 'Original enumeration changed' unless normal_report.include?('8 combinations')
abort 'Synchronised fingerprint should differ from all combinations' if fingerprints(sync) == fingerprints(normal)
puts 'CLI checks passed'
