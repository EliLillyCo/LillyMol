#!/usr/bin/env ruby
require 'open3'
exe, indir = ARGV
input = File.join(indir, 'input.smi')
def run(exe, input, args, success: true)
  out, err, status = Open3.capture3(exe, *args, input)
  abort "#{args.inspect}: #{err}" unless status.success? == success
  out
end
base = ['-s', 'C-Cl', '-k', '2', '-m', '1', '-M', '20']
plain = run(exe, input, base)
abort 'Default behaviour added caps' if plain.include?('*')
capped = run(exe, input, [*base, '-B', 'chiraldummy'])
expected = {'enantiomer1' => 'Br[C@@H](F)*', 'enantiomer2' => 'Br[C@H](F)*'}
expected.each do |name, smiles|
  abort "Lost orientation for #{name}: #{capped}" unless capped.include?("#{smiles} #{name} FRAGID=")
end
# Multiple potential cuts at a centre must never produce a second cap.
multiple = run(exe, input, ['-B', 'chiraldummy', '-s', 'C-[Cl,F,Br]', '-k', '3', '-m', '1', '-M', '20'])
multiple.lines.grep(/FRAGID=/).each do |line|
  abort "Multiple caps: #{line}" if line.split.first.count('*') > 1
end
recap = run(exe, input, [*base, '-B', 'chiraldummy', '-B', 'recap'])
expected.each_value { |smiles| abort 'RECAP lost cap orientation' unless recap.include?(smiles) }
text = run(exe, input, [*base, '-B', 'chiraldummy', '-C', 'def'])
expected.each do |name, smiles|
  abort 'Text complement lacks cap' unless text.include?("Cl #{name} #{smiles} COMP #{name}")
end
proto = run(exe, input, [*base, '-B', 'chiraldummy', '-C', 'def', '-B', 'proto'])
expected.each_value { |smiles| abort 'Proto complement lacks cap' unless proto.include?(%Q(comp: "#{smiles}")) }
limited = run(exe, input, ['-B', 'chiraldummy', '-s', 'C-Cl', '-k', '2', '-m', '1', '-M', '3'])
abort 'Caps not counted in maximum fragment size' if limited.lines.grep(/FRAGID=/).any? { |line| line.include?('*') }
# A cut elsewhere in a capped fragment can still recurse, and null provenance
# must be handled safely during parent-index mapping.
File.write('chain.smi', "Cl[C@H](F)CCOC chain\n")
run(exe, 'chain.smi', ['-B', 'chiraldummy', '-s', 'C-Cl', '-s', 'C-O', '-k', '3', '-m', '1', '-M', '20'])
%w[BRB BRF BRS dearom].each do |option|
  run(exe, input, [*base, '-B', 'chiraldummy', '-B', option], success: false)
end
run(exe, input, [*base, '-B', 'chiraldummy', '-I', 'enva'], success: false)
run(exe, input, [*base, '-B', 'chiraldummy', '-C', 'fragatype', '-P', 'UST:Z'], success: false)
puts 'chiral dummy checks passed'
