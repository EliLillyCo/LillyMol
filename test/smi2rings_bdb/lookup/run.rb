#!/usr/bin/env ruby
require 'open3'
exe, indir = ARGV
one = File.join(indir, 'one.smi')
two = File.join(indir, 'two.smi')
query = File.join(indir, 'query.smi')
def run(exe, args, success = true)
  stdout, stderr, status = Open3.capture3(exe, *args)
  abort "#{args.inspect}: #{stderr}" unless status.success? == success
  stdout
end
def counts(text)
  text.lines.to_h do |line|
    fields = line.split
    [fields[1], fields.find { |field| field.start_with?('URS:') }.delete_prefix('URS:').to_i]
  end
end
labels = ['-j', 'iso', '-j', 'ring']
run(exe, ['-d', 'STORE', '-d', 'one.bdb', '-n', *labels, one])
run(exe, ['-d', 'STORE', '-d', 'two.bdb', '-n', *labels, two])
lookup = ['-d', 'LOOKUP', '-n', *labels]
single = run(exe, [*lookup, '-d', 'one.bdb', query])
expected = {'query' => 1, 'fallback' => 2, 'absent' => 0}
abort "Wrong single-database counts: #{single}" unless counts(single) == expected
multiple = run(exe, [*lookup, '-d', 'one.bdb', '-d', 'two.bdb', query])
abort "More-specific match lost to fallback: #{multiple}" unless counts(multiple) == expected
reversed = run(exe, [*lookup, '-d', 'two.bdb', '-d', 'one.bdb', query])
abort 'Database ordering changed counts' unless counts(reversed) == expected
# The highest count wins across databases; counts are not summed.
abort 'Counts were summed' unless counts(multiple)['fallback'] == 2
proto = run(exe, [*lookup, '-d', 'one.bdb', '-d', 'two.bdb', '-Y', 'proto', query])
abort 'Proto output lost attachment-specific match' unless proto.lines.first.include?('urs: 1')
# In-memory aggregation must produce the same lookups as direct storage.
run(exe, ['-d', 'STORE', '-d', 'cached.bdb', '-n', '-Y', 'hash', *labels, one])
cached = run(exe, [*lookup, '-d', 'cached.bdb', query])
abort 'Cached build differs from direct build' unless cached == single
# Existing lookup filters and bare-ring-only mode remain available.
filtered = run(exe, [*lookup, '-d', 'one.bdb', '-d', 'two.bdb', '-f', '2', query])
abort 'URS filter changed' unless counts(filtered) == {'fallback' => 2}
bare = run(exe, ['-d', 'LOOKUP', '-d', 'one.bdb', '-d', 'two.bdb', '-n', '-j', 'ring', query])
abort 'Bare-ring counts changed' unless counts(bare) == {'query' => 2, 'fallback' => 2, 'absent' => 0}
puts 'ring database lookup checks passed'
