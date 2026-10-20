#!/usr/bin/env ruby
require 'open3'
exe, indir = ARGV
input = File.join(indir, 'input.gfp')
fp = ['-F', 'FPDSC,w=0.2', '-F', 'NCSELW,nc,w=0.8']
ids = %w[PBCHM599243 PBCHM773628 duplicate]
def run(exe, args, success = true)
  stdout, stderr, status = Open3.capture3(exe, *args)
  abort "#{args.inspect}: #{stderr}" unless status.success? == success
  [stdout, stderr]
end
def records(text)
  text.split("|\n").reject(&:empty?).map { |record| record + "|\n" }
end
def target(record)
  record[/^PCN<([^>]+)>/, 1]
end
base, = run(exe, [*fp, '-n', '2', input])
all = records(base)
abort 'Wrong baseline targets' unless all.map { |r| target(r) } == ids
single, = run(exe, [*fp, '-n', '2', '-O', ids[0], input])
abort 'Single query differs from full run' unless records(single) == [all[0]]
abort 'Nonselected fingerprints must remain candidate neighbours' unless single.include?('PCN<duplicate>')
selections = [ ['-O', ids[2], '-O', ids[0]], ['-O', "#{ids[2]},#{ids[0]}"],
               ['-O', "#{ids[2]},#{ids[0]}", '-O', ids[0], '-O', ids[2]] ]
selections.each do |selection|
  selected, report = run(exe, [*fp, '-n', '2', '-v', *selection, input])
  abort 'Selection changed output/order or duplicated a query' unless records(selected) == [all[0], all[2]]
  abort 'Missing selected-query statistics' unless report.include?('2 selected queries from 3 fingerprints')
end
# Check the same paths for thresholds, omission, index output, and Tversky.
[ ['-T', '0.1'], ['-T', '0.1', '-z'], ['-o', '-x', '-n', '2'],
  ['-n', '2', '-A', 'AVDIST', '-N', 'NCNT', '-j', '3'],
  ['-n', '2', '-V', 'a=0.3', '-V', 'b=1.7'] ].each do |options|
  full, = run(exe, [*fp, *options, input])
  selected, = run(exe, [*fp, *options, '-O', "#{ids[2]},#{ids[0]}", input])
  expected = records(full).select { |record| [ids[0], ids[2]].include?(target(record)) }
  abort "Selected output differs: #{options.inspect}" unless records(selected) == expected
end
# Compare complete framed proto records, including CRCs, without requiring Python
# protobuf bindings or a separate decoder in the regression environment.
def framed_records(fname)
  data = File.binread(fname)
  offset = 0
  result = []
  while offset < data.bytesize
    length = data.byteslice(offset, 8).unpack1('Q<')
    size = length + 16
    record = data.byteslice(offset, size)
    abort 'Incomplete TFDataRecord' unless record && record.bytesize == size
    result << record
    offset += size
  end
  result
end
run(exe, [*fp, '-n', '2', '-S', 'full.tfrecord', input])
run(exe, [*fp, '-n', '2', '-S', 'selected.tfrecord', '-O', "#{ids[2]},#{ids[0]}", input])
full_protos = framed_records('full.tfrecord')
abort 'Wrong full proto count' unless full_protos.size == 3
abort 'Selected protos differ from full run' unless framed_records('selected.tfrecord') == [full_protos[0], full_protos[2]]
# All invalid selections must fail before opening an existing -S file.
File.binwrite('preserve.tfrecord', 'sentinel')
[ ['-O', 'missing'], ['-O', "#{ids[0]},missing"], ['-O', ''],
  ['-O', "#{ids[0]},"], ['-O', ",#{ids[0]}"], ['-O', "#{ids[0]},,#{ids[2]}"],
  ['-O', ids[0], '-p'] ].each do |selection|
  run(exe, [*fp, '-S', 'preserve.tfrecord', *selection, input], false)
  abort 'Output was opened before validating selection' unless File.binread('preserve.tfrecord') == 'sentinel'
end
File.write('duplicate_ids.gfp', File.read(input) + records(File.read(input))[0])
_, error = run(exe, [*fp, '-O', ids[0], 'duplicate_ids.gfp'], false)
abort 'Missing ambiguity diagnostic' unless error.include?('matches multiple fingerprints')
run(exe, [*fp, '-O', ids[1], 'duplicate_ids.gfp']) # Unselected duplicate IDs remain allowed.
run(exe, [*fp, '-s', '1', '-O', ids[2], input], false)
puts 'selected-query checks passed'
