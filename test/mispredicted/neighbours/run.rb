#!/usr/bin/env ruby
require 'open3'
exe, indir = ARGV
nn = File.join(indir, 'neighbours.tfrecord')
exp = File.join(indir, 'experimental.txt')
pred = File.join(indir, 'predictions.txt')
def run(exe, args, success = true)
  stdout, stderr, status = Open3.capture3(exe, *args)
  abort "#{args.inspect}: #{stderr}" unless status.success? == success
  [stdout, stderr]
end
baseline, = run(exe, ['-E', exp, '-d', pred])
with_nn, coverage = run(exe, ['-E', exp, '-d', '-v', '-N', nn, '-Y', 'summary.txt', '-T', 'results', pred])
lines = with_nn.lines.map(&:split)
abort 'Wrong NN header' unless lines.shift == %w[ID EXPT AVDIFF NN_DIST]
abort 'Predictions or ranking changed' unless lines.map { |line| line[0...-1] } == baseline.lines.drop(1).map(&:split)
expected = {'a' => '0.2', 'b' => '.', 'c' => '.', 'd' => '0', 'e' => '.'}
lines.each { |row| abort "Wrong distance #{row.inspect}" unless row[-1] == expected[row[0]] }
abort 'Wrong coverage counts' unless coverage.include?('2 prediction identifiers with distances, 1 absent IDs, 2 empty neighbour lists')
summary = File.readlines('summary.txt').map(&:split)
abort 'Wrong summary header' unless summary.shift[-1] == 'NN_DIST'
abort 'Malformed summary rows' unless summary.size == 5 && summary.all? { |row| row.size == 10 && row[-1] == expected[row[0]] }
protos = File.readlines('results.textproto')
abort 'Wrong proto row count' unless protos.size == 5
protos.each do |line|
  id = line[/id: "([^"]+)"/, 1]
  distance = line[/nearest_neighbour_distance: ([\d.eE+-]+)/, 1]
  if expected[id] == '.'
    abort 'Missing distance should be absent from proto' if distance
  else
    abort 'Wrong proto distance' unless distance && (distance.to_f - expected[id].to_f).abs < 1.0e-6
  end
end
text, = run(exe, ['-E', exp, '-N', nn, pred])
abort 'Missing labelled distance in normal output' unless text.lines.all? { |line| line.include?(' NN_DIST ') }
cls = File.join(indir, 'classes.txt')
cpred = File.join(indir, 'predicted_classes.txt')
class_base, = run(exe, ['-C', '-E', cls, '-d', cpred])
class_nn, = run(exe, ['-C', '-E', cls, '-d', '-N', nn, cpred])
abort 'Classification header missing NN_DIST' unless class_nn.lines.first.split[-1] == 'NN_DIST'
abort 'Classification ranking changed' unless class_nn.lines.drop(1).map { |line| line.split[0...-1] } == class_base.lines.drop(1).map(&:split)
class_nn.lines.drop(1).each do |line|
  row = line.split
  abort 'Wrong classification distance' unless row[-1] == expected[row[0]]
end
run(exe, ['-E', exp, '-N', 'missing.tfrecord', pred], false)
run(exe, ['-E', exp, '-N', nn, '-N', nn, pred], false)
before = File.binread(nn)
run(exe, ['-E', exp, '-N', nn, '-Y', nn, pred], false)
abort 'Neighbour input was overwritten' unless File.binread(nn) == before
puts 'nearest-neighbour checks passed'
