#!/usr/bin/env ruby
# Sampling is random: check filenames and partition consistency, not record order.
require 'fileutils'
require 'open3'
require 'set'

exe, indir = ARGV
all_ids = File.readlines(File.join(indir, 'input.smi')).map { |line| line.split[1] }.to_set
scenarios = [
  ['default_multiple', 2, nil, false],
  ['offset_multiple', 2, 7, false],
  ['offset_chronological', 2, 7, true],
  ['default_single', 1, nil, false],
  ['offset_single', 1, 7, false],
  ['zero_single', 1, 0, false]
]
scenarios.each do |name, count, start, chronological|
  FileUtils.mkdir_p(name)
  cmd = [exe, '-R', 'train', '-E', 'test', '-N', count.to_s, '-n', '10',
         '-s', '1', '-S', '.dat', '-M', File.join(indir, 'input.smi')]
  cmd += ['-b', start.to_s] unless start.nil?
  cmd << '-C' if chronological
  cmd << File.join(indir, 'activity.txt')
  stdout, stderr, status = Open3.capture3(*cmd, chdir: name)
  abort "#{name}: #{stderr}" unless status.success?
  numbered = count > 1 || !start.nil?
  numbers = (start || 0)...((start || 0) + count)
  expected = []
  numbers.each do |number|
    %w[train test].each do |stem|
      expected << (numbered ? "#{stem}#{number}.dat" : stem)
      expected << "#{stem}#{number}.smi"
    end
  end
  if chronological
    %w[train test].each do |stem|
      expected += ["#{stem}#{(start || 0) + count}.dat", "#{stem}#{(start || 0) + count}.smi"]
    end
  end
  actual = Dir.children(name).sort
  abort "#{name}: expected #{expected.sort}, got #{actual}" unless actual == expected.sort
  numbers.each do |number|
    selections = []
    %w[train test].each do |stem|
      filename = numbered ? "#{stem}#{number}.dat" : stem
      lines = File.readlines(File.join(name, filename), chomp: true)
      abort "#{name}: missing header" unless lines.shift == 'ID activity'
      ids = lines.map { |line| line.split[0] }
      smiles_ids = File.readlines(File.join(name, "#{stem}#{number}.smi")).map { |line| line.split[1] }
      abort "#{name}: activity/SMILES mismatch" unless ids == smiles_ids
      abort "#{name}: wrong selection size" unless ids.size == 10 && ids.to_set.size == 10
      selections << ids.to_set
    end
    abort "#{name}: overlapping train/test" unless (selections[0] & selections[1]).empty?
    abort "#{name}: missing records" unless (selections[0] | selections[1]) == all_ids
  end
  puts "#{name}: passed"
end

# A single split can still send its training records to stdout when -b is set.
FileUtils.mkdir_p('stdout_single')
cmd = [exe, '-E', 'test', '-N', '1', '-b', '7', '-n', '10', '-s', '1',
       '-S', '.dat', File.join(indir, 'activity.txt')]
stdout, stderr, status = Open3.capture3(*cmd, chdir: 'stdout_single')
abort stderr unless status.success?
abort 'stdout_single: unexpected files' unless Dir.children('stdout_single') == ['test7.dat']
lines = stdout.lines.map(&:chomp)
abort 'stdout_single: missing header' unless lines.shift == 'ID activity'
train_ids = lines.map { |line| line.split[0] }.to_set
test_ids = File.readlines('stdout_single/test7.dat').drop(1).map { |line| line.split[0] }.to_set
abort 'stdout_single: wrong selection size' unless train_ids.size == 10 && test_ids.size == 10
abort 'stdout_single: overlapping selections' unless (train_ids & test_ids).empty?
abort 'stdout_single: missing records' unless (train_ids | test_ids) == all_ids
puts 'stdout_single: passed'
