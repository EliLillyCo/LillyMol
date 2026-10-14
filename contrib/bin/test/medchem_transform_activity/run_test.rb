#!/usr/bin/env ruby

require 'fileutils'
require 'json'
require 'open3'
require 'tmpdir'

root = File.expand_path('../../../..', __dir__)
script = File.join(root, 'contrib/bin/medchem_transform_activity.rb')
fixture = __dir__

def run_or_die(cmd)
  return if system(*cmd)

  warn "Command failed: #{cmd.join(' ')}"
  exit 1
end

def compare_text(fname, expected_fname, label)
  expected = File.read(expected_fname)
  got = File.read(fname)
  return if got == expected

  warn "#{label} mismatch"
  warn 'Expected:'
  warn expected
  warn 'Got:'
  warn got
  exit 1
end

def compare_json(fname, expected_fname, label)
  expected = JSON.parse(File.read(expected_fname))
  got = JSON.parse(File.read(fname))
  return if got == expected

  warn "#{label} mismatch"
  warn 'Expected:'
  warn JSON.pretty_generate(expected)
  warn 'Got:'
  warn JSON.pretty_generate(got)
  exit 1
end

Dir.mktmpdir('medchem_transform_activity_test') do |tmpdir|
  output = File.join(tmpdir, 'output.txt')
  profile_json = File.join(tmpdir, 'profile.json')
  dbname = File.join(tmpdir, 'structures.bdb')

  run_or_die([
    RbConfig.ruby,
    script,
    '-A', File.join(fixture, 'activity.txt'),
    '-d', dbname,
    '-o', output,
    '--profile-json', profile_json,
    '--reactions', File.join(fixture, 'REACTIONS'),
    '--buildsmidb', File.join(fixture, 'fake_buildsmidb.rb'),
    '--medchem-wizard', File.join(fixture, 'fake_medchem_wizard.rb'),
    '--in-database', File.join(fixture, 'fake_in_database.rb'),
    '--model-command', File.join(fixture, 'fake_model.rb'),
    File.join(fixture, 'input.smi')
  ])

  compare_text(output, File.join(fixture, 'expected.txt'), 'Unified output')
  compare_json(profile_json, File.join(fixture, 'expected_profile.json'), 'Unified profile JSON')
end

Dir.mktmpdir('medchem_transform_activity_split_test') do |tmpdir|
  found = File.join(tmpdir, 'found.smi')
  notfound = File.join(tmpdir, 'notfound.smi')
  predictions = File.join(tmpdir, 'predictions.txt')
  output = File.join(tmpdir, 'output.txt')
  profile_json = File.join(tmpdir, 'profile.json')
  dbname = File.join(tmpdir, 'structures.bdb')

  run_or_die([
    RbConfig.ruby,
    script,
    'generate',
    '-d', dbname,
    '--found', found,
    '--notfound', notfound,
    '--buildsmidb', File.join(fixture, 'fake_buildsmidb.rb'),
    '--medchem-wizard', File.join(fixture, 'fake_medchem_wizard.rb'),
    '--in-database', File.join(fixture, 'fake_in_database.rb'),
    File.join(fixture, 'input.smi')
  ])

  run_or_die([RbConfig.ruby, File.join(fixture, 'fake_model.rb'), notfound, predictions])

  run_or_die([
    RbConfig.ruby,
    script,
    'analyse',
    '-A', File.join(fixture, 'activity.txt'),
    '--found', found,
    '--predictions', predictions,
    '-o', output,
    '--profile-json', profile_json,
    '--reactions', File.join(fixture, 'REACTIONS')
  ])

  compare_text(output, File.join(fixture, 'expected.txt'), 'Split output')
  compare_json(profile_json, File.join(fixture, 'expected_profile.json'), 'Split profile JSON')

  # Restart consumes existing results, even when generation executables are unavailable.
  run_or_die([
    RbConfig.ruby, script, 'restart',
    '-A', File.join(fixture, 'activity.txt'), '--found', found,
    '--predictions', predictions, '-o', output, '--profile-json', profile_json,
    '--reactions', File.join(fixture, 'REACTIONS'),
    '--buildsmidb', File.join(tmpdir, 'unavailable_builder'),
    '--medchem-wizard', File.join(tmpdir, 'unavailable_wizard'),
    '--in-database', File.join(tmpdir, 'unavailable_lookup')
  ])
  compare_text(output, File.join(fixture, 'expected.txt'), 'Restart output')
  compare_json(profile_json, File.join(fixture, 'expected_profile.json'), 'Restart profile JSON')
end

# Reject malformed wizard names before they reach lookup/scoring.
[
  "CO mol1__MCTA__add oxygen\n",
  "CO mol1__MCTA__add__MCTA__oxygen\n",
  "CO mol1__MCTA__\n",
  "CO mol1 add_oxygen\n"
].each do |input|
  _out, _err, status = Open3.capture3(RbConfig.ruby, script,
                                     '--validate-generated-names', stdin_data: input)
  raise "Accepted malformed wizard name: #{input}" if status.success?
end

Dir.mktmpdir('medchem_transform_activity_legacy_test') do |tmpdir|
  found = File.join(tmpdir, 'found.smi')
  predictions = File.join(tmpdir, 'predictions.txt')
  output = File.join(tmpdir, 'output.txt')
  File.write(found, "CC mol1 add_carbon mol2\nCCC mol2 add_carbon mol3\n")
  File.write(predictions, "id transformation prediction\nmol1 add_oxygen 3\nmol2 add_nitrogen 1.5\n")
  run_or_die([RbConfig.ruby, script, 'analyse', '-A', File.join(fixture, 'activity.txt'),
              '--found', found, '--predictions', predictions, '-o', output])
  compare_text(output, File.join(fixture, 'expected.txt'), 'Legacy analysis')
end

# Explicit stems are normalized, including names with an existing other suffix.
Dir.mktmpdir('medchem_transform_activity_suffix_test') do |tmpdir|
  found_stem = File.join(tmpdir, 'found')
  notfound_stem = File.join(tmpdir, 'notfound.txt')
  predictions = File.join(tmpdir, 'predictions.txt')
  output = File.join(tmpdir, 'output.txt')
  run_or_die([
    RbConfig.ruby, script, 'generate',
    '-d', File.join(tmpdir, 'structures.bdb'), '-o', output,
    '--found', found_stem, '--notfound', notfound_stem,
    '--predictions', predictions,
    '--buildsmidb', File.join(fixture, 'fake_buildsmidb.rb'),
    '--medchem-wizard', File.join(fixture, 'fake_medchem_wizard.rb'),
    '--in-database', File.join(fixture, 'fake_in_database.rb'),
    File.join(fixture, 'input.smi')
  ])
  raise 'Found suffix missing' unless File.file?("#{found_stem}.smi")
  raise 'Notfound suffix missing' unless File.file?("#{notfound_stem}.smi")
  run_or_die([RbConfig.ruby, File.join(fixture, 'fake_model.rb'), "#{notfound_stem}.smi", predictions])
  run_or_die([RbConfig.ruby, script, 'analyse', '-A', File.join(fixture, 'activity.txt'),
              '--found', found_stem, '--predictions', predictions, '-o', output])
  compare_text(output, File.join(fixture, 'expected.txt'), 'Normalized filenames')
end

puts 'medchem_transform_activity test passed'
