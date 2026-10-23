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
    '--msort-parallel', File.join(fixture, 'fake_msort_parallel.rb'),
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
    '--msort-parallel', File.join(fixture, 'fake_msort_parallel.rb'),
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
    '--in-database', File.join(tmpdir, 'unavailable_lookup'),
    '--msort-parallel', File.join(tmpdir, 'unavailable_sorter')
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
    '--msort-parallel', File.join(fixture, 'fake_msort_parallel.rb'),
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

# Sorting preserves identifiers and the original lookup output on failure.
Dir.mktmpdir('medchem sorting test ') do |tmpdir|
  wizard = File.join(tmpdir, 'wizard.rb')
  products = (1..10).to_a.reverse.map { |n| "#{'C' * n} mol1__MCTA__sort_#{n}\n" }.join
  File.write(wizard, "#!/usr/bin/env ruby\nprint #{products.dump}\n")
  File.chmod(0o755, wizard)
  found = File.join(tmpdir, 'found.smi')
  notfound = File.join(tmpdir, 'notfound.smi')
  log = File.join(tmpdir, 'sort.log')
  command = [RbConfig.ruby, script, 'generate', '-d', File.join(tmpdir, 'db'),
             '--found', found, '--notfound', notfound,
             '--buildsmidb', File.join(fixture, 'fake_buildsmidb.rb'),
             '--medchem-wizard', wizard,
             '--in-database', File.join(fixture, 'fake_in_database.rb'),
             '--msort-parallel', File.join(fixture, 'fake_msort_parallel.rb')]
  input = File.join(fixture, 'input.smi')
  _out, err, status = Open3.capture3({'MCTA_SORT_LOG' => log}, *command, input)
  raise err unless status.success?
  expected = products.lines.reverse.join
  raise 'Notfound molecules not sorted' unless File.read(notfound) == expected
  raise 'Default sorting thread count is not 8' unless File.readlines(log, chomp: true)[3] == '8'
  _out, err, status = Open3.capture3({'MCTA_SORT_LOG' => log}, *command,
                                    '--sort-threads', '1', input)
  raise err unless status.success?
  raise 'Custom sorting thread count ignored' unless File.readlines(log, chomp: true)[3] == '1'
  _out, _err, status = Open3.capture3({'MCTA_SORT_FAIL' => '1'}, *command, input)
  raise 'Sorter failure ignored' if status.success?
  raise 'Sorter failure damaged lookup output' unless File.read(notfound) == products
  raise 'Sorting temporary file leaked' unless Dir.glob(File.join(tmpdir, 'mcta_sorted_*')).empty?
  %w[MCTA_SORT_TRUNCATE MCTA_SORT_EXTRA].each do |variable|
    _out, err, status = Open3.capture3({variable => '1'}, *command, input)
    raise 'Sorter size mismatch accepted' if status.success?
    raise 'Missing size mismatch diagnostic' unless err.include?('Sorted file size mismatch')
    raise 'Size mismatch damaged lookup output' unless File.read(notfound) == products
    raise 'Size mismatch leaked temporary file' unless Dir.glob(File.join(tmpdir, 'mcta_sorted_*')).empty?
  end
  _out, _err, status = Open3.capture3(*command, '--sort-threads', '0', input)
  raise 'Invalid sort thread count accepted' if status.success?
  # Samples are grouped by reaction, across different starting identifiers.
  products = "CCCCCC mol3__MCTA__add_oxygen\nCCCCC mol3__MCTA__add_nitrogen\n" \
             "CCCC mol2__MCTA__add_oxygen\nCCC mol2__MCTA__add_nitrogen\n" \
             "CC mol1__MCTA__add_oxygen\nC mol1__MCTA__add_nitrogen\n"
  File.write(wizard, "#!/usr/bin/env ruby\nprint #{products.dump}\n")
  sample_command = [*command, '--unique-rows', File.join(fixture, 'fake_unique_rows.rb'),
                    '--max-per-reaction', '2']
  _out, err, status = Open3.capture3(*sample_command, input)
  raise err unless status.success?
  expected = products.lines.reverse.first(4).join
  raise 'Wrong smallest reaction samples' unless File.read(notfound) == expected
  _out, _err, status = Open3.capture3({'MCTA_SAMPLE_FAIL' => '1'}, *sample_command, input)
  raise 'Sampling failure ignored' if status.success?
  raise 'Sampling failure damaged sorted file' unless File.read(notfound) == products.lines.reverse.join
  raise 'Sampling temporary file leaked' unless Dir.glob(File.join(tmpdir, 'mcta_sampled_*')).empty?
  %w[0 -1 invalid].each do |limit|
    _out, _err, status = Open3.capture3(*command, '--max-per-reaction', limit, input)
    raise 'Invalid sample limit accepted' if status.success?
  end
  # Run mode must score the sampled file, not all generated products.
  run_command = command.dup
  run_command[2] = 'run'
  _out, err, status = Open3.capture3(*run_command, '--max-per-reaction', '1',
                                    '--unique-rows', File.join(fixture, 'fake_unique_rows.rb'),
                                    '--model-command', File.join(fixture, 'fake_model.rb'),
                                    '-A', File.join(fixture, 'activity.txt'),
                                    '-o', File.join(tmpdir, 'table.txt'), input)
  raise err unless status.success?
  raise 'Missing scored reaction rows' unless File.readlines(File.join(tmpdir, 'table.txt')).size == 3
  raise 'Run mode used wrong samples' unless File.read(notfound) == products.lines.reverse.first(2).join
  # No candidates: do not invoke the sorter at all.
  File.write(wizard, "#!/usr/bin/env ruby\n")
  File.delete(log)
  _out, err, status = Open3.capture3({'MCTA_SORT_LOG' => log}, *command, '--max-per-reaction', '2',
                                    '--unique-rows', File.join(tmpdir, 'unavailable_sampler'), input)
  raise err unless status.success?
  raise 'Empty notfound file was sorted' if File.exist?(log)
  raise 'Empty lookup output changed' unless File.empty?(notfound)
end

puts 'medchem_transform_activity test passed'
