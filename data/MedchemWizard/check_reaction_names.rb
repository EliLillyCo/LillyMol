#!/usr/bin/env ruby

require 'optparse'
require 'pathname'

def usage(parser)
  $stderr << parser << "\n"
  exit 1
end

def extract_reaction_name(path)
  contents = path.read

  if (match = contents.match(/^\s*name:\s*"([^"]+)"/))
    return match[1]
  end

  if (match = contents.match(/\(A\s+C\s+Comment\s+"([^"]+)"\)/))
    return match[1]
  end

  nil
end

def expected_name(fname)
  File.basename(fname).sub(/\.(rxn|textproto)\z/, '')
end

def each_reaction(reactions)
  reactions.each_line.with_index(1) do |line, lineno|
    line = line.strip
    next if line.empty? || line.start_with?('#')

    is_proto = line.start_with?('PROTO:')
    fname = is_proto ? line.sub(/\APROTO:/, '') : line
    yield lineno, fname, is_proto
  end
end

options = {
  reactions: nil,
  show_matches: false
}

parser = OptionParser.new do |opts|
  opts.banner = "Usage: #{File.basename($PROGRAM_NAME)} [options]"

  opts.on('-R', '--reactions FILE', 'REACTIONS file to check') do |value|
    options[:reactions] = value
  end

  opts.on('--show-matches', 'Also report matching records') do
    options[:show_matches] = true
  end

  opts.on('-h', '--help', 'Show usage') do
    puts opts
    exit 0
  end
end

parser.parse!

script_dir = Pathname.new(__dir__)
reactions = Pathname.new(options[:reactions] || script_dir.join('REACTIONS'))
usage(parser) unless reactions.file?

root = reactions.dirname
records = 0
missing = 0
no_name = 0
mismatches = 0
matches = 0
name_to_files = Hash.new { |hash, key| hash[key] = [] }

puts "status line file expected_name reaction_name"

each_reaction(reactions) do |lineno, fname, _is_proto|
  records += 1
  path = root.join(fname)

  unless path.file?
    missing += 1
    puts "missing #{lineno} #{fname} #{expected_name(fname)} ."
    next
  end

  name = extract_reaction_name(path)
  unless name
    no_name += 1
    puts "no_name #{lineno} #{fname} #{expected_name(fname)} ."
    next
  end

  name_to_files[name] << fname
  expected = expected_name(fname)

  if expected == name
    matches += 1
    puts "match #{lineno} #{fname} #{expected} #{name}" if options[:show_matches]
  else
    mismatches += 1
    puts "mismatch #{lineno} #{fname} #{expected} #{name}"
  end
end

duplicate_names = name_to_files.select { |_name, files| files.size > 1 }
duplicate_names.each do |name, files|
  puts "duplicate_name . #{files.join(',')} . #{name}"
end

$stderr << "records #{records}\n"
$stderr << "matches #{matches}\n"
$stderr << "mismatches #{mismatches}\n"
$stderr << "missing #{missing}\n"
$stderr << "no_name #{no_name}\n"
$stderr << "duplicate_names #{duplicate_names.size}\n"

exit(missing.zero? && no_name.zero? && mismatches.zero? && duplicate_names.empty? ? 0 : 1)
